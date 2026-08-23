#include <charconv>

#include <poisfft.h>

#include <Igor/Defer.hpp>
#include <Igor/Logging.hpp>
#include <Igor/Math.hpp>
#include <Igor/Timer.hpp>

#include "BoundaryConditions.hpp"
#include "Grid.hpp"
#include "IO.hpp"
#include "VTKWriter.hpp"

using Float              = double;

constexpr Float x_min    = 0.0;
constexpr Float x_max    = 1.0;
constexpr Float y_min    = 0.0;
constexpr Float y_max    = 1.0;

constexpr Float rho      = 1.0;
constexpr Float mu       = 1.0;
constexpr Float Uin      = 1.0;

constexpr Float tend     = 5e-2;
constexpr Float dt_write = tend / 100.0;

// Fully-developed (Poiseuille) inflow profile: parabolic, zero at the walls, peak value Uin at the
// centreline. A wall-consistent inlet avoids the corner singularity a top-hat profile would create.
constexpr auto inlet_u(Float y) -> Float {
  constexpr Float H = y_max - y_min;
  const Float s     = (y - y_min) / H;  // Normalised wall distance in [0, 1].
  return Uin * 4.0 * s * (1.0 - s);
}

// =================================================================================================
template <typename Float, Layout LAYOUT>
constexpr void calc_div(const Grid<Float, LAYOUT>& grid,
                        const FaceVector<Float, LAYOUT> uf,
                        Scalar<Float, LAYOUT> div) {
  grid.foreach_i(FOREACH_FUNC {
    div(i, j) = (uf.right(i, j) - uf.left(i, j)) / grid.dx() +  //
                (uf.top(i, j) - uf.bottom(i, j)) / grid.dy();
  });
}

// =================================================================================================
template <typename Float, Layout LAYOUT>
constexpr void interpolate(const Grid<Float, LAYOUT>& grid,
                           const FaceVector<Float, LAYOUT> uf,
                           Vector<Float, LAYOUT> ui) {
  grid.foreach_i(FOREACH_FUNC {
    ui.x(i, j) = (uf.right(i, j) + uf.left(i, j)) / 2.0;
    ui.y(i, j) = (uf.top(i, j) + uf.bottom(i, j)) / 2.0;
  });
}

// =================================================================================================
template <typename Float, Layout LAYOUT>
constexpr void calc_flux(const Grid<Float, LAYOUT>& grid,
                         const FaceVector<Float, LAYOUT> u,
                         const Scalar<Float, LAYOUT> p,
                         Scalar<Float, LAYOUT> FUX,
                         VertexScalar<Float, LAYOUT> FUY,
                         VertexScalar<Float, LAYOUT> FVX,
                         Scalar<Float, LAYOUT> FVY) {
  std::fill_n(FUX.data(), FUX.size(), std::numeric_limits<Float>::quiet_NaN());
  std::fill_n(FUY.data(), FUY.size(), std::numeric_limits<Float>::quiet_NaN());
  std::fill_n(FVX.data(), FVX.size(), std::numeric_limits<Float>::quiet_NaN());
  std::fill_n(FVY.data(), FVY.size(), std::numeric_limits<Float>::quiet_NaN());

  grid.foreach_a(FOREACH_FUNC {
    const auto ui   = (u.right(i, j) + u.left(i, j)) / 2.0;
    const auto dudx = (u.right(i, j) - u.left(i, j)) / grid.dx();
    FUX(i, j)       = -Igor::sqr(ui) - p(i, j) / rho + 2.0 * mu / rho * dudx;

    const auto vi   = (u.top(i, j) + u.bottom(i, j)) / 2.0;
    const auto dvdy = (u.top(i, j) - u.bottom(i, j)) / grid.dy();
    FVY(i, j)       = -Igor::sqr(vi) - p(i, j) / rho + 2.0 * mu / rho * dvdy;
    ;
  });

  grid.foreach_vertex_i(FOREACH_FUNC {
    const auto ui   = (u.x(i, j) + u.x(i, j - 1)) / 2.0;
    const auto dudy = (u.x(i, j) - u.x(i, j - 1)) / grid.dy();

    const auto vi   = (u.y(i, j) + u.y(i - 1, j)) / 2.0;
    const auto dvdx = (u.y(i, j) - u.y(i - 1, j)) / grid.dx();

    FUY(i, j)       = -ui * vi + mu / rho * (dudy + dvdx);
    FVX(i, j)       = -ui * vi + mu / rho * (dudy + dvdx);
  });
}

// =================================================================================================
template <typename Float, Layout LAYOUT>
constexpr void update_u(const Grid<Float, LAYOUT>& grid,
                        Float dt,
                        const Scalar<Float, LAYOUT> FUX,
                        const VertexScalar<Float, LAYOUT> FUY,
                        const VertexScalar<Float, LAYOUT> FVX,
                        const Scalar<Float, LAYOUT> FVY,
                        const FaceVector<Float, LAYOUT> u_old,
                        FaceVector<Float, LAYOUT> u) {
  grid.template foreach_face_i<Dimension::X>(FOREACH_FUNC {
    u.x(i, j) = u_old.x(i, j) + dt * ((FUX(i, j) - FUX(i - 1, j)) / grid.dx() +
                                      (FUY(i, j + 1) - FUY(i, j)) / grid.dy());
  });

  grid.template foreach_face_i<Dimension::Y>(FOREACH_FUNC {
    u.y(i, j) = u_old.y(i, j) + dt * ((FVX(i + 1, j) - FVX(i, j)) / grid.dx() +
                                      (FVY(i, j) - FVY(i, j - 1)) / grid.dy());
  });
}

// =================================================================================================
template <typename Float, Layout LAYOUT>
void correct_outflow(const Grid<Float, LAYOUT>& grid, FaceVector<Float, LAYOUT> u) {
  IGOR_ASSERT(grid.nghost() == 1, "Expected exactly one ghost cell but got {}", grid.nghost());

  // Rescale the outflow so that it matches the inflow exactly (global continuity).
  Float Qin  = 0.0;
  Float Qout = 0.0;
  for (Index j = 0; j < u.x.ny(); ++j) {
    Qin  += u.x(0, j) * grid.dy();
    Qout += u.x(u.x.nx() - 1, j) * grid.dy();
  }
  const Float corr = (Qin - Qout) / (static_cast<Float>(u.x.ny()) * grid.dy());
  for (Index j = 0; j < u.x.ny(); ++j) {
    u.x(u.x.nx() - 1, j) += corr;
  }
}

// =================================================================================================
template <typename Float, Layout LAYOUT>
constexpr void shift_dp_to_zero(const Grid<Float, LAYOUT>& grid, Scalar<Float, LAYOUT> dp) {
  Float avg_dp = 0.0;
  grid.template foreach_a<Exec::SERIAL>([=, &avg_dp](Index i, Index j) { avg_dp += dp(i, j); });
  avg_dp /= static_cast<Float>(grid.nx() * grid.ny());
  grid.foreach_a(FOREACH_FUNC { dp(i, j) -= avg_dp; });
}

auto main(int argc, char** argv) -> int {
  const auto usage_str = Igor::detail::format("Usage: {} <grid size>", argv[0]);
  if (argc < 2) {
    Igor::Error("{}", usage_str);
    return 1;
  }

  Index N = 0;
  if (std::from_chars(argv[1], argv[1] + std::strlen(argv[1]), N).ec != std::errc{} || N <= 0) {
    Igor::Error("{}", usage_str);
    Igor::Error("  Invalid grid size `{}`", argv[1]);
    return 1;
  }

  const auto output_dir = get_output_directory();
  if (!init_output_directory(output_dir)) { return 1; }

  Grid<Float, Layout::C> grid(x_min, x_max, N, y_min, y_max, N, 1);

  auto u_old = grid.alloc_face_vector();
  auto u     = grid.alloc_face_vector();

  auto FUX   = grid.alloc_scalar();
  auto FUY   = grid.alloc_vertex_scalar();
  auto FVX   = grid.alloc_vertex_scalar();
  auto FVY   = grid.alloc_scalar();

  auto ui    = grid.alloc_vector();
  auto p     = grid.alloc_scalar();  // Pressure (accumulated across steps).
  auto dp    = grid.alloc_scalar();  // Pressure correction of the current step.
  auto div   = grid.alloc_scalar();

  Float dt   = 0.0;
  Float t    = 0.0;

  // = Linear solver ===============================================================================
  const std::array<int, 2> ns   = {grid.nx(), grid.ny()};
  const std::array<Float, 2> Ls = {grid.x_max() - grid.x_min(), grid.y_max() - grid.y_min()};
  const std::array<int, 4> BCs  = {
      PoisFFT::NEUMANN_STAG, PoisFFT::NEUMANN_STAG, PoisFFT::NEUMANN_STAG, PoisFFT::NEUMANN_STAG};
  PoisFFT::Solver<2, Float> solver(ns.data(), Ls.data(), BCs.data(), PoisFFT::FINITE_DIFFERENCE_2);
  const std::array<int, 2> ngs = {grid.nghost(), grid.nghost()};
  // = Linear solver ===============================================================================

  const VelocityBConds<Float> bconds{
      .left   = Dirichlet<Float>{.U = [](Float y, Float /*t*/) { return inlet_u(y); }, .V = 0.0},
      .right  = Neumann{.clipped = true},
      .bottom = Dirichlet<Float>{.U = 0.0, .V = 0.0},
      .top    = Dirichlet<Float>{.U = 0.0, .V = 0.0},
  };

  // Initial condition: the analytic fully-developed profile.
  grid.foreach_face_i<Dimension::X>(FOREACH_FUNC {
    u.x(i, j) = 0.0;  // inlet_u(grid.ym(j));
  });
  grid.foreach_face_i<Dimension::Y>(FOREACH_FUNC { u.y(i, j) = 0.0; });
  apply_velocity_bconds(grid, bconds, u);
  interpolate(grid, u, ui);

  VTKWriter writer(output_dir, grid);
  writer.add_field("u", ui);
  writer.add_field("p", p);
  writer.add_field("div", div);
  if (!writer.write(t)) { return 1; }

  IGOR_TIME_SCOPE("Solver")
  while (t < tend) {
    // Time-step size: convective CFL and the (2D) explicit-diffusion limit dt <= h^2 / (4 nu).
    Float u_max = 0.0;
    grid.foreach_face_i<Dimension::X, Exec::SERIAL>(
        [=, &u_max](Index i, Index j) { u_max = std::max(std::abs(u.x(i, j)), u_max); });
    grid.foreach_face_i<Dimension::Y, Exec::SERIAL>(
        [=, &u_max](Index i, Index j) { u_max = std::max(std::abs(u.y(i, j)), u_max); });
    const auto h = std::min(grid.dx(), grid.dy());
    dt           = std::min({
        0.5 * h / u_max,
        0.2 * Igor::sqr(h) * rho / mu,
        dt_write,
        tend - t,
    });

    copy(u, u_old);

    // 1) Predictor: advance convection + diffusion, pressure
    calc_flux(grid, u, p, FUX, FUY, FVX, FVY);
    update_u(grid, dt, FUX, FUY, FVX, FVY, u_old, u);
    apply_velocity_bconds(grid, bconds, u);
    correct_outflow(grid, u);

    // 2) Pressure correction: solve  laplace(dp) = (rho/dt) div(uf).
    calc_div(grid, u, div);
    grid.foreach_i(FOREACH_FUNC { div(i, j) *= rho / dt; });
    solver.execute(dp.data(), div.data(), ngs.data(), ngs.data());
    apply_neumann_bconds(grid, dp);
    shift_dp_to_zero(grid, dp);

    // 3) Project: correct the face velocity (compact gradient) and the cell velocity (central
    //    gradient of the updated pressure), and accumulate the pressure.
    grid.foreach_a(FOREACH_FUNC { p(i, j) += dp(i, j); });
    grid.foreach_face_i<Dimension::X>(
        FOREACH_FUNC { u.x(i, j) -= (dt / rho) * (dp(i, j) - dp(i - 1, j)) / grid.dx(); });
    grid.foreach_face_i<Dimension::Y>(
        FOREACH_FUNC { u.y(i, j) -= (dt / rho) * (dp(i, j) - dp(i, j - 1)) / grid.dy(); });

    interpolate(grid, u, ui);
    calc_div(grid, u, div);

    t += dt;
    if (should_save(t, dt, dt_write, tend)) {
      if (!writer.write(t)) { return 1; }
    }
  }

  Igor::Info("Ok.");
}
