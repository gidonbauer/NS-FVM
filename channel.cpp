#include <poisfft.h>

#include <Igor/Defer.hpp>
#include <Igor/Logging.hpp>
#include <Igor/Math.hpp>
#include <Igor/Timer.hpp>

#include "Grid.hpp"
#include "IO.hpp"
#include "VTKWriter.hpp"

using Float              = double;

constexpr Float x_min    = 0.0;
constexpr Float x_max    = 1.0;
constexpr Float y_min    = 0.0;
constexpr Float y_max    = 1.0;
constexpr Index N        = 64;

constexpr Float rho      = 1.0;
constexpr Float mu       = 1.0;
constexpr Float Uin      = 1.0;

constexpr Float tend     = 0.1;
constexpr Float dt_write = tend / 100.0;

// constexpr Index NUM_SUB_ITER = 5;

template <typename Float, Layout LAYOUT>
constexpr void calc_flux(const Grid<Float, LAYOUT>& grid,
                         const Vector<Float, LAYOUT> u,
                         FaceVector<Float, LAYOUT> FU,
                         FaceVector<Float, LAYOUT> FV) {
  grid.template foreach_face_i<Dimension::X>(FOREACH_FUNC {
    const auto ui     = (u.x(i, j) + u.x(i - 1, j)) / 2.0;
    const auto vi     = (u.y(i, j) + u.y(i - 1, j)) / 2.0;

    const auto dudx   = (u.x(i, j) - u.x(i - 1, j)) / grid.dx();
    const auto dvdx   = (u.y(i, j) - u.y(i - 1, j)) / grid.dx();

    const auto dudy_l = (u.x(i - 1, j + 1) - u.x(i - 1, j - 1)) / (2.0 * grid.dy());
    const auto dudy_r = (u.x(i, j + 1) - u.x(i, j - 1)) / (2.0 * grid.dy());
    const auto dudy   = (dudy_r + dudy_l) / 2.0;

    FU.x(i, j)        = Igor::sqr(ui) + 2.0 * mu / rho * dudx;
    FV.x(i, j)        = ui * vi + mu / rho * (dudy + dvdx);
  });
  grid.template foreach_face_i<Dimension::Y>(FOREACH_FUNC {
    const auto ui     = (u.x(i, j) + u.x(i, j - 1)) / 2.0;
    const auto vi     = (u.y(i, j) + u.y(i, j - 1)) / 2.0;

    const auto dudy   = (u.x(i, j) - u.x(i, j - 1)) / grid.dy();
    const auto dvdy   = (u.y(i, j) - u.y(i, j - 1)) / grid.dy();

    const auto dvdx_b = (u.y(i + 1, j - 1) - u.y(i - 1, j - 1)) / (2.0 * grid.dx());
    const auto dvdx_t = (u.y(i + 1, j) - u.y(i - 1, j)) / (2.0 * grid.dx());
    const auto dvdx   = (dvdx_t + dvdx_b) / 2.0;

    FU.y(i, j)        = ui * vi + mu / rho * (dudy + dvdx);
    FV.y(i, j)        = Igor::sqr(vi) + 2.0 * mu / rho * dvdy;
  });
}

template <typename Float, Layout LAYOUT>
constexpr void boundary_conditions(Vector<Float, LAYOUT> u) {
  // Left: Dirichlet
  for (Index i = 0; i < u.x.nx(); ++i) {
    for (Index j = -u.x.nghost(); j < 0; ++j) {
      u.x(i, j) = Uin;
      u.y(i, j) = 0.0;
    }
  }

  // Right: homogeneous Neumann
  for (Index i = 0; i < u.x.nx(); ++i) {
    for (Index j = u.x.ny(); j < u.x.ny() + u.x.nghost(); ++j) {
      u.x(i, j) = u.x(i, 2 * u.x.ny() - 1 - j);
      u.y(i, j) = u.y(i, 2 * u.y.ny() - 1 - j);
    }
  }

  // Top & bottom: homogeneous Dirichlet
  for (Index j = 0; j < u.x.ny(); ++j) {
    for (Index i = -u.x.nghost(); i < 0; ++i) {
      u.x(i, j) = 0.0;
      u.y(i, j) = 0.0;
    }
    for (Index i = u.x.nx(); i < u.x.nx() + u.x.nghost(); ++i) {
      u.x(i, j) = 0.0;
      u.y(i, j) = 0.0;
    }
  }
}

template <typename Float, Layout LAYOUT>
constexpr void update_u(const Grid<Float, LAYOUT>& grid,
                        Float dt,
                        const FaceVector<Float, LAYOUT> FU,
                        const FaceVector<Float, LAYOUT> FV,
                        const Vector<Float, LAYOUT> u_old,
                        Vector<Float, LAYOUT> u) {
  grid.foreach_i(FOREACH_FUNC {
    u.x(i, j) = u_old.x(i, j) - dt * ((FU.right(i, j) - FU.left(i, j)) / grid.dx() +
                                      (FU.top(i, j) - FU.bottom(i, j)) / grid.dy());
    u.y(i, j) = u_old.y(i, j) - dt * ((FV.right(i, j) - FV.left(i, j)) / grid.dx() +
                                      (FV.top(i, j) - FV.bottom(i, j)) / grid.dy());
  });
}

template <typename Float, Layout LAYOUT>
void calc_div(const Grid<Float, LAYOUT>& grid, Vector<Float, LAYOUT> u, Scalar<Float, LAYOUT> div) {
  grid.foreach_i(FOREACH_FUNC {
    div(i, j) = (u.x(i + 1, j) - u.x(i - 1, j)) / (2.0 * grid.dx()) +  //
                (u.y(i, j + 1) - u.y(i, j - 1)) / (2.0 * grid.dy());
  });
}

auto main() -> int {
  const auto output_dir = get_output_directory();
  if (!init_output_directory(output_dir)) { return 1; }

  Grid<Float, Layout::C> grid(x_min, x_max, N, y_min, y_max, N, 1);

  auto u_old = grid.alloc_vector();
  IGOR_DEFER(grid.free(u_old););

  auto u = grid.alloc_vector();
  IGOR_DEFER(grid.free(u););

  auto FU = grid.alloc_face_vector();
  IGOR_DEFER(grid.free(FU););

  auto FV = grid.alloc_face_vector();
  IGOR_DEFER(grid.free(FV););

  auto p = grid.alloc_scalar();
  IGOR_DEFER(grid.free(p););

  auto div = grid.alloc_scalar();
  IGOR_DEFER(grid.free(div););

  Float dt = 0.0;
  Float t  = 0.0;

  // = Linear solver ===============================================================================
  const std::array<int, 2> ns   = {grid.nx() + 2 * grid.nghost(), grid.ny() + 2 * grid.nghost()};
  const std::array<Float, 2> Ls = {(grid.x_max() - grid.x_min()) + 2.0 * grid.dx(),
                                   (grid.y_max() - grid.y_min()) + 2.0 * grid.dy()};
  const std::array<int, 4> BCs  = {
      PoisFFT::NEUMANN, PoisFFT::NEUMANN, PoisFFT::NEUMANN, PoisFFT::NEUMANN};
  PoisFFT::Solver<2, Float> solver(ns.data(), Ls.data(), BCs.data(), PoisFFT::SPECTRAL);
  // = Linear solver ===============================================================================

  grid.foreach_i(FOREACH_FUNC {
    u.x(i, j) = Uin;
    u.y(i, j) = 0.0;
  });
  boundary_conditions(u);

  VTKWriter writer(output_dir, grid);
  writer.add_field("u", u);
  writer.add_field("p", p);
  if (!writer.write(t)) { return 1; }

  IGOR_TIME_SCOPE("Solver")
  while (t < tend) {
    Float u_max = 0.0;
    grid.foreach_i<Exec::SERIAL>([=, &u_max](Index i, Index j) {
      u_max = std::max({std::abs(u.x(i, j)), std::abs(u.y(i, j)), u_max});
    });
    dt = std::min(0.5 * std::min(grid.dx(), grid.dy()) / u_max, dt_write);
    dt = std::min(dt, tend - t);

    copy(u, u_old);

#if 0
    for (Index sub_iter = 0; sub_iter < NUM_SUB_ITER; ++sub_iter) {
      grid.foreach_a(FOREACH_FUNC { u.x(i, j) = 0.5 * (u.x(i, j) + u_old.x(i, j)); });

      calc_reconstruction(grid, u, uL, uR);
      calc_flux(grid, uL, uR, F);
      update_u(grid, dt, F, u_old, u);
      boundary_conditions(u);
    }
#else
    calc_flux(grid, u, FU, FV);
    update_u(grid, dt / 2.0, FU, FV, u_old, u);
    boundary_conditions(u);

    calc_div(grid, u, div);
    grid.foreach_i(FOREACH_FUNC { div(i, j) /= dt; });
    solver.execute(p.data(), div.data());
    grid.foreach_i(FOREACH_FUNC {
      u.x(i, j) -= dt * (p(i + 1, j) - p(i - 1, j)) / (2.0 * grid.dx());
      u.y(i, j) -= dt * (p(i, j + 1) - p(i, j - 1)) / (2.0 * grid.dy());
    });

    calc_flux(grid, u, FU, FV);
    update_u(grid, dt, FU, FV, u_old, u);
    boundary_conditions(u);

    calc_div(grid, u, div);
    grid.foreach_i(FOREACH_FUNC { div(i, j) /= dt; });
    solver.execute(p.data(), div.data());
    grid.foreach_i(FOREACH_FUNC {
      u.x(i, j) -= dt * (p(i + 1, j) - p(i - 1, j)) / (2.0 * grid.dx());
      u.y(i, j) -= dt * (p(i, j + 1) - p(i, j - 1)) / (2.0 * grid.dy());
    });
#endif

    t += dt;
    if (should_save(t, dt, dt_write, tend)) {
      if (!writer.write(t)) { return 1; }
    }
  }

  Igor::Info("Ok.");
}
