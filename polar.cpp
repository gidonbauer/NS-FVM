#include <charconv>
#include <numbers>

#include <poisfft.h>

#include <Igor/Timer.hpp>

#include "BoundaryConditions.hpp"
#include "Common.hpp"
#include "Grid.hpp"
#include "IO.hpp"
#include "Mac.hpp"
#include "Monitor.hpp"
#include "MultigridPoisson.hpp"
#include "VTKWriter.hpp"

// = Setup =========================================================================================
using Float               = double;

constexpr Float theta_min = 0.0;
constexpr Float theta_max = 2.0 * std::numbers::pi_v<Float>;
constexpr Float r_min     = 1.0;
constexpr Float r_max     = 2.0;

constexpr Float rho       = 1.0;
constexpr Float mu        = 1.0;
constexpr Float CFL       = 0.7;
constexpr Float tend      = 3.0;

constexpr Float dt_write  = tend / 100.0;

constexpr Float U_wall    = 1.0;
// = Setup =========================================================================================

// =================================================================================================
template <typename Float>
constexpr auto uth_analytical(Float r) -> Float {
  constexpr auto k1 = -1.0 / 3.0;
  constexpr auto k2 = 4.0 / 3.0;
  return k1 * r + k2 / r;
}
// =================================================================================================

auto main(int argc, char** argv) -> int {
  const auto usage_str = Igor::detail::format("Usage: {} <grid size>", argv[0]);
  if (argc < 2) {
    Igor::Error("{}", usage_str);
    return 1;
  }

  Index NY = 0;
  if (std::from_chars(argv[1], argv[1] + std::strlen(argv[1]), NY).ec != std::errc{} || NY <= 0) {
    Igor::Error("{}", usage_str);
    Igor::Error("  Invalid grid size `{}`", argv[1]);
    return 1;
  }
  Index NX              = NY;

  const auto output_dir = get_output_directory();
  if (!init_output_directory(output_dir)) { return 1; }

  Grid<Float, Layout::C> grid(theta_min, theta_max, NX, r_min, r_max, NY, 1, Coordinates::POLAR);

  auto u_old      = grid.alloc_face_vector();
  auto u          = grid.alloc_face_vector();

  auto FUX        = grid.alloc_scalar();
  auto FUY        = grid.alloc_vertex_scalar();
  auto FVX        = grid.alloc_vertex_scalar();
  auto FVY        = grid.alloc_scalar();

  auto ui         = grid.alloc_vector();
  auto p          = grid.alloc_scalar();
  auto dp         = grid.alloc_scalar();
  auto div        = grid.alloc_scalar();

  auto ua         = grid.alloc_vector();

  Float dt        = 0.0;
  Float t         = 0.0;

  Index mg_cycles = 0;
  Float mg_res    = 0.0;

  // = Linear solver ===============================================================================
  // const std::array<int, 2> ns   = {grid.nx(), grid.ny()};
  // const std::array<Float, 2> Ls = {grid.x_max() - grid.x_min(), grid.y_max() - grid.y_min()};
  // const std::array<int, 4> BCs  = {
  //     PoisFFT::NEUMANN_STAG, PoisFFT::NEUMANN_STAG, PoisFFT::NEUMANN_STAG,
  //     PoisFFT::NEUMANN_STAG};
  // PoisFFT::Solver<2, Float> solver(ns.data(), Ls.data(), BCs.data(),
  // PoisFFT::FINITE_DIFFERENCE_2); const std::array<int, 2> ngs = {grid.nghost(), grid.nghost()};
  MultigridSolver solver(grid);
  // = Linear solver ===============================================================================

  const VelocityBConds<Float> bconds{
      .left   = Periodic{},
      .right  = Periodic{},
      .bottom = Dirichlet<Float>{.U = U_wall, .V = 0.0},
      .top    = Dirichlet<Float>{.U = 0.0, .V = 0.0},
  };

  // Initial condition: the analytic fully-developed profile.
  grid.foreach_face_i<Dimension::X>(FOREACH_FUNC { u.x(i, j) = 0.0; });
  grid.foreach_face_i<Dimension::Y>(FOREACH_FUNC { u.y(i, j) = 0.0; });
  apply_velocity_bconds(grid, bconds, u);
  interpolate(grid, u, ui);

  grid.foreach_a(FOREACH_FUNC {
    ua.x(i, j) = uth_analytical(grid.ym(j));
    ua.y(i, j) = 0.0;
  });

  VTKWriter writer(output_dir, grid);
  writer.add_field("u", ui);
  writer.add_field("ua", ua);
  writer.add_field("p", p);
  writer.add_field("div", div);
  if (!writer.write(t)) { return 1; }

  Stats p_stats   = stats(grid, p);
  Stats u_stats   = stats(grid, u.x);
  Stats v_stats   = stats(grid, u.y);
  Stats div_stats = stats(grid, div);

  Float p_max     = std::max(std::abs(p_stats.min), std::abs(p_stats.max));
  Float u_max     = std::max(std::abs(u_stats.min), std::abs(u_stats.max));
  Float v_max     = std::max(std::abs(v_stats.min), std::abs(v_stats.max));
  Float div_max   = std::max(std::abs(div_stats.min), std::abs(div_stats.max));

  Monitor<Float> monitor(output_dir + "/monitor.log");
  monitor.add_variable(&t, "t");
  monitor.add_variable(&dt, "dt");
  monitor.add_variable(&p_max, "absmax(p)");
  monitor.add_variable(&u_max, "absmax(u)");
  monitor.add_variable(&v_max, "absmax(v)");
  monitor.add_variable(&div_max, "absmax(div)");
  monitor.add_variable(&mg_res, "res(MG)");
  monitor.add_variable(&mg_cycles, "cycles(MG)");
  monitor.write();

  IGOR_TIME_SCOPE("Solver")
  while (t < tend) {
    // Time-step size: convective CFL and the (2D) explicit-diffusion limit dt <= h^2 / (4 nu).
    dt = adjust_dt(grid, u, rho, mu, CFL);
    dt = std::min({dt, dt_write, tend - t});

    copy(u, u_old);

    for (Index sub_iter = 0; sub_iter < 2; ++sub_iter) {
      const auto local_dt = sub_iter == 0 ? dt / 2.0 : dt;

      // 1) Predictor
      calc_flux(grid, u, p, rho, mu, FUX, FUY, FVX, FVY);
      update_u(grid, local_dt, FUX, FUY, FVX, FVY, u_old, u);
      apply_velocity_bconds(grid, bconds, u);

      // 2) Pressure correction
      calc_div(grid, u, div);
      grid.foreach_i(FOREACH_FUNC { div(i, j) *= rho / local_dt; });
      // TODO: This is not correct for polar coordinates
      // solver.execute(dp.data(), div.data(), ngs.data(), ngs.data());
      if (!solver.solve(dp, div, 1e-3)) {
        Igor::Warn("Multigrid solver did not converge after {} cycles: res = {:.8e}",
                   solver.num_cycles(),
                   solver.res());
      }
      mg_cycles = solver.num_cycles();
      mg_res    = solver.res();
      apply_neumann_bconds(grid, dp);

      // 3) Project
      correct_velocity(grid, dp, rho, local_dt, u, p);
    }

    interpolate(grid, u, ui);
    calc_div(grid, u, div);

    p_stats    = stats(grid, p);
    u_stats    = stats(grid, u.x);
    v_stats    = stats(grid, u.y);
    div_stats  = stats(grid, div);

    p_max      = std::max(std::abs(p_stats.min), std::abs(p_stats.max));
    u_max      = std::max(std::abs(u_stats.min), std::abs(u_stats.max));
    v_max      = std::max(std::abs(v_stats.min), std::abs(v_stats.max));
    div_max    = std::max(std::abs(div_stats.min), std::abs(div_stats.max));

    t         += dt;
    if (should_save(t, dt, dt_write, tend)) {
      if (!writer.write(t)) { return 1; }
    }
    monitor.write();
  }

  Float L1 = 0.0;
  grid.foreach_range<Exec::SERIAL>(
      u.x.nx() / 2, u.x.nx() / 2 + 1, 0, u.x.ny(), [=, &L1](Index i, Index j) {
        const auto uth_exp  = uth_analytical(grid.ym(j));
        L1                 += std::abs(uth_exp - u.x(i, j)) * grid.dy();
      });
  Igor::Info("NY = {}", NY);
  Igor::Info("L1-error = {:.8e}", L1);

  Igor::Info("Ok.");
}
