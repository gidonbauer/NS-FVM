#include <charconv>
#include <numbers>

#include <Igor/Timer.hpp>

#include "Advection-Diffusion.hpp"
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
constexpr Float pi        = std::numbers::pi_v<Float>;

constexpr Float theta_min = 0.0;
constexpr Float theta_max = pi;
constexpr Float r_min     = 1.0;
constexpr Float r_max     = 2.0;

constexpr Float Uavg      = 1.0;
constexpr Float D         = 1e-1;

constexpr Float rho       = 1.0;
constexpr Float mu        = 1.0;
constexpr Float CFL       = 0.7;
constexpr Float tend      = 0.25;

constexpr Float dt_write  = tend / 100.0;
// = Setup =========================================================================================

// =================================================================================================
template <typename Float>
constexpr auto uth_analytical(Float r) -> Float {
  return Uavg / 0.328189 *
         (4.0 * std::numbers::ln2_v<Float> * (r - 1.0 / r) - 3.0 * r * std::log(r));
}
// =================================================================================================

// =================================================================================================
template <typename Float, Layout LAYOUT>
void correct_outflow(const Grid<Float, LAYOUT>& grid, FaceVector<Float, LAYOUT> u) {
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

  Grid<Float> grid(theta_min, theta_max, 2 * N, r_min, r_max, N, 3, Coordinates::POLAR);

  auto u_old      = grid.alloc_face_vector();
  auto u          = grid.alloc_face_vector();
  auto ui         = grid.alloc_vector();

  auto FUX        = grid.alloc_scalar();
  auto FUY        = grid.alloc_vertex_scalar();
  auto FVX        = grid.alloc_vertex_scalar();
  auto FVY        = grid.alloc_scalar();

  auto div        = grid.alloc_scalar();
  auto p          = grid.alloc_scalar();
  auto dp         = grid.alloc_scalar();

  auto T_old      = grid.alloc_scalar();
  auto T          = grid.alloc_scalar();
  auto FT         = grid.alloc_face_vector();

  Float dt        = 0.0;
  Float t         = 0.0;

  Index mg_cycles = 0;
  Float mg_res    = 0.0;

  const BConds<Float> u_bconds{
      .left   = Dirichlet<Float>{.val = [](Float r, Float /*t*/) { return uth_analytical(r); }},
      .right  = Neumann{},
      .bottom = Dirichlet<Float>{.val = 0.0},
      .top    = Dirichlet<Float>{.val = 0.0},
  };

  const BConds<Float> v_bconds{
      .left   = Dirichlet<Float>{.val = 0.0},
      .right  = Neumann{},
      .bottom = Dirichlet<Float>{.val = 0.0},
      .top    = Dirichlet<Float>{.val = 0.0},
  };

  const BConds<Float> T_bconds{
      .left   = Neumann{},
      .right  = Neumann{},
      .bottom = Dirichlet<Float>{.val = [](Float theta,
                                           Float /*t*/) { return 10.0 * (theta < pi / 2.0); }},
      .top    = Dirichlet<Float>{.val = [](Float theta,
                                           Float /*t*/) { return 10.0 * (theta < pi / 2.0); }},
  };

  const BConds<Float> dp_bconds{
      .left   = Neumann{},
      .right  = Neumann{},
      .bottom = Neumann{},
      .top    = Neumann{},
  };

  MultigridSolver solver(grid, dp_bconds);

  // Initial condition: the analytic fully-developed profile.
  grid.foreach_face_i<Dimension::X>(FOREACH_FUNC { u.x(i, j) = Uavg; });
  grid.foreach_face_i<Dimension::Y>(FOREACH_FUNC { u.y(i, j) = 0.0; });
  apply_velocity_bconds(grid, u_bconds, v_bconds, u);
  interpolate(grid, u, ui);

  VTKWriter writer(output_dir, grid);
  writer.add_field("u", ui);
  writer.add_field("p", p);
  writer.add_field("T", T);
  writer.add_field("div", div);
  if (!writer.write(t)) { return 1; }

  Stats p_stats   = stats(grid, p);
  Stats u_stats   = stats(grid, u.x);
  Stats v_stats   = stats(grid, u.y);
  Stats div_stats = stats(grid, div);
  Stats T_stats   = stats(grid, T);
  Float div_max   = std::max(std::abs(div_stats.min), std::abs(div_stats.max));

  Monitor<Float> monitor(output_dir + "/monitor.log");
  monitor.add_variable(&t, "t");
  monitor.add_variable(&dt, "dt");
  monitor.add_variable(&p_stats.max, "abs(p)");
  monitor.add_variable(&u_stats.max, "abs(u)");
  monitor.add_variable(&v_stats.max, "abs(v)");
  monitor.add_variable(&T_stats.min, "min(T)");
  monitor.add_variable(&T_stats.max, "max(T)");
  monitor.add_variable(&div_max, "absmax(div)");
  monitor.add_variable(&mg_res, "res(MG)");
  monitor.add_variable(&mg_cycles, "cycles(MG)");
  monitor.write();

  IGOR_TIME_SCOPE("Solver")
  while (t < tend) {
    dt = std::min({
        adjust_dt(grid, u, rho, mu, CFL),
        advection_adjust_dt(grid, D, CFL),
        dt_write,
        tend - t,
    });

    copy(u, u_old);
    copy(T, T_old);

    for (Index sub_iter = 0; sub_iter < 2; ++sub_iter) {
      const auto local_dt = sub_iter == 0 ? dt / 2.0 : dt;

      // 1) Predictor
      calc_flux(grid, u, p, rho, mu, FUX, FUY, FVX, FVY);
      update_u(grid, local_dt, FUX, FUY, FVX, FVY, u_old, u);
      apply_velocity_bconds(grid, u_bconds, v_bconds, u);
      correct_outflow(grid, u);

      // 2) Pressure correction
      calc_div(grid, u, div);
      grid.foreach_i(FOREACH_FUNC { div(i, j) *= rho / local_dt; });
      if (!solver.solve(dp, div, 1e-4)) {
        Igor::Warn("Multigrid solver did not converge after {} cycles: res = {:.8e}",
                   solver.num_cycles(),
                   solver.res());
      }
      mg_cycles = solver.num_cycles();
      mg_res    = solver.res();
      apply_bconds(grid, dp_bconds, dp, t);

      // 3) Project
      correct_velocity(grid, dp, rho, local_dt, u, p);

      // Update temperature
      advection_calc_flux(grid, u, T, D, FT);
      advection_update_s(grid, local_dt, FT, T_old, T);
      apply_bconds(grid, T_bconds, T, t);
    }

    interpolate(grid, u, ui);
    calc_div(grid, u, div);

    p_stats    = stats(grid, p);
    u_stats    = stats(grid, u.x);
    v_stats    = stats(grid, u.y);
    T_stats    = stats(grid, T);
    div_stats  = stats(grid, div);
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
  Igor::Info("L1({}) = {:.8e}", N, L1);

  Igor::Info("Ok.");
}
