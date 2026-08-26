#include <charconv>
#include <numbers>

#include <poisfft.h>

#include <Igor/Defer.hpp>
#include <Igor/Logging.hpp>
#include <Igor/Math.hpp>
#include <Igor/Timer.hpp>

#include "BoundaryConditions.hpp"
#include "Common.hpp"
#include "Grid.hpp"
#include "IO.hpp"
#include "Mac.hpp"
#include "Monitor.hpp"
#include "VTKWriter.hpp"

using Float              = double;

constexpr Float x_min    = 0.0;
constexpr Float x_max    = 1.0;
constexpr Float y_min    = 0.0;
constexpr Float y_max    = 1.0;

constexpr Float rho      = 1.0;
constexpr Float mu       = 1e-2;

constexpr Float CFL      = 0.5;
constexpr Float tend     = 1.0;
constexpr Float dt_write = tend / 100.0;

// =================================================================================================
auto F(Float t) -> Float {
  constexpr auto pi = std::numbers::pi_v<Float>;
  return std::exp(-2.0 * (mu / rho) * (2.0 * pi) * (2.0 * pi) * t);
}
auto u_analytical(Float x, Float y, Float t) -> Float {
  constexpr auto pi = std::numbers::pi_v<Float>;
  return std::sin(2.0 * pi * x) * std::cos(2.0 * pi * y) * F(t);
}
auto v_analytical(Float x, Float y, Float t) -> Float {
  constexpr auto pi = std::numbers::pi_v<Float>;
  return -std::cos(2.0 * pi * x) * std::sin(2.0 * pi * y) * F(t);
}

// =================================================================================================
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
      .left   = Periodic{},
      .right  = Periodic{},
      .bottom = Periodic{},
      .top    = Periodic{},
  };

  grid.foreach_face_i<Dimension::X>(
      FOREACH_FUNC { u.x(i, j) = u_analytical(grid.x(i), grid.ym(j), 0.0); });
  grid.foreach_face_i<Dimension::Y>(
      FOREACH_FUNC { u.y(i, j) = v_analytical(grid.xm(i), grid.y(j), 0.0); });
  apply_velocity_bconds(grid, bconds, u);
  interpolate(grid, u, ui);

  VTKWriter writer(output_dir, grid);
  writer.add_field("u", ui);
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
      solver.execute(dp.data(), div.data(), ngs.data(), ngs.data());
      apply_neumann_bconds(grid, dp);
      shift_dp_to_zero(grid, dp);

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

  Float L1_u = 0.0;
  grid.foreach_face_i<Dimension::X, Exec::SERIAL>([=, &L1_u](Index i, Index j) {
    const auto u_exp  = u_analytical(grid.x(i), grid.ym(j), t);
    L1_u             += std::abs(u_exp - u.x(i, j)) * grid.dv();
  });
  Float L1_v = 0.0;
  grid.foreach_face_i<Dimension::Y, Exec::SERIAL>([=, &L1_v](Index i, Index j) {
    const auto v_exp  = v_analytical(grid.xm(i), grid.y(j), t);
    L1_v             += std::abs(v_exp - u.y(i, j)) * grid.dv();
  });

  Igor::Info("L1(u) = {:.8e}", L1_u);
  Igor::Info("L1(v) = {:.8e}", L1_v);

  Igor::Info("Ok.");
}
