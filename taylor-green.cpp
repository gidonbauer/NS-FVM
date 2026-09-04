#include <charconv>
#include <numbers>

#include <poisfft.h>

#include <Igor/Defer.hpp>
#include <Igor/Logging.hpp>
#include <Igor/Math.hpp>
#include <Igor/Timer.hpp>

#include "Advection-Diffusion.hpp"
#include "BoundaryConditions.hpp"
#include "Common.hpp"
#include "Grid.hpp"
#include "IO.hpp"
#include "Mac.hpp"
#include "Monitor.hpp"
#include "MultigridPoisson.hpp"
#include "Temperature.hpp"
#include "VTKWriter.hpp"

using Float           = double;
constexpr auto pi     = std::numbers::pi_v<Float>;

constexpr Float x_min = 0.0;
constexpr Float x_max = 1.0;
constexpr Float y_min = 0.0;
constexpr Float y_max = 1.0;

// Properties of water @ 20°C and 1 atm
// constexpr Float rho   = 998.21;     // [kg/m^3]
// constexpr Float mu    = 1.0014e-4;  // [Pa*s]
// constexpr Float kappa = 0.59803;    // [W/(m*K)]
// constexpr Float cV    = 4.1566e3;   // [J/(kg*K)]

// Properties of nitrogen @ 20°C and 1 atm
constexpr Float rho      = 1.1648;      // [kg/m^3]
constexpr Float mu       = 1.7573e-05;  // [Pa*s]
constexpr Float kappa    = 0.025473;    // [W/(m*K)]
constexpr Float cV       = 0.74307e3;   // [J/(kg*K)]

constexpr Float CFL      = 0.5;
constexpr Float tend     = 1e3;
constexpr Float dt_write = tend / 100.0;

// =================================================================================================
auto F(Float t) -> Float { return std::exp(-2.0 * (mu / rho) * (2.0 * pi) * (2.0 * pi) * t); }
auto u_analytical(Float x, Float y, Float t) -> Float {
  return 1e-2 * std::sin(2.0 * pi * x) * std::cos(2.0 * pi * y) * F(t);
}
auto v_analytical(Float x, Float y, Float t) -> Float {
  return 1e-2 * -std::cos(2.0 * pi * x) * std::sin(2.0 * pi * y) * F(t);
}

// =================================================================================================
auto main(int argc, char** argv) -> int {
  const auto usage_str =
      Igor::detail::format("Usage: {} [--linear-solver=MG|FFT] <grid size>", argv[0]);

  Index N        = 0;
  bool multigrid = false;
  for (int i = 1; i < argc; ++i) {
    std::string_view arg(argv[i]);
    if (arg.starts_with("--linear-solver")) {
      arg.remove_prefix(std::strlen("--linear-solver"));
      if (arg.starts_with('=')) {
        arg.remove_prefix(1);
      } else {
        if (i + 1 >= argc) {
          Igor::Error("Expected argument for `{}`", argv[i]);
          return 1;
        }
        i   += 1;
        arg  = argv[i];
      }
      if (arg == "MG") {
        multigrid = true;
      } else if (arg == "FFT") {
        multigrid = false;
      } else {
        Igor::Error("{}", usage_str);
        Igor::Error("  Invalid linear solver `{}`", arg);
        return 1;
      }
    } else {
      if (std::from_chars(arg.data(), arg.data() + arg.size(), N).ec != std::errc{} || N <= 0) {
        Igor::Error("{}", usage_str);
        Igor::Error("  Invalid grid size `{}`", arg);
        return 1;
      }
    }
  }

  if (N == 0) {
    Igor::Error("{}", usage_str);
    Igor::Error("  Missing grid size.");
    return 1;
  }

  const auto output_dir = get_output_directory();
  if (!init_output_directory(output_dir)) { return 1; }

  Grid<Float, Layout::C> grid(x_min, x_max, N, y_min, y_max, N, 3);

  auto u_old        = grid.alloc_face_vector();
  auto u            = grid.alloc_face_vector();

  auto FUX          = grid.alloc_scalar();
  auto FUY          = grid.alloc_vertex_scalar();
  auto FVX          = grid.alloc_vertex_scalar();
  auto FVY          = grid.alloc_scalar();

  auto ui           = grid.alloc_vector();
  auto p            = grid.alloc_scalar();
  auto dp           = grid.alloc_scalar();
  auto div          = grid.alloc_scalar();

  auto T_old        = grid.alloc_scalar();
  auto T            = grid.alloc_scalar();
  auto Tsrc         = grid.alloc_scalar();
  auto FT           = grid.alloc_face_vector();

  Float dt          = 0.0;
  Float t           = 0.0;

  Index mg_cycles   = 0;
  Float mg_residual = 0.0;

  const BConds<Float> bconds{
      .left   = Periodic{},
      .right  = Periodic{},
      .bottom = Periodic{},
      .top    = Periodic{},
  };

  // = Linear solver ===============================================================================
  const std::array<int, 2> ns   = {grid.nx(), grid.ny()};
  const std::array<Float, 2> Ls = {grid.x_max() - grid.x_min(), grid.y_max() - grid.y_min()};
  // const std::array<int, 4> BCs  = {
  //     PoisFFT::NEUMANN_STAG, PoisFFT::NEUMANN_STAG, PoisFFT::NEUMANN_STAG,
  //     PoisFFT::NEUMANN_STAG};
  const std::array<int, 4> BCs = {
      PoisFFT::PERIODIC, PoisFFT::PERIODIC, PoisFFT::PERIODIC, PoisFFT::PERIODIC};
  PoisFFT::Solver<2, Float> fft_solver(
      ns.data(), Ls.data(), BCs.data(), PoisFFT::FINITE_DIFFERENCE_2);
  const std::array<int, 2> ngs = {grid.nghost(), grid.nghost()};

  // ~~~~~

  MultigridSolver mg_solver(grid, bconds);
  // = Linear solver ===============================================================================

  grid.foreach_face_i<Dimension::X>(
      FOREACH_FUNC { u.x(i, j) = u_analytical(grid.x(i), grid.ym(j), 0.0); });
  grid.foreach_face_i<Dimension::Y>(
      FOREACH_FUNC { u.y(i, j) = v_analytical(grid.xm(i), grid.y(j), 0.0); });
  apply_velocity_bconds(grid, bconds, bconds, u);
  interpolate(grid, u, ui);

  fill(T, 293.15);
  grid.foreach_i(FOREACH_FUNC {
    const auto r_sqr  = Igor::sqr(grid.xm(i) - 0.5) + Igor::sqr(grid.ym(j) - 0.5);
    T(i, j)          += std::exp(-20.0 * r_sqr);
  });

  VTKWriter writer(output_dir, grid);
  writer.add_field("u", ui);
  writer.add_field("p", p);
  writer.add_field("div", div);
  writer.add_field("T", T);
  writer.add_field("Tsrc", Tsrc);
  if (!writer.write(t)) { return 1; }

  Stats p_stats    = stats(grid, p);
  Stats u_stats    = stats(grid, u.x);
  Stats v_stats    = stats(grid, u.y);
  Stats T_stats    = stats(grid, T);
  Stats Tsrc_stats = stats(grid, Tsrc);
  Stats div_stats  = stats(grid, div);
  Float div_max    = std::max(std::abs(div_stats.min), std::abs(div_stats.max));

  Monitor<Float> monitor(output_dir + "/monitor.log");
  monitor.add_variable(&t, "t");
  monitor.add_variable(&dt, "dt");
  monitor.add_variable(&p_stats.max, "max(p)");
  monitor.add_variable(&u_stats.max, "max(u)");
  monitor.add_variable(&v_stats.max, "max(v)");
  monitor.add_variable(&T_stats.min, "min(T)");
  monitor.add_variable(&T_stats.max, "max(T)");
  monitor.add_variable(&Tsrc_stats.min, "min(Tsrc)");
  monitor.add_variable(&Tsrc_stats.max, "max(Tsrc)");
  monitor.add_variable(&div_max, "absmax(div)");
  if (multigrid) {
    monitor.add_variable(&mg_residual, "res(MG)");
    monitor.add_variable(&mg_cycles, "cycles(MG)");
  }
  monitor.write();

  IGOR_TIME_SCOPE("Solver")
  while (t < tend) {
    dt = std::min({
        adjust_dt(grid, u, rho, mu, CFL),
        advection_adjust_dt(grid, kappa / (rho * cV), CFL),
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
      apply_velocity_bconds(grid, bconds, bconds, u);

      // 2) Pressure correction
      calc_div(grid, u, div);
      grid.foreach_i(FOREACH_FUNC { div(i, j) *= rho / local_dt; });
      if (multigrid) {
        mg_solver.solve(dp, div, 1e-8);
        mg_cycles   = mg_solver.num_cycles();
        mg_residual = mg_solver.res();
      } else {
        fft_solver.execute(dp.data(), div.data(), ngs.data(), ngs.data());
      }
      apply_bconds(grid, bconds, dp, t);

      // 3) Project
      correct_velocity(grid, dp, rho, local_dt, u, p);

      // Update species
      calc_viscous_temperature_src(grid, u, mu, rho, cV, Tsrc);
      advection_calc_flux(grid, u, T, kappa / (rho * cV), FT);
      advection_update_s(grid, local_dt, FT, Tsrc, T_old, T);
      apply_bconds(grid, bconds, T, t);
    }

    interpolate(grid, u, ui);
    calc_div(grid, u, div);

    p_stats     = stats(grid, p);
    u_stats     = stats(grid, u.x);
    v_stats     = stats(grid, u.y);
    T_stats     = stats(grid, T);
    Tsrc_stats  = stats(grid, Tsrc);
    div_stats   = stats(grid, div);
    div_max     = std::max(std::abs(div_stats.min), std::abs(div_stats.max));

    t          += dt;
    if (should_save(t, dt, dt_write, tend)) {
      if (!writer.write(t)) { return 1; }
    }
    monitor.write();
  }

  Float L1_u = 0.0;
  grid.foreach_face_i<Dimension::X, Exec::SERIAL>([=, &L1_u](Index i, Index j) {
    const auto u_exp  = u_analytical(grid.x(i), grid.ym(j), t);
    L1_u             += std::abs(u_exp - u.x(i, j)) * grid.dv(i, j);
  });
  Float L1_v = 0.0;
  grid.foreach_face_i<Dimension::Y, Exec::SERIAL>([=, &L1_v](Index i, Index j) {
    const auto v_exp  = v_analytical(grid.xm(i), grid.y(j), t);
    L1_v             += std::abs(v_exp - u.y(i, j)) * grid.dv(i, j);
  });

  Igor::Info("L1(u) = {:.8e}", L1_u);
  Igor::Info("L1(v) = {:.8e}", L1_v);

  Igor::Info("Ok.");
}
