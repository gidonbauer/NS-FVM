#include <charconv>
#include <numbers>

#include <Igor/Timer.hpp>

#include "BoundaryConditions.hpp"
#include "Common.hpp"
#include "Grid.hpp"
#include "HDFWriter.hpp"
#include "IO.hpp"
#include "Mac.hpp"
#include "Monitor.hpp"
#include "MultigridPoisson.hpp"

// = Setup =========================================================================================
using Float               = double;

constexpr Float theta_min = 0.0;
constexpr Float theta_max = 2.0 * std::numbers::pi_v<Float>;
constexpr Float r_min     = 1.0;
constexpr Float r_max     = 10.0;

constexpr Float Uinf      = 1.0;
constexpr Float rho       = 1.0;
constexpr Float mu        = 1e-3;

constexpr Float Re        = Uinf * rho * r_min / mu;

constexpr Float CFL       = 0.7;
constexpr Float tend      = 25.0;
constexpr Float dt_write  = tend / 10.0;
// = Setup =========================================================================================

// =================================================================================================
template <typename Float, Layout LAYOUT>
constexpr void custom_velocity_top_boundary(const Grid<Float, LAYOUT>& grid,
                                            FaceVector<Float, LAYOUT> u) {

  // - U -----------
  grid.foreach_range(
      -u.x.nghost(), u.x.nx() + u.x.nghost(), 0, 1, FOREACH_FUNC {
        if (i >= u.x.nx() * 3 / 8 && i <= u.x.nx() * 5 / 8) {
          // Linear extrapolation
          const auto sN    = u.x(i, u.x.ny() - 1);
          const auto theta = grid.x(i);
          const auto v     = Uinf * -std::sin(theta);
          for (j = u.x.ny(); j < u.x.ny() + u.x.nghost(); ++j) {
            u.x(i, j) = sN + 2.0 * (v - sN) * (j - u.x.ny() + 1);
          }
        } else {
          for (j = u.x.ny(); j < u.x.ny() + u.x.nghost(); ++j) {
            u.x(i, j) = u.x(i, 2 * u.x.ny() - j - 1);
          }
        }
      });

  // - V -----------
  grid.foreach_range(
      -u.y.nghost(), u.y.nx() + u.y.nghost(), 0, 1, FOREACH_FUNC {
        if (i >= u.y.nx() * 3 / 8 && i <= u.y.nx() * 5 / 8) {
          // Linear extrapolation
          const auto sN    = u.y(i, u.y.ny() - 2);
          const auto theta = grid.xm(i);
          const auto v     = Uinf * std::cos(theta);
          for (j = u.y.ny() - 1; j < u.y.ny() + u.y.nghost(); ++j) {
            u.y(i, j) = sN + (v - sN) * (j - u.y.ny() + 2);
          }
        } else {
          for (j = u.y.ny(); j < u.y.ny() + u.y.nghost(); ++j) {
            u.y(i, j) = u.y(i, 2 * u.y.ny() - j - 1);
          }
        }
      });
}

// =================================================================================================
[[nodiscard]] auto parse_index(std::string_view str, Index& out) noexcept -> bool {
  const auto* end = str.data() + str.size();
  const auto res  = std::from_chars(str.data(), end, out);
  return res.ec == std::errc{} && res.ptr == end;
}

[[nodiscard]] auto pop_arg(int& argc, char**& argv) -> char* {
  IGOR_ASSERT(argc > 0, "No arguments to pop.");
  argc -= 1;
  argv += 1;
  return argv[-1];
}

[[nodiscard]] constexpr auto strip_dashes(std::string_view arg) noexcept -> std::string_view {
  if (arg.starts_with("--")) { return arg.substr(2); }
  if (arg.starts_with("-")) { return arg.substr(1); }
  return {};
}

// =================================================================================================
auto main(int argc, char** argv) -> int {
  Index N              = -1;
  Index min_size       = 2;
  Index num_pre        = 0;
  Index num_post       = 4;
  const auto* prog     = pop_arg(argc, argv);
  const auto usage_str = Igor::detail::format(
      "Usage: {} [--min=<min>] [--pre=<pre>] [--post=<post>] <grid size>", prog);

  while (argc > 0) {
    const std::string_view arg = pop_arg(argc, argv);

    // Grid size
    if (!arg.starts_with('-')) {
      if (!parse_index(arg, N) || N <= 0) {
        Igor::Error("{}", usage_str);
        Igor::Error("  Invalid grid size `{}`", argv[1]);
        return 1;
      }
      continue;
    }

    const std::string_view flag = strip_dashes(arg);
    if (flag.empty()) {
      Igor::Error("{}", usage_str);
      Igor::Error("  Expected a flag but got `{}`", arg);
      return 1;
    }

    const auto eq               = flag.find('=');
    const std::string_view name = flag.substr(0, eq);

    if (name == "h" || name == "help") {
      Igor::Info("{}", usage_str);
      return 0;
    }

    std::string_view value;
    if (eq != std::string_view::npos) {
      value = flag.substr(eq + 1);
    } else if (argc > 0) {
      value = pop_arg(argc, argv);
    } else {
      Igor::Error("{}", usage_str);
      Igor::Error("  Flag `{}` expects a value", arg);
      return 1;
    }

    bool ok = true;
    if (name == "pre") {
      ok = parse_index(value, num_pre);
    } else if (name == "post") {
      ok = parse_index(value, num_post);
    } else if (name == "min") {
      ok = parse_index(value, min_size);
    } else {
      Igor::Error("{}", usage_str);
      Igor::Error("  Unknown flag `{}`", arg);
      return 1;
    }

    if (!ok) {
      Igor::Error("{}", usage_str);
      Igor::Error("  Invalid value `{}` for flag `{}`", value, name);
      return 1;
    }
  }

  if (N < 0) {
    Igor::Error("{}", usage_str);
    Igor::Error("  Did not provide grid size.");
    return 1;
  }

  Igor::Info("Re   = {}", Re);
  Igor::Info("Uinf = {}", Uinf);

  const auto output_dir = get_output_directory("bench/output");
  if (!init_output_directory(output_dir)) { return 1; }

  Grid<Float> grid(theta_min, theta_max, N, r_min, r_max, N, 1, Coordinates::POLAR);

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

  Float dt        = 0.0;
  Float t         = 0.0;

  Float iter_time = 0.0;

  const BConds<Float> u_bconds{
      .left   = Periodic{},
      .right  = Periodic{},
      .bottom = Dirichlet<Float>{.val = 0.0},
      .top =
          Dirichlet<Float>{.val = [](Float theta, Float /*t*/) { return Uinf * -std::sin(theta); }},
  };
  const BConds<Float> v_bconds{
      .left   = Periodic{},
      .right  = Periodic{},
      .bottom = Dirichlet<Float>{.val = 0.0},
      .top =
          Dirichlet<Float>{.val = [](Float theta, Float /*t*/) { return Uinf * std::cos(theta); }},
  };
  const BConds<Float> dp_bconds{
      .left   = Periodic{},
      .right  = Periodic{},
      .bottom = Neumann{},
      .top    = Neumann{},
  };

  MultigridSolver solver(grid, dp_bconds, min_size, num_pre, num_post);
  Index mg_cycles        = 0;
  Index mg_num_iter_pre  = solver.num_iter_pre();
  Index mg_num_iter_post = solver.num_iter_post();
  Float mg_res           = 0.0;

  grid.foreach_face_i<Dimension::X>(FOREACH_FUNC {
    const auto theta = grid.x(i);
    u.x(i, j)        = Uinf * -std::sin(theta);
  });
  grid.foreach_face_i<Dimension::Y>(FOREACH_FUNC {
    const auto theta = grid.xm(i);
    u.y(i, j)        = Uinf * std::cos(theta);
  });
  apply_velocity_bconds(grid, u_bconds, v_bconds, u);
  custom_velocity_top_boundary(grid, u);
  interpolate(grid, u, ui);

  HDFWriter writer(output_dir, grid);
  writer.add_field("u", ui);
  writer.add_field("p", p);
  writer.add_field("div", div);
  if (!writer.write(t)) { return 1; }

  Stats p_stats   = stats(grid, p);
  Stats u_stats   = stats(grid, u.x);
  Stats v_stats   = stats(grid, u.y);
  Stats div_stats = stats(grid, div);
  Float div_max   = std::max(std::abs(div_stats.min), std::abs(div_stats.max));

  Monitor<Float> monitor(output_dir + "/monitor.log");
  monitor.add_variable(&t, "t");
  monitor.add_variable(&dt, "dt");
  monitor.add_variable(&p_stats.max, "max(p)");
  monitor.add_variable(&u_stats.max, "max(u_theta)");
  monitor.add_variable(&v_stats.max, "max(u_r)");
  monitor.add_variable(&div_max, "absmax(div)");
  monitor.add_variable(&mg_res, "res(MG)");
  monitor.add_variable(&mg_cycles, "cycles(MG)");
  monitor.add_variable(&mg_num_iter_pre, "iter_pre(MG)");
  monitor.add_variable(&mg_num_iter_post, "iter_post(MG)");
  monitor.add_variable(&iter_time, "time(iter) [s]");
  monitor.write();

  IGOR_TIME_SCOPE("Solver")
  while (t < tend) {
    const auto t_begin = std::chrono::high_resolution_clock::now();

    dt                 = std::min({
        adjust_dt(grid, u, rho, mu, CFL),
        tend - t,
    });

    copy(u, u_old);

    for (Index sub_iter = 0; sub_iter < 2; ++sub_iter) {
      const auto local_dt = sub_iter == 0 ? dt / 2.0 : dt;

      // 1) Predictor
      calc_flux(grid, u, p, rho, mu, FUX, FUY, FVX, FVY);
      update_u(grid, local_dt, FUX, FUY, FVX, FVY, u_old, u);
      apply_velocity_bconds(grid, u_bconds, v_bconds, u);
      custom_velocity_top_boundary(grid, u);

      // 2) Pressure correction
      calc_div(grid, u, div);
      grid.foreach_i(FOREACH_FUNC { div(i, j) *= rho / local_dt; });
      if (!solver.solve(dp, div, 1e-3 / Igor::sqr(dt))) {
        Igor::Warn("t={:.8f}: Multigrid solver did not converge after {} cycles: res = {:.8e}",
                   t,
                   solver.num_cycles(),
                   solver.res());
      }
      mg_cycles        = solver.num_cycles();
      mg_num_iter_pre  = solver.num_iter_pre();
      mg_num_iter_post = solver.num_iter_post();
      mg_res           = solver.res();
      apply_bconds(grid, dp_bconds, dp, t);

      // 3) Project
      correct_velocity(grid, dp, rho, local_dt, u, p);
    }

    interpolate(grid, u, ui);
    calc_div(grid, u, div);

    p_stats    = stats(grid, p);
    u_stats    = stats(grid, u.x);
    v_stats    = stats(grid, u.y);
    div_stats  = stats(grid, div);
    div_max    = std::max(std::abs(div_stats.min), std::abs(div_stats.max));

    t         += dt;
    if (should_save(t, dt, dt_write, tend)) {
      if (!writer.write(t)) { return 1; }
    }

    iter_time =
        std::chrono::duration<Float>(std::chrono::high_resolution_clock::now() - t_begin).count();
    monitor.write();
  }

  Igor::Info("Ok.");
}
