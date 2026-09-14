#include <algorithm>
#include <cerrno>
#include <charconv>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <numbers>
#include <numeric>
#include <string>
#include <string_view>
#include <vector>

#include <Igor/Defer.hpp>
#include <Igor/Logging.hpp>

#include "Mac.hpp"
#include "MultigridPoisson.hpp"

using Float             = double;
constexpr Layout LAYOUT = Layout::C;

using GridType          = Grid<Float, LAYOUT>;
using ScalarType        = decltype(std::declval<GridType>().alloc_scalar());
using FaceVectorType    = decltype(std::declval<GridType>().alloc_face_vector());
using SolverType        = MultigridSolver<Float, LAYOUT>;

struct Problem {
  std::string name;
  ScalarType rhs;
  ScalarType sol0;
  Float tol;
  Float dt;
  Index num_iter_post;
};

struct Timing {
  Float min;
  Float median;
  Float mean;
  Float max;
};

struct BenchResult {
  Timing time;
  Index cycles;
  Index num_iter_pre;
  Index num_iter_post;
  Float res;
  bool converged;
};

[[nodiscard]] auto summarize(std::vector<Float> samples) noexcept -> Timing {
  IGOR_ASSERT(!samples.empty(), "Cannot summarize an empty set of samples.");
  std::sort(samples.begin(), samples.end());
  const auto n   = static_cast<Float>(samples.size());
  const auto sum = std::accumulate(samples.begin(), samples.end(), Float{0.0});
  return {
      .min    = samples.front(),
      .median = samples[samples.size() / 2],
      .mean   = sum / n,
      .max    = samples.back(),
  };
}

[[nodiscard]] auto
bench_solve(SolverType& solver, const Problem& problem, ScalarType sol, Index reps) -> BenchResult {
  copy(problem.sol0, sol);
  solver.num_iter_post() = problem.num_iter_post;
  static_cast<void>(solver.solve(sol, problem.rhs, problem.tol));

  std::vector<Float> samples;
  samples.reserve(static_cast<size_t>(reps));

  BenchResult result{};
  for (Index rep = 0; rep < reps; ++rep) {
    copy(problem.sol0, sol);
    solver.num_iter_post() = problem.num_iter_post;

    const auto t_begin     = std::chrono::steady_clock::now();
    const bool ok          = solver.solve(sol, problem.rhs, problem.tol);
    const auto t_end       = std::chrono::steady_clock::now();

    samples.push_back(std::chrono::duration<Float>(t_end - t_begin).count());
    result.cycles        = solver.num_cycles();
    result.num_iter_pre  = solver.num_iter_pre();
    result.num_iter_post = solver.num_iter_post();
    result.res           = solver.res();
    result.converged     = ok;
  }

  result.time = summarize(std::move(samples));
  return result;
}

void report(const GridType& grid,
            const Problem& problem,
            const BenchResult& result,
            std::FILE* file) {
  constexpr Float TO_MS = 1e3;

  // NOLINTBEGIN
  std::fprintf(file, "%s (dt = %.4e, tol = %.4e)\n", problem.name.c_str(), problem.dt, problem.tol);
  std::fprintf(file, "  cycles        = %d\n", result.cycles);
  std::fprintf(file, "  num_iter_pre  = %d\n", result.num_iter_pre);
  std::fprintf(file, "  num_iter_post = %d\n", result.num_iter_post);
  std::fprintf(file, "  res           = %.4e\n", result.res);
  std::fprintf(file, "  converged     = %s\n", result.converged ? "yes" : "no");
  std::fprintf(file, "  wall-time:\n");
  std::fprintf(file, "    min    = %.4f ms\n", result.time.min * TO_MS);
  std::fprintf(file, "    median = %.4f ms\n", result.time.median * TO_MS);
  std::fprintf(file, "    mean   = %.4f ms\n", result.time.mean * TO_MS);
  std::fprintf(file, "    max    = %.4f ms\n", result.time.max * TO_MS);

  if (result.cycles == 0) {
    std::fprintf(file,
                 "    [WARN] zero cycles: tolerance met by the initial guess, timing is a residual "
                 "evaluation only. Lower tol_scale to load the solver.\n");
    return;
  }

  const auto cells        = static_cast<Float>(grid.nx()) * static_cast<Float>(grid.ny());
  const auto cycles       = static_cast<Float>(result.cycles);
  const Float per_cycle   = result.time.median / cycles;
  const Float cell_cycles = cells * cycles / result.time.median * 1e-6;
  std::fprintf(file, "  wall-time per cycle = %.4f ms\n", per_cycle * TO_MS);
  std::fprintf(file, "  throughput          = %.4f Mcell-cycles/s\n", cell_cycles);
  std::fputc('\n', file);
  // NOLINTEND
}

namespace CartesianCase {

constexpr Float x_min = 0.0;
constexpr Float x_max = 1.0;
constexpr Float y_min = 0.0;
constexpr Float y_max = 1.0;

constexpr Float rho   = 1.0;   // [kg/m^3]
constexpr Float mu    = 1e-5;  // [Pa*s]

constexpr Float L     = x_max - x_min;
constexpr Float Re    = 600.0;
constexpr Float Uwall = Re * mu / (rho * L);

constexpr Float CFL   = 0.5;

void run(Index N, Float abstol, Index reps, Index spinup, std::FILE* file) {
  GridType grid(x_min, x_max, N, y_min, y_max, N, 1);

  auto u_old = grid.alloc_face_vector();
  auto u     = grid.alloc_face_vector();

  auto FUX   = grid.alloc_scalar();
  auto FUY   = grid.alloc_vertex_scalar();
  auto FVX   = grid.alloc_vertex_scalar();
  auto FVY   = grid.alloc_scalar();

  auto p     = grid.alloc_scalar();
  auto dp    = grid.alloc_scalar();
  auto div   = grid.alloc_scalar();

  Float dt   = 0.0;
  Float t    = 0.0;

  const BConds<Float> u_bconds{
      .left   = Dirichlet<Float>{.val = 0.0},
      .right  = Dirichlet<Float>{.val = 0.0},
      .bottom = Dirichlet<Float>{.val = 0.0},
      .top    = Dirichlet<Float>{.val = Uwall},
  };
  const BConds<Float> v_bconds{
      .left   = Dirichlet<Float>{.val = 0.0},
      .right  = Dirichlet<Float>{.val = 0.0},
      .bottom = Dirichlet<Float>{.val = 0.0},
      .top    = Dirichlet<Float>{.val = 0.0},
  };
  const BConds<Float> dp_bconds{
      .left   = Neumann{},
      .right  = Neumann{},
      .bottom = Neumann{},
      .top    = Neumann{},
  };
  SolverType solver(grid, dp_bconds);

  Problem first{
      .name          = "Cartesian, initial step, cold start",
      .rhs           = grid.alloc_scalar(),
      .sol0          = grid.alloc_scalar(),
      .tol           = 0.0,
      .dt            = 0.0,
      .num_iter_post = solver.num_iter_post(),
  };
  Problem developed{
      .name          = "Cartesian, developed profile, warm start",
      .rhs           = grid.alloc_scalar(),
      .sol0          = grid.alloc_scalar(),
      .tol           = 0.0,
      .dt            = 0.0,
      .num_iter_post = solver.num_iter_post(),
  };

  for (Index iter = 0; iter < spinup; ++iter) {
    dt = adjust_dt(grid, u, rho, mu, CFL);

    copy(u, u_old);

    for (Index sub_iter = 0; sub_iter < 2; ++sub_iter) {
      const auto local_dt = sub_iter == 0 ? dt / 2.0 : dt;

      // 1) Predictor
      calc_flux(grid, u, p, rho, mu, FUX, FUY, FVX, FVY);
      update_u(grid, local_dt, FUX, FUY, FVX, FVY, u_old, u);
      apply_velocity_bconds(grid, u_bconds, v_bconds, u, t);

      // 2) Pressure correction
      calc_div(grid, u, div);
      grid.foreach_i(FOREACH_FUNC { div(i, j) *= rho / local_dt; });
      const auto tol = abstol / Igor::sqr(local_dt);
      if (iter == 0 && sub_iter == 1) {
        copy(div, first.rhs);
        fill(first.sol0, 0.0);
        first.tol = tol;
        first.dt  = local_dt;
      }
      if (iter + 1 == spinup && sub_iter == 1) {
        copy(div, developed.rhs);
        fill(developed.sol0, 0.0);
        developed.tol           = tol;
        developed.dt            = local_dt;
        developed.num_iter_post = solver.num_iter_post();
      }
      solver.solve(dp, div, tol);
      apply_bconds(grid, dp_bconds, dp, t);

      // 3) Project
      correct_velocity(grid, dp, rho, local_dt, u, p);
    }

    t += dt;
  }

  auto sol = grid.alloc_scalar();
  for (const Problem& problem : {first, developed}) {
    const auto result = bench_solve(solver, problem, sol, reps);
    report(grid, problem, result, file);
  }
}

}  // namespace CartesianCase

namespace PolarCase {

constexpr Float theta_min = 0.0;
constexpr Float theta_max = 2.0 * std::numbers::pi_v<Float>;
constexpr Float r_min     = 1.0;
constexpr Float r_max     = 10.0;

constexpr Float Uinf      = 1.0;
constexpr Float rho       = 1.0;
constexpr Float mu        = 1e-3;

constexpr Float CFL       = 0.7;

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

void run(Index N, Float abstol, Index reps, Index spinup, std::FILE* file) {
  GridType grid(theta_min, theta_max, N, r_min, r_max, N, 1, Coordinates::POLAR);

  auto u_old = grid.alloc_face_vector();
  auto u     = grid.alloc_face_vector();

  auto FUX   = grid.alloc_scalar();
  auto FUY   = grid.alloc_vertex_scalar();
  auto FVX   = grid.alloc_vertex_scalar();
  auto FVY   = grid.alloc_scalar();

  auto p     = grid.alloc_scalar();
  auto dp    = grid.alloc_scalar();
  auto div   = grid.alloc_scalar();

  Float dt   = 0.0;
  Float t    = 0.0;

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
  SolverType solver(grid, dp_bconds);

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

  Problem first{
      .name          = "Polar, initial step, cold start",
      .rhs           = grid.alloc_scalar(),
      .sol0          = grid.alloc_scalar(),
      .tol           = 0.0,
      .dt            = 0.0,
      .num_iter_post = solver.num_iter_post(),
  };
  Problem developed{
      .name          = "Polar, developed profile, warm start",
      .rhs           = grid.alloc_scalar(),
      .sol0          = grid.alloc_scalar(),
      .tol           = 0.0,
      .dt            = 0.0,
      .num_iter_post = solver.num_iter_post(),
  };

  for (Index iter = 0; iter < spinup; ++iter) {
    dt = adjust_dt(grid, u, rho, mu, CFL);

    copy(u, u_old);

    for (Index sub_iter = 0; sub_iter < 2; ++sub_iter) {
      const auto local_dt = sub_iter == 0 ? dt / 2.0 : dt;

      // 1) Predictor
      calc_flux(grid, u, p, rho, mu, FUX, FUY, FVX, FVY);
      update_u(grid, local_dt, FUX, FUY, FVX, FVY, u_old, u);
      apply_velocity_bconds(grid, u_bconds, v_bconds, u, t);
      custom_velocity_top_boundary(grid, u);

      // 2) Pressure correction
      calc_div(grid, u, div);
      grid.foreach_i(FOREACH_FUNC { div(i, j) *= rho / local_dt; });
      const auto tol = abstol / Igor::sqr(local_dt);
      if (iter == 0 && sub_iter == 1) {
        copy(div, first.rhs);
        fill(first.sol0, 0.0);
        first.tol = tol;
        first.dt  = local_dt;
      }
      if (iter + 1 == spinup && sub_iter == 1) {
        copy(div, developed.rhs);
        fill(developed.sol0, 0.0);
        developed.tol           = tol;
        developed.dt            = local_dt;
        developed.num_iter_post = solver.num_iter_post();
      }
      solver.solve(dp, div, tol);
      apply_bconds(grid, dp_bconds, dp, t);

      // 3) Project
      correct_velocity(grid, dp, rho, local_dt, u, p);
    }

    t += dt;
  }

  auto sol = grid.alloc_scalar();
  for (const Problem& problem : {first, developed}) {
    const auto result = bench_solve(solver, problem, sol, reps);
    report(grid, problem, result, file);
  }
}

}  // namespace PolarCase

// = Command line parsing ==========================================================================
[[nodiscard]] auto parse_index(std::string_view str, Index& out) noexcept -> bool {
  const auto* end = str.data() + str.size();
  const auto res  = std::from_chars(str.data(), end, out);
  return res.ec == std::errc{} && res.ptr == end && out > 0;
}

[[nodiscard]] auto parse_float(std::string_view str, Float& out) -> bool {
  const std::string buf(str);
  if (buf.empty()) { return false; }

  char* end    = nullptr;
  const auto v = std::strtod(buf.c_str(), &end);
  if (end != buf.c_str() + buf.size() || !(v > 0.0)) { return false; }

  out = v;
  return true;
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
  Index N      = 512;
  Float abstol = 1e-4;
  Index reps   = 20;
  Index spinup = 10;
  std::string filename;

  const auto* prog     = pop_arg(argc, argv);
  const auto usage_str = Igor::detail::format(
      "Usage: {} [-N=<N>] [--abstol=<abstol>] [--reps=<reps>] [--spinup=<spinup>] [-o=<output>]\n"
      "  -N        Grid size per direction                         (default: {})\n"
      "  --abstol  Absolute tolerance; the solver is asked for     (default: {:g})\n"
      "            `abstol / dt^2`, as the drivers do\n"
      "  --reps    Timed repetitions per problem                   (default: {})\n"
      "  --spinup  Timesteps to run before capturing the problems  (default: {})\n"
      "  -o        Write the report to this file                   (default: stdout)\n"
      "  -h        Show this message\n"
      "\n"
      "Every flag takes its value as `-flag=value` or as `-flag value`, and accepts either one or\n"
      "two leading dashes.",
      prog,
      N,
      abstol,
      reps,
      spinup);

  while (argc > 0) {
    const std::string_view arg  = pop_arg(argc, argv);
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
    if (name == "N") {
      ok = parse_index(value, N);
    } else if (name == "abstol") {
      ok = parse_float(value, abstol);
    } else if (name == "reps") {
      ok = parse_index(value, reps);
    } else if (name == "spinup") {
      ok = parse_index(value, spinup);
    } else if (name == "o" || name == "output") {
      ok = !value.empty();
      if (ok) { filename = value; }
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

  std::FILE* out = stdout;
  if (!filename.empty()) {
    out = std::fopen(filename.c_str(), "w");  // NOLINT
    if (out == nullptr) {
      Igor::Error("Could not open file `{}`: {}", filename, std::strerror(errno));
      return 1;
    }
  }
  IGOR_DEFER(if (out != stdout) { std::fclose(out); });  // NOLINT

  // NOLINTBEGIN
  std::fprintf(out, "N        = %d\n", N);
  std::fprintf(out, "abstol   = %g\n", abstol);
  std::fprintf(out, "reps     = %d\n", reps);
  std::fprintf(out, "spinup   = %d\n", spinup);
  std::fprintf(out,
               "parallel = %s\n",
#ifdef NS_FVM_PARALLEL
               "yes"
#else
               "no"
#endif
  );
  std::fputc('\n', out);
  // NOLINTEND

  PolarCase::run(N, abstol, reps, spinup, out);
  CartesianCase::run(N, abstol, reps, spinup, out);

  if (out != stdout) { Igor::Info("Wrote benchmark to `{}`", filename); }
}
