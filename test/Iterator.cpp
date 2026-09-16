#include <functional>
#include <numbers>

#include <Igor/Math.hpp>
#include <Igor/Timer.hpp>

#include "BoundaryConditions.hpp"
#include "Common.hpp"
#include "Grid.hpp"

using Float        = double;
constexpr Float pi = std::numbers::pi_v<Float>;

constexpr auto approx_eq(Float lhs, Float rhs, Float abstol) -> bool {
  return std::abs(lhs - rhs) <= abstol;
}

// =================================================================================================
template <Exec EXEC>
auto test_minmax_reduce() -> bool {
  const Index N = 1 << 14;
  Grid<Float> grid(0.0, 1.0, N, 0.0, 1.0, N, 1);
  auto s = grid.alloc_scalar();

  grid.foreach_i<EXEC>(FOREACH_FUNC { s(i, j) = static_cast<Float>(i * N + j); });
  s(-1, -1) = -1.0;
  s(N, N)   = std::numeric_limits<Float>::max();

  struct MinMax {
    Float min, max;
  };
  bool any_failed = false;

  // Internal
  {
    const auto res = grid.transform_reduce_i<EXEC>(
        MinMax{
            .min = std::numeric_limits<Float>::max(),
            .max = -std::numeric_limits<Float>::max(),
        },
        FOREACH_FUNC { return MinMax{.min = s(i, j), .max = s(i, j)}; },
        [](MinMax lhs, const MinMax& rhs) {
          lhs.min = std::min(lhs.min, rhs.min);
          lhs.max = std::max(lhs.max, rhs.max);
          return lhs;
        });

    constexpr Float min_expected = 0.0;
    constexpr auto max_expected  = static_cast<Float>((N - 1) * N + (N - 1));
    if (res.min != min_expected) {
      Igor::Error("Incorrect min value in interior, expected {} but got {}", min_expected, res.min);
      any_failed = true;
    }
    if (res.max != max_expected) {
      Igor::Error("Incorrect max value in interior, expected {} but got {}", max_expected, res.max);
      any_failed = true;
    }
  }

  // All
  {
    const auto res = grid.transform_reduce_a<EXEC>(
        MinMax{
            .min = std::numeric_limits<Float>::max(),
            .max = -std::numeric_limits<Float>::max(),
        },
        FOREACH_FUNC { return MinMax{.min = s(i, j), .max = s(i, j)}; },
        [](MinMax lhs, const MinMax& rhs) {
          lhs.min = std::min(lhs.min, rhs.min);
          lhs.max = std::max(lhs.max, rhs.max);
          return lhs;
        });

    constexpr Float min_expected = -1.0;
    constexpr auto max_expected  = std::numeric_limits<Float>::max();
    if (res.min != min_expected) {
      Igor::Error(
          "Incorrect min value in entire field, expected {} but got {}", min_expected, res.min);
      any_failed = true;
    }
    if (res.max != max_expected) {
      Igor::Error(
          "Incorrect max value in entire field, expected {} but got {}", max_expected, res.max);
      any_failed = true;
    }
  }

  return !any_failed;
}

// =================================================================================================
auto test_residual() -> bool {
  const Index N = 1 << 14;
  Grid<Float> grid(0.0, 1.0, N, 0.0, 1.0, N, 1);
  auto sol            = grid.alloc_scalar();
  auto rhs            = grid.alloc_scalar();
  auto res            = grid.alloc_scalar();

  const Float inv_dx2 = 1.0 / Igor::sqr(grid.dx());
  const Float inv_dy2 = 1.0 / Igor::sqr(grid.dy());

  grid.foreach_i(FOREACH_FUNC {
    const auto x = grid.xm(i);
    const auto y = grid.ym(j);
    sol(i, j)    = std::cos(2.0 * pi * x) * std::cos(2.0 * pi * y);
  });
  apply_neumann_bconds(grid, sol);

  grid.foreach_i(FOREACH_FUNC {
    rhs(i, j) = (sol(i - 1, j) - 2.0 * sol(i, j) + sol(i + 1, j)) * inv_dx2 +
                (sol(i, j - 1) - 2.0 * sol(i, j) + sol(i, j + 1)) * inv_dy2;
  });
  rhs(0, 0)     += 1.0;

  Float max_res  = grid.transform_reduce_i(
      0.0,
      FOREACH_FUNC {
        const Float L = (sol(i - 1, j) - 2.0 * sol(i, j) + sol(i + 1, j)) * inv_dx2 +
                        (sol(i, j - 1) - 2.0 * sol(i, j) + sol(i, j + 1)) * inv_dy2;
        res(i, j)     = rhs(i, j) - L;
        return std::abs(res(i, j));
      },
      [](Float lhs, Float rhs) { return std::max(lhs, rhs); });

  if (!approx_eq(max_res, 1.0, 1e-12)) {
    Igor::Error("Expected max_res to be 1.0 but is {:.12e}.", max_res);
    return false;
  }

  return true;
}

// =================================================================================================
template <Exec EXEC>
auto test_sum() -> bool {
  const Index N = 1 << 14;
  Grid<Float> grid(0.0, 1.0, N, 0.0, 1.0, N, 1);
  auto s = grid.alloc_scalar();

  grid.foreach_i<EXEC>(FOREACH_FUNC { s(i, j) = 1.0; });
  Float sum_plus_one =
      grid.transform_reduce_i<EXEC>(1.0, FOREACH_FUNC { return s(i, j); }, std::plus<>{});

  constexpr Float sum_plus_one_expected = N * N + 1.0;
  if (sum_plus_one != sum_plus_one_expected) {
    Igor::Error("Expected sum_plus_one to be {} but is {}.", sum_plus_one_expected, sum_plus_one);
    return false;
  }

  return true;
}

// =================================================================================================
auto test_stats() -> bool {
  constexpr Index N     = 1 << 14;
  constexpr Float x_min = -1.0;
  constexpr Float x_max = 1.0;
  constexpr Float y_min = -1.0;
  constexpr Float y_max = 1.0;
  Grid<Float> grid(x_min, x_max, N, y_min, y_max, N, 1);
  auto s = grid.alloc_scalar();

  fill(s, 1.0);
  s(0, 0)                   = 0.0;
  s(s.nx() - 1, s.ny() - 1) = 2.0;

  const auto s_stats        = stats(grid, s);

  bool any_failed           = false;
  if (const auto exp_vol = (x_max - x_min) * (y_max - y_min);
      !approx_eq(s_stats.volume, exp_vol, 1e-12)) {
    Igor::Error("Calculated volume {:.12e} does not match expected volume {:.12e}",
                s_stats.volume,
                exp_vol);
    any_failed = true;
  }
  if (const auto exp_min = 0.0; !approx_eq(s_stats.min, exp_min, 1e-12)) {
    Igor::Error(
        "Calculated minimum {:.12e} does not match expected minimum {:.12e}", s_stats.min, exp_min);
    any_failed = true;
  }
  if (const auto exp_max = 2.0; !approx_eq(s_stats.max, exp_max, 1e-12)) {
    Igor::Error(
        "Calculated maximum {:.12e} does not match expected maximum {:.12e}", s_stats.max, exp_max);
    any_failed = true;
  }
  if (const auto exp_sum = (x_max - x_min) * (y_max - y_min);
      !approx_eq(s_stats.sum, exp_sum, 1e-12)) {
    Igor::Error("Calculated sum {:.12e} does not match expected sum {:.12e}", s_stats.sum, exp_sum);
    any_failed = true;
  }
  if (const auto exp_stddev = std::sqrt(2.0 / (grid.nx() * grid.ny()));
      !approx_eq(s_stats.stddev, exp_stddev, 1e-12)) {
    Igor::Error(
        "Calculated standard deviation {:.12e} does not match expected standard deviation {:.12e}",
        s_stats.stddev,
        exp_stddev);
    any_failed = true;
  }

  return !any_failed;
}

// =================================================================================================
auto main() -> int {
  bool any_failed = false;

  IGOR_TIME_SCOPE("MinMax serial")
  if (!test_minmax_reduce<Exec::SERIAL>()) {
    Igor::Error("test_minmax_reduce(serial) failed.");
    any_failed = true;
  }

  IGOR_TIME_SCOPE("MinMax parallel")
  if (!test_minmax_reduce<Exec::PARALLEL>()) {
    Igor::Error("test_minmax_reduce(parallel) failed.");
    any_failed = true;
  }

  IGOR_TIME_SCOPE("Residual")
  if (!test_residual()) {
    Igor::Error("test_residual failed.");
    any_failed = true;
  }

  IGOR_TIME_SCOPE("Sum serial")
  if (!test_sum<Exec::SERIAL>()) {
    Igor::Error("test_sum(serial) failed.");
    any_failed = true;
  }

  IGOR_TIME_SCOPE("Sum parallel")
  if (!test_sum<Exec::PARALLEL>()) {
    Igor::Error("test_sum(parallel) failed.");
    any_failed = true;
  }

  IGOR_TIME_SCOPE("Stats")
  if (!test_stats()) {
    Igor::Error("test_stats failed.");
    any_failed = true;
  }

  return any_failed ? 1 : 0;
}
