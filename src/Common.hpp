#pragma once

#include <Igor/Math.hpp>

#include "Grid.hpp"

template <typename Float>
struct Stats {
  Float min;
  Float max;
  Float sum;
  Float stddev;
  Float volume;
};

template <typename Float, Layout LAYOUT>
constexpr auto stats(const Grid<Float, LAYOUT>& grid, const Scalar<Float, LAYOUT> f)
    -> Stats<Float> {
  Float min    = std::numeric_limits<Float>::max();
  Float max    = -std::numeric_limits<Float>::max();
  Float sum    = 0.0;
  Float sum2   = 0.0;
  Float volume = 0.0;
  grid.foreach_i([=, &min, &max, &sum, &sum2, &volume](Index i, Index j) {
    volume += grid.dv(i, j);
    sum    += grid.dv(i, j) * f(i, j);
    sum2   += grid.dv(i, j) * Igor::sqr(f(i, j));
    max     = std::max(max, f(i, j));
    min     = std::min(min, f(i, j));
  });

  Stats s{
      .min    = min,
      .max    = max,
      .sum    = sum,
      .stddev = 0.0,
      .volume = volume,
  };
  if (volume > 0.0) { sum2 -= sum * sum / volume; }
  if (sum2 > 0.0) { s.stddev = std::sqrt(sum2 / volume); }
  return s;
}
