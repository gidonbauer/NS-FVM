#pragma once

#include <Igor/Math.hpp>

#include "Grid.hpp"

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
  Stats s = grid.transform_reduce_i(
      Stats<Float>{
          .min    = std::numeric_limits<Float>::max(),
          .max    = -std::numeric_limits<Float>::max(),
          .sum    = 0.0,
          .stddev = 0.0,
          .volume = 0.0,
      },
      FOREACH_FUNC->Stats<Float> {
        return Stats<Float>{
            .min    = f(i, j),
            .max    = f(i, j),
            .sum    = grid.dv(i, j) * f(i, j),
            .stddev = grid.dv(i, j) * Igor::sqr(f(i, j)),
            .volume = grid.dv(i, j),
        };
      },
      [](Stats<Float> lhs, const Stats<Float>& rhs) -> Stats<Float> {
        lhs.min     = std::min(lhs.min, rhs.min);
        lhs.max     = std::max(lhs.max, rhs.max);
        lhs.sum    += rhs.sum;
        lhs.stddev += rhs.stddev;
        lhs.volume += rhs.volume;
        return lhs;
      });

  if (s.volume > 0.0) { s.stddev -= s.sum * s.sum / s.volume; }
  if (s.stddev > 0.0) { s.stddev = std::sqrt(s.stddev / s.volume); }
  return s;
}

template <typename Float, Layout LAYOUT>
constexpr auto stats(const Grid<Float, LAYOUT>& grid,
                     const Scalar<Float, LAYOUT> f,
                     const Scalar<Float, LAYOUT> metric) -> Stats<Float> {
  Stats s = grid.transform_reduce_i(
      Stats<Float>{
          .min    = std::numeric_limits<Float>::max(),
          .max    = -std::numeric_limits<Float>::max(),
          .sum    = 0.0,
          .stddev = 0.0,
          .volume = 0.0,
      },
      FOREACH_FUNC->Stats<Float> {
        return Stats<Float>{
            .min    = f(i, j),
            .max    = f(i, j),
            .sum    = metric(i, j) * grid.dv(i, j) * f(i, j),
            .stddev = metric(i, j) * grid.dv(i, j) * Igor::sqr(f(i, j)),
            .volume = metric(i, j) * grid.dv(i, j),
        };
      },
      [](Stats<Float> lhs, const Stats<Float>& rhs) -> Stats<Float> {
        lhs.min     = std::min(lhs.min, rhs.min);
        lhs.max     = std::max(lhs.max, rhs.max);
        lhs.sum    += rhs.sum;
        lhs.stddev += rhs.stddev;
        lhs.volume += rhs.volume;
        return lhs;
      });

  if (s.volume > 0.0) { s.stddev -= s.sum * s.sum / s.volume; }
  if (s.stddev > 0.0) { s.stddev = std::sqrt(s.stddev / s.volume); }
  return s;
}
