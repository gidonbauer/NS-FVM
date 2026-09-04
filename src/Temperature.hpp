#pragma once

#include <Igor/Math.hpp>

#include "Grid.hpp"

template <typename Float, Layout LAYOUT>
constexpr void calc_viscous_temperature_src(const Grid<Float, LAYOUT>& grid,
                                            const FaceVector<Float, LAYOUT> u,
                                            Float mu,
                                            Float rho,
                                            Float cV,
                                            Scalar<Float, LAYOUT> src) {
  // mu/(rho cV) * [nabla u : nabla u + (nabla u)^T : nabla u]
  grid.foreach_i(FOREACH_FUNC {
    const auto dudx = (u.right(i, j) - u.left(i, j)) / grid.dx();
    const auto dvdy = (u.top(i, j) - u.bottom(i, j)) / grid.dy();
    const auto dudy = ((u.right(i, j + 1) + u.left(i, j + 1)) / 2.0 -
                       (u.right(i, j - 1) + u.left(i, j - 1)) / 2.0) /
                      (2.0 * grid.dy());
    const auto dvdx = ((u.top(i + 1, j) + u.bottom(i + 1, j)) / 2.0 -
                       (u.top(i - 1, j) + u.bottom(i - 1, j)) / 2.0) /
                      (2.0 * grid.dx());
    src(i, j) = mu / (rho * cV) *
                (2.0 * Igor::sqr(dudx) + Igor::sqr(dudy) + Igor::sqr(dvdx) + 2.0 * Igor::sqr(dvdy) +
                 2.0 * dvdx * dudy);
  });
}
