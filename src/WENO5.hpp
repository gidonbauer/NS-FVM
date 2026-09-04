#pragma once

#include <Igor/Math.hpp>

#include "Grid.hpp"

template <typename Float>
constexpr auto weno_left_biased_reconstruction(Float um2, Float um1, Float u0, Float up1, Float up2)
    -> Float {
  constexpr Float eps  = 1e-6;
  constexpr Float dL_0 = 1.0 / 10.0;
  constexpr Float dL_1 = 6.0 / 10.0;
  constexpr Float dL_2 = 3.0 / 10.0;

  const auto q0        = (2.0 * um2 - 7.0 * um1 + 11.0 * u0) / 6.0;
  const auto q1        = (-1.0 * um1 + 5.0 * u0 + 2.0 * up1) / 6.0;
  const auto q2        = (2.0 * u0 + 5.0 * up1 - 1.0 * up2) / 6.0;

  const auto b0        = (13.0 / 12.0) * Igor::sqr(um2 - 2.0 * um1 + u0) +
                         (1.0 / 4.0) * Igor::sqr(um2 - 4.0 * um1 + 3.0 * u0);
  const auto b1        = (13.0 / 12.0) * Igor::sqr(um1 - 2.0 * u0 + up1) +  //
                         (1.0 / 4.0) * Igor::sqr(um1 - up1);
  const auto b2        = (13.0 / 12.0) * Igor::sqr(u0 - 2.0 * up1 + up2) +
                         (1.0 / 4.0) * Igor::sqr(3.0 * u0 - 4.0 * up1 + up2);

  const auto a0        = dL_0 / Igor::sqr(b0 + eps);
  const auto a1        = dL_1 / Igor::sqr(b1 + eps);
  const auto a2        = dL_2 / Igor::sqr(b2 + eps);
  const auto asum      = a0 + a1 + a2;

  return (a0 * q0 + a1 * q1 + a2 * q2) / asum;
}

template <typename Float>
constexpr auto
weno_right_biased_reconstruction(Float um1, Float u0, Float up1, Float up2, Float up3) -> Float {
  constexpr Float eps  = 1e-6;
  constexpr Float dR_0 = 3.0 / 10.0;
  constexpr Float dR_1 = 6.0 / 10.0;
  constexpr Float dR_2 = 1.0 / 10.0;

  const Float q0       = (-1.0 * um1 + 5.0 * u0 + 2.0 * up1) / 6.0;
  const Float q1       = (2.0 * u0 + 5.0 * up1 - 1.0 * up2) / 6.0;
  const Float q2       = (11.0 * up1 - 7.0 * up2 + 2.0 * up3) / 6.0;

  const auto b0        = (13.0 / 12.0) * Igor::sqr(um1 - 2.0 * u0 + up1) +
                         (1.0 / 4.0) * Igor::sqr(um1 - 4.0 * u0 + 3.0 * up1);
  const auto b1        = (13.0 / 12.0) * Igor::sqr(u0 - 2.0 * up1 + up2) +  //
                         (1.0 / 4.0) * Igor::sqr(u0 - up2);
  const auto b2        = (13.0 / 12.0) * Igor::sqr(up1 - 2.0 * up2 + up3) +
                         (1.0 / 4.0) * Igor::sqr(3.0 * up1 - 4.0 * up2 + up3);

  const auto a0        = dR_0 / Igor::sqr(b0 + eps);
  const auto a1        = dR_1 / Igor::sqr(b1 + eps);
  const auto a2        = dR_2 / Igor::sqr(b2 + eps);
  const auto asum      = a0 + a1 + a2;

  return (a0 * q0 + a1 * q1 + a2 * q2) / asum;
}

template <typename Float, Layout LAYOUT>
constexpr void weno_reconstruction(const Grid<Float, LAYOUT>& grid,
                                   const Scalar<Float, LAYOUT> s,
                                   FaceVector<Float, LAYOUT> sL,
                                   FaceVector<Float, LAYOUT> sR) {
  IGOR_ASSERT(grid.nghost() >= 3, "WENO-5 scheme requires at least 3 ghost cells.");

  grid.template foreach_face_i<Dimension::X>(FOREACH_FUNC {
    sL.x(i, j) = weno_left_biased_reconstruction(
        s(i - 3, j), s(i - 2, j), s(i - 1, j), s(i, j), s(i + 1, j));
    sR.x(i, j) = weno_right_biased_reconstruction(
        s(i - 2, j), s(i - 1, j), s(i, j), s(i + 1, j), s(i + 2, j));
  });

  grid.template foreach_face_i<Dimension::Y>(FOREACH_FUNC {
    sL.y(i, j) = weno_left_biased_reconstruction(
        s(i, j - 3), s(i, j - 2), s(i, j - 1), s(i, j), s(i, j + 1));
    sR.y(i, j) = weno_right_biased_reconstruction(
        s(i, j - 2), s(i, j - 1), s(i, j), s(i, j + 1), s(i, j + 2));
  });
}
