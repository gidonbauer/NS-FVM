#pragma once

#include <Igor/Math.hpp>

#include "Grid.hpp"

namespace Cartesian {

template <typename Float, Layout LAYOUT>
constexpr void advection_calc_flux(const Grid<Float, LAYOUT>& grid,
                                   const FaceVector<Float, LAYOUT> u,
                                   const Scalar<Float, LAYOUT> s,
                                   const FaceVector<Float, LAYOUT> sL,
                                   const FaceVector<Float, LAYOUT> sR,
                                   Float D,
                                   FaceVector<Float, LAYOUT> F) {
  grid.template foreach_face_i<Dimension::X>(FOREACH_FUNC {
    // const auto si   = (s(i, j) + s(i - 1, j)) / 2.0;
    // const auto si = u.x(i, j) >= 0.0 ? s(i - 1, j) : s(i, j);
    const auto si   = u.x(i, j) >= 0.0 ? sL.x(i, j) : sR.x(i, j);
    const auto dsdx = (s(i, j) - s(i - 1, j)) / grid.dx();
    F.x(i, j)       = -si * u.x(i, j) + D * dsdx;
  });

  grid.template foreach_face_i<Dimension::Y>(FOREACH_FUNC {
    // const auto si   = (s(i, j) + s(i, j - 1)) / 2.0;
    // const auto si = u.y(i, j) >= 0.0 ? s(i, j - 1) : s(i, j);
    const auto si   = u.y(i, j) >= 0.0 ? sL.y(i, j) : sR.y(i, j);
    const auto dsdy = (s(i, j) - s(i, j - 1)) / grid.dy();
    F.y(i, j)       = -si * u.y(i, j) + D * dsdy;
  });
}

template <typename Float, Layout LAYOUT>
constexpr void advection_update_s(const Grid<Float, LAYOUT>& grid,
                                  Float dt,
                                  const FaceVector<Float, LAYOUT> F,
                                  const Scalar<Float, LAYOUT> s_old,
                                  Scalar<Float, LAYOUT> s) {
  grid.foreach_i(FOREACH_FUNC {
    s(i, j) = s_old(i, j) + dt * ((F.right(i, j) - F.left(i, j)) / grid.dx() +
                                  (F.top(i, j) - F.bottom(i, j)) / grid.dy());
  });
}

}  // namespace Cartesian

// = WENO-5 ========================================================================================
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
// = WENO-5 ========================================================================================

template <typename Float, Layout LAYOUT>
constexpr void advection_calc_flux(const Grid<Float, LAYOUT>& grid,
                                   const FaceVector<Float, LAYOUT> u,
                                   const Scalar<Float, LAYOUT> s,
                                   Float D,
                                   FaceVector<Float, LAYOUT> F) {
  static auto sL = grid.alloc_face_vector();
  static auto sR = grid.alloc_face_vector();
  weno_reconstruction(grid, s, sL, sR);

  switch (grid.coords()) {
    case Coordinates::CARTESIAN: return Cartesian::advection_calc_flux(grid, u, s, sL, sR, D, F);
    case Coordinates::POLAR:     Igor::Todo("Calculate advection flux in polar coordinates.");
  }
  Igor::Panic("Unreachable");
}

template <typename Float, Layout LAYOUT>
constexpr void advection_update_s(const Grid<Float, LAYOUT>& grid,
                                  Float dt,
                                  const FaceVector<Float, LAYOUT> F,
                                  const Scalar<Float, LAYOUT> s_old,
                                  Scalar<Float, LAYOUT> s) {
  switch (grid.coords()) {
    case Coordinates::CARTESIAN: return Cartesian::advection_update_s(grid, dt, F, s_old, s);
    case Coordinates::POLAR:     Igor::Todo("Update species in polar coordinates.");
  }
  Igor::Panic("Unreachable");
}

template <typename Float, Layout LAYOUT>
constexpr auto advection_adjust_dt(const Grid<Float, LAYOUT>& grid, Float D, Float CFL) -> Float {
  const auto h = std::min(grid.dx(), grid.dy());
  return CFL * Igor::sqr(h) / D;
}
