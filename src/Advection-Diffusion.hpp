#pragma once

#include <Igor/Math.hpp>

#include "Grid.hpp"
#include "WENO5.hpp"

// =================================================================================================
// = Cartesian =====================================================================================
// =================================================================================================
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
    const auto si   = u.x(i, j) >= 0.0 ? sL.x(i, j) : sR.x(i, j);
    const auto dsdx = (s(i, j) - s(i - 1, j)) / grid.dx();
    F.x(i, j)       = -si * u.x(i, j) + D * dsdx;
  });

  grid.template foreach_face_i<Dimension::Y>(FOREACH_FUNC {
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

template <typename Float, Layout LAYOUT>
constexpr void advection_update_s(const Grid<Float, LAYOUT>& grid,
                                  Float dt,
                                  const FaceVector<Float, LAYOUT> F,
                                  const Scalar<Float, LAYOUT> src,
                                  const Scalar<Float, LAYOUT> s_old,
                                  Scalar<Float, LAYOUT> s) {
  grid.foreach_i(FOREACH_FUNC {
    s(i, j) = s_old(i, j) + dt * ((F.right(i, j) - F.left(i, j)) / grid.dx() +
                                  (F.top(i, j) - F.bottom(i, j)) / grid.dy() + src(i, j));
  });
}

}  // namespace Cartesian
// =================================================================================================
// = Cartesian =====================================================================================
// =================================================================================================

// =================================================================================================
// = Polar =========================================================================================
// =================================================================================================
namespace Polar {

template <typename Float, Layout LAYOUT>
constexpr void advection_calc_flux(const Grid<Float, LAYOUT>& grid,
                                   const FaceVector<Float, LAYOUT> u,
                                   const Scalar<Float, LAYOUT> s,
                                   const FaceVector<Float, LAYOUT> sL,
                                   const FaceVector<Float, LAYOUT> sR,
                                   Float D,
                                   FaceVector<Float, LAYOUT> F) {
  grid.template foreach_face_i<Dimension::X>(FOREACH_FUNC {
    const auto si   = u.x(i, j) >= 0.0 ? sL.x(i, j) : sR.x(i, j);
    const auto dsdx = (s(i, j) - s(i - 1, j)) / grid.dx();
    const auto r    = grid.ym(j);
    F.x(i, j)       = -si * u.x(i, j) + D * dsdx / r;
  });

  grid.template foreach_face_i<Dimension::Y>(FOREACH_FUNC {
    const auto si   = u.y(i, j) >= 0.0 ? sL.y(i, j) : sR.y(i, j);
    const auto dsdy = (s(i, j) - s(i, j - 1)) / grid.dy();
    F.y(i, j)       = -si * u.y(i, j) + D * dsdy;
  });
}

template <typename Float, Layout LAYOUT>
constexpr void advection_update_s(const Grid<Float, LAYOUT>& grid,
                                  Float dt,
                                  const FaceVector<Float, LAYOUT> F,
                                  const Scalar<Float, LAYOUT> src,
                                  const Scalar<Float, LAYOUT> s_old,
                                  Scalar<Float, LAYOUT> s) {
  grid.foreach_i(FOREACH_FUNC {
    const auto dFthdth = (F.right(i, j) - F.left(i, j)) / grid.dx();
    const auto dFrdr   = (F.top(i, j) - F.bottom(i, j)) / grid.dy();
    const auto Fr      = (F.top(i, j) + F.bottom(i, j)) / 2.0;
    const auto r       = grid.ym(j);
    s(i, j)            = s_old(i, j) + dt * (dFrdr + dFthdth / r + Fr / r + src(i, j));
  });
}

template <typename Float, Layout LAYOUT>
constexpr void advection_update_s(const Grid<Float, LAYOUT>& grid,
                                  Float dt,
                                  const FaceVector<Float, LAYOUT> F,
                                  const Scalar<Float, LAYOUT> s_old,
                                  Scalar<Float, LAYOUT> s) {
  grid.foreach_i(FOREACH_FUNC {
    const auto dFthdth = (F.right(i, j) - F.left(i, j)) / grid.dx();
    const auto dFrdr   = (F.top(i, j) - F.bottom(i, j)) / grid.dy();
    const auto Fr      = (F.top(i, j) + F.bottom(i, j)) / 2.0;
    const auto r       = grid.ym(j);
    s(i, j)            = s_old(i, j) + dt * (dFrdr + dFthdth / r + Fr / r);
  });
}

}  // namespace Polar
// =================================================================================================
// = Polar =========================================================================================
// =================================================================================================

template <typename Float, Layout LAYOUT>
constexpr auto advection_adjust_dt(const Grid<Float, LAYOUT>& grid, Float D, Float CFL) -> Float {
  const auto h = std::min(grid.dx(), grid.dy());
  return CFL * 0.25 * Igor::sqr(h) / D;
}

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
    case Coordinates::POLAR:     return Polar::advection_calc_flux(grid, u, s, sL, sR, D, F);
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
    case Coordinates::POLAR:     return Polar::advection_update_s(grid, dt, F, s_old, s);
  }
  Igor::Panic("Unreachable");
}

template <typename Float, Layout LAYOUT>
constexpr void advection_update_s(const Grid<Float, LAYOUT>& grid,
                                  Float dt,
                                  const FaceVector<Float, LAYOUT> F,
                                  const Scalar<Float, LAYOUT> src,
                                  const Scalar<Float, LAYOUT> s_old,
                                  Scalar<Float, LAYOUT> s) {
  switch (grid.coords()) {
    case Coordinates::CARTESIAN: return Cartesian::advection_update_s(grid, dt, F, src, s_old, s);
    case Coordinates::POLAR:     return Polar::advection_update_s(grid, dt, F, src, s_old, s);
  }
  Igor::Panic("Unreachable");
}
