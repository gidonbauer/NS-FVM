#pragma once

/* * * * * * * * * * * *\
*      x => theta       *
*      y => r           *
*      u => u_theta     *
*      v => u_r         *
\* * * * * * * * * * * */

#include <Igor/Math.hpp>

#include "Grid.hpp"

namespace Polar {

// =================================================================================================
template <typename Float, Layout LAYOUT>
constexpr void calc_div(const Grid<Float, LAYOUT>& grid,
                        const FaceVector<Float, LAYOUT> uf,
                        Scalar<Float, LAYOUT> div) {
  grid.foreach_i(FOREACH_FUNC {
    const auto duthdth = (uf.right(i, j) - uf.left(i, j)) / grid.dx();
    const auto durdr   = (uf.top(i, j) - uf.bottom(i, j)) / grid.dy();
    const auto ur      = (uf.top(i, j) + uf.bottom(i, j)) / 2.0;
    const auto r       = grid.ym(j);
    div(i, j)          = durdr + duthdth / r + ur / r;
  });
}

// =================================================================================================
template <typename Float, Layout LAYOUT>
constexpr void calc_mom_flux(const Grid<Float, LAYOUT>& grid,
                             const FaceVector<Float, LAYOUT> u,
                             const Scalar<Float, LAYOUT> p,
                             Float rho,
                             Float mu,
                             Scalar<Float, LAYOUT> FUX,
                             VertexScalar<Float, LAYOUT> FUY,
                             VertexScalar<Float, LAYOUT> FVX,
                             Scalar<Float, LAYOUT> FVY) {
  const auto nu = mu / rho;
  grid.foreach_a(FOREACH_FUNC {
    const auto uth     = (u.right(i, j) + u.left(i, j)) / 2.0;
    const auto ur      = (u.top(i, j) + u.bottom(i, j)) / 2.0;
    const auto duthdth = (u.right(i, j) - u.left(i, j)) / grid.dx();
    const auto durdr   = (u.top(i, j) - u.bottom(i, j)) / grid.dy();
    const auto r       = grid.ym(j);

    FUX(i, j)          = -Igor::sqr(uth) - p(i, j) / rho + 2.0 * nu * (duthdth + ur) / r;
    FVY(i, j)          = -Igor::sqr(ur) - p(i, j) / rho + 2.0 * nu * durdr;
  });

  grid.foreach_vertex_i(FOREACH_FUNC {
    const auto uth    = (u.x(i, j) + u.x(i, j - 1)) / 2.0;
    const auto ur     = (u.y(i, j) + u.y(i - 1, j)) / 2.0;
    const auto duthdr = (u.x(i, j) - u.x(i, j - 1)) / grid.dy();
    const auto durdth = (u.y(i, j) - u.y(i - 1, j)) / grid.dx();
    const auto r      = grid.y(j);

    FUY(i, j)         = -uth * ur + nu * (duthdr + durdth / r - uth / r);
    FVX(i, j)         = -uth * ur + nu * (duthdr + durdth / r - uth / r);
  });
}

// =================================================================================================
template <typename Float, Layout LAYOUT>
constexpr void update_u(const Grid<Float, LAYOUT>& grid,
                        Float dt,
                        const Scalar<Float, LAYOUT> FUX,
                        const VertexScalar<Float, LAYOUT> FUY,
                        const VertexScalar<Float, LAYOUT> FVX,
                        const Scalar<Float, LAYOUT> FVY,
                        const FaceVector<Float, LAYOUT> u_old,
                        FaceVector<Float, LAYOUT> u) {
  grid.template foreach_face_i<Dimension::X>(FOREACH_FUNC {
    const auto dTththdth = (FUX(i, j) - FUX(i - 1, j)) / grid.dx();
    const auto dTrthdr   = (FUY(i, j + 1) - FUY(i, j)) / grid.dy();
    const auto Trth      = (FUY(i, j + 1) + FUY(i, j)) / 2.0;
    const auto Tthr      = Trth;
    const auto r         = grid.ym(j);

    u.x(i, j)            = u_old.x(i, j) + dt * (dTththdth / r + dTrthdr + (Trth + Tthr) / r);
  });

  grid.template foreach_face_i<Dimension::Y>(FOREACH_FUNC {
    const auto dTthrdth = (FVX(i + 1, j) - FVX(i, j)) / grid.dx();
    const auto dTrrdr   = (FVY(i, j) - FVY(i, j - 1)) / grid.dy();
    const auto Tthth    = (FUX(i, j) + FUX(i, j - 1)) / 2.0;
    const auto Trr      = (FVY(i, j) + FVY(i, j - 1)) / 2.0;
    const auto r        = grid.y(j);

    u.y(i, j)           = u_old.y(i, j) + dt * (dTthrdth / r + dTrrdr + (Trr - Tthth) / r);
  });
}

// =================================================================================================
template <typename Float, Layout LAYOUT>
constexpr void correct_velocity(const Grid<Float, LAYOUT>& grid,
                                const Scalar<Float, LAYOUT> dp,
                                Float rho,
                                Float dt,
                                FaceVector<Float, LAYOUT> u,
                                Scalar<Float, LAYOUT> p) {
  grid.foreach_a(FOREACH_FUNC { p(i, j) += dp(i, j); });
  grid.template foreach_face_i<Dimension::X>(FOREACH_FUNC {
    const auto r  = grid.ym(j);
    u.x(i, j)    -= (dt / rho) * (dp(i, j) - dp(i - 1, j)) / (r * grid.dx());
  });
  grid.template foreach_face_i<Dimension::Y>(
      FOREACH_FUNC { u.y(i, j) -= (dt / rho) * (dp(i, j) - dp(i, j - 1)) / grid.dy(); });
}

}  // namespace Polar

namespace ALEPolar {

using Polar::calc_div;
using Polar::correct_velocity;
using Polar::update_u;

// =================================================================================================
// We assume that the ALE mesh movement is only in the axis directions, meaning only in r- or theta-
// direction.
// Additionally, we assmue that the grid velocity is the same everywhere. This means that the
// Jacobian matrix is the identity matrix and its determinant is one. The flux then almost identical
// to the non-ALE case, the only difference are that the radius `r` changes, this is done in the
// `Grid` class, and the addition of the grid veclocity `w`. The fluxes `FUY` and `FVX` are then no
// longer the same.
template <typename Float, Layout LAYOUT>
constexpr void calc_mom_flux(const Grid<Float, LAYOUT>& grid,
                             const FaceVector<Float, LAYOUT> u,
                             const Scalar<Float, LAYOUT> p,
                             Float rho,
                             Float mu,
                             Vec2<Float> w,
                             Scalar<Float, LAYOUT> FUX,
                             VertexScalar<Float, LAYOUT> FUY,
                             VertexScalar<Float, LAYOUT> FVX,
                             Scalar<Float, LAYOUT> FVY) {
  const auto nu = mu / rho;
  grid.foreach_a(FOREACH_FUNC {
    const auto uth     = (u.right(i, j) + u.left(i, j)) / 2.0;
    const auto ur      = (u.top(i, j) + u.bottom(i, j)) / 2.0;
    const auto duthdth = (u.right(i, j) - u.left(i, j)) / grid.dx();
    const auto durdr   = (u.top(i, j) - u.bottom(i, j)) / grid.dy();
    const auto r       = grid.ym(j);

    FUX(i, j) = -Igor::sqr(uth) - p(i, j) / rho + 2.0 * nu * (duthdth + ur) / r + uth * w.theta();
    FVY(i, j) = -Igor::sqr(ur) - p(i, j) / rho + 2.0 * nu * durdr + ur * w.r();
  });

  grid.foreach_vertex_i(FOREACH_FUNC {
    const auto uth    = (u.x(i, j) + u.x(i, j - 1)) / 2.0;
    const auto ur     = (u.y(i, j) + u.y(i - 1, j)) / 2.0;
    const auto duthdr = (u.x(i, j) - u.x(i, j - 1)) / grid.dy();
    const auto durdth = (u.y(i, j) - u.y(i - 1, j)) / grid.dx();
    const auto r      = grid.y(j);

    FUY(i, j)         = -uth * ur + nu * (duthdr + durdth / r - uth / r) + ur * w.r();
    FVX(i, j)         = -uth * ur + nu * (duthdr + durdth / r - uth / r) + uth * w.theta();
  });
}

}  // namespace ALEPolar
