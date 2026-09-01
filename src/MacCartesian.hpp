#pragma once

#include <Igor/Math.hpp>

#include "Grid.hpp"

namespace Cartesian {

// =================================================================================================
template <typename Float, Layout LAYOUT>
constexpr void calc_div(const Grid<Float, LAYOUT>& grid,
                        const FaceVector<Float, LAYOUT> uf,
                        Scalar<Float, LAYOUT> div) {
  grid.foreach_i(FOREACH_FUNC {
    div(i, j) = (uf.right(i, j) - uf.left(i, j)) / grid.dx() +  //
                (uf.top(i, j) - uf.bottom(i, j)) / grid.dy();
  });
}

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
template <typename Float, Layout LAYOUT>
constexpr void calc_flux(const Grid<Float, LAYOUT>& grid,
                         const FaceVector<Float, LAYOUT> u,
                         const Scalar<Float, LAYOUT> p,
                         Float rho,
                         Float mu,
                         Scalar<Float, LAYOUT> FUX,
                         VertexScalar<Float, LAYOUT> FUY,
                         VertexScalar<Float, LAYOUT> FVX,
                         Scalar<Float, LAYOUT> FVY) {
  grid.foreach_a(FOREACH_FUNC {
    const auto ui   = (u.right(i, j) + u.left(i, j)) / 2.0;
    const auto dudx = (u.right(i, j) - u.left(i, j)) / grid.dx();
    FUX(i, j)       = -Igor::sqr(ui) - p(i, j) / rho + 2.0 * mu / rho * dudx;

    const auto vi   = (u.top(i, j) + u.bottom(i, j)) / 2.0;
    const auto dvdy = (u.top(i, j) - u.bottom(i, j)) / grid.dy();
    FVY(i, j)       = -Igor::sqr(vi) - p(i, j) / rho + 2.0 * mu / rho * dvdy;
  });

  grid.foreach_vertex_i(FOREACH_FUNC {
    const auto ui   = (u.x(i, j) + u.x(i, j - 1)) / 2.0;
    const auto dudy = (u.x(i, j) - u.x(i, j - 1)) / grid.dy();

    const auto vi   = (u.y(i, j) + u.y(i - 1, j)) / 2.0;
    const auto dvdx = (u.y(i, j) - u.y(i - 1, j)) / grid.dx();

    FUY(i, j)       = -ui * vi + mu / rho * (dudy + dvdx);
    FVX(i, j)       = -ui * vi + mu / rho * (dudy + dvdx);
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
    u.x(i, j) = u_old.x(i, j) + dt * ((FUX(i, j) - FUX(i - 1, j)) / grid.dx() +
                                      (FUY(i, j + 1) - FUY(i, j)) / grid.dy());
  });

  grid.template foreach_face_i<Dimension::Y>(FOREACH_FUNC {
    u.y(i, j) = u_old.y(i, j) + dt * ((FVX(i + 1, j) - FVX(i, j)) / grid.dx() +
                                      (FVY(i, j) - FVY(i, j - 1)) / grid.dy());
  });
}

// =================================================================================================
template <typename Float, Layout LAYOUT>
constexpr void shift_dp_to_zero(const Grid<Float, LAYOUT>& grid, Scalar<Float, LAYOUT> dp) {
  Float avg_dp = 0.0;
  grid.template foreach_a<Exec::SERIAL>([=, &avg_dp](Index i, Index j) { avg_dp += dp(i, j); });
  avg_dp /= static_cast<Float>(grid.nx() * grid.ny());
  grid.foreach_a(FOREACH_FUNC { dp(i, j) -= avg_dp; });
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
  grid.template foreach_face_i<Dimension::X>(
      FOREACH_FUNC { u.x(i, j) -= (dt / rho) * (dp(i, j) - dp(i - 1, j)) / grid.dx(); });
  grid.template foreach_face_i<Dimension::Y>(
      FOREACH_FUNC { u.y(i, j) -= (dt / rho) * (dp(i, j) - dp(i, j - 1)) / grid.dy(); });
}

}  // namespace Cartesian
