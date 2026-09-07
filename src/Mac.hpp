#pragma once

#include <Igor/Math.hpp>

#include "Grid.hpp"
#include "MacCartesian.hpp"
#include "MacPolar.hpp"

// =================================================================================================
template <typename Float, Layout LAYOUT>
constexpr void calc_div(const Grid<Float, LAYOUT>& grid,
                        const FaceVector<Float, LAYOUT> uf,
                        Scalar<Float, LAYOUT> div) {
  switch (grid.coords()) {
    case Coordinates::CARTESIAN: return Cartesian::calc_div(grid, uf, div);
    case Coordinates::POLAR:     return Polar::calc_div(grid, uf, div);
  }
  Igor::Panic("Unreachable");
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
  switch (grid.coords()) {
    case Coordinates::CARTESIAN:
      return Cartesian::calc_flux(grid, u, p, rho, mu, FUX, FUY, FVX, FVY);
    case Coordinates::POLAR: return Polar::calc_flux(grid, u, p, rho, mu, FUX, FUY, FVX, FVY);
  }
  Igor::Panic("Unreachable");
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
  switch (grid.coords()) {
    case Coordinates::CARTESIAN: return Cartesian::update_u(grid, dt, FUX, FUY, FVX, FVY, u_old, u);
    case Coordinates::POLAR:     return Polar::update_u(grid, dt, FUX, FUY, FVX, FVY, u_old, u);
  }
  Igor::Panic("Unreachable");
}

// =================================================================================================
template <typename Float, Layout LAYOUT>
constexpr void correct_velocity(const Grid<Float, LAYOUT>& grid,
                                const Scalar<Float, LAYOUT> dp,
                                Float rho,
                                Float dt,
                                FaceVector<Float, LAYOUT> u,
                                Scalar<Float, LAYOUT> p) {
  switch (grid.coords()) {
    case Coordinates::CARTESIAN: return Cartesian::correct_velocity(grid, dp, rho, dt, u, p);
    case Coordinates::POLAR:     return Polar::correct_velocity(grid, dp, rho, dt, u, p);
  }
  Igor::Panic("Unreachable");
}

// =================================================================================================
template <typename Float, Layout LAYOUT>
constexpr auto adjust_dt(const Grid<Float, LAYOUT>& grid,
                         const FaceVector<Float, LAYOUT> u,
                         Float rho,
                         Float mu,
                         Float CFL) noexcept -> Float {
  Float ux_max = 0.0;
  Float uy_max = 0.0;
  grid.template foreach_face_i<Dimension::X, Exec::SERIAL>(
      [=, &ux_max](Index i, Index j) { ux_max = std::max(std::abs(u.x(i, j)), ux_max); });
  grid.template foreach_face_i<Dimension::Y, Exec::SERIAL>(
      [=, &uy_max](Index i, Index j) { uy_max = std::max(std::abs(u.y(i, j)), uy_max); });

  // Correction for polar coordinates
  const auto hx = grid.coords() == Coordinates::POLAR ? grid.ym(0) * grid.dx() : grid.dx();
  const auto hy = grid.dy();

  // Advection: dt * (|u|/hx + |v|/hy) <= CFL
  const auto adv = ux_max / hx + uy_max / hy;
  // Diffusion: dt * 2 * nu * (1/hx^2 + 1/hy^2) <= CFL
  const auto diff         = 2.0 * (mu / rho) * (1.0 / Igor::sqr(hx) + 1.0 / Igor::sqr(hy));

  constexpr auto no_limit = std::numeric_limits<Float>::max();
  return std::min(adv > 0.0 ? CFL / adv : no_limit,  //
                  diff > 0.0 ? CFL / diff : no_limit);
}
