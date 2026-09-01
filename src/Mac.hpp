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
constexpr void interpolate(const Grid<Float, LAYOUT>& grid,
                           const FaceVector<Float, LAYOUT> uf,
                           Vector<Float, LAYOUT> ui) {
  switch (grid.coords()) {
    case Coordinates::CARTESIAN: return Cartesian::interpolate(grid, uf, ui);
    case Coordinates::POLAR:     return Polar::interpolate(grid, uf, ui);
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
  Float u_max = 0.0;
  grid.template foreach_face_i<Dimension::X, Exec::SERIAL>(
      [=, &u_max](Index i, Index j) { u_max = std::max(std::abs(u.x(i, j)), u_max); });
  grid.template foreach_face_i<Dimension::Y, Exec::SERIAL>(
      [=, &u_max](Index i, Index j) { u_max = std::max(std::abs(u.y(i, j)), u_max); });
  const auto h = std::min(grid.dx(), grid.dy());
  return std::min({
      CFL * h / u_max,
      CFL * 0.25 * Igor::sqr(h) * rho / mu,
  });
}
