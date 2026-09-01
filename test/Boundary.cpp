#include <cmath>
#include <limits>
#include <numbers>

#include "BoundaryConditions.hpp"

constexpr auto pi = std::numbers::pi_v<double>;

template <typename Float>
auto U(Float x, Float y) -> Float {
  return std::sin(2.0 * pi * x) * std::cos(2.0 * pi * y);
}

template <typename Float>
auto V(Float x, Float y) -> Float {
  return -std::cos(2.0 * pi * x) * std::sin(2.0 * pi * y);
}

template <typename Float, Layout LAYOUT>
void fill_velocity(const Grid<Float, LAYOUT>& grid, FaceVector<Float, LAYOUT> u) {
  grid.template foreach_face_i<Dimension::X>(
      FOREACH_FUNC { u.x(i, j) = U(grid.x(i), grid.ym(j)); });
  grid.template foreach_face_i<Dimension::Y>(
      FOREACH_FUNC { u.y(i, j) = V(grid.xm(i), grid.y(j)); });
}

template <typename Float>
constexpr auto
approx_eq(Float lhs, Float rhs, Float tol = std::sqrt(std::numeric_limits<Float>::epsilon())) {
  return std::abs(lhs - rhs) <= tol;
}

auto test_dirichlet() -> bool {
  bool any_failed = false;

  // = NGHOST: 1 ===================================================================================
  {
    constexpr Index NGHOST = 1;
    Grid<double> grid(0.0, 1.0, 64, 0.0, 1.0, 64, NGHOST);
    auto u = grid.alloc_face_vector();
    fill_velocity(grid, u);

    VelocityBConds<double> bconds{
        .left   = Dirichlet<double>{.U = 1.0, .V = 1.0},
        .right  = Dirichlet<double>{.U = 2.0, .V = 2.0},
        .bottom = Dirichlet<double>{.U = 3.0, .V = 3.0},
        .top    = Dirichlet<double>{.U = 4.0, .V = 4.0},
    };
    apply_velocity_bconds(grid, bconds, u);

    // Left: U
    grid.foreach_range(0, 1, 0, u.x.ny(), [=, &any_failed](Index i, Index j) {
      constexpr double U = 1.0;
      if (u.x(0, j) != U) {
        Igor::Error("NGHOST={}: ({}, {}): Expected u.x={} but got {}", NGHOST, i, j, U, u.x(0, j));
        any_failed = true;
      }

      // const auto ui = (u.x(1, j) + u.x(-1, j)) / 2.0;
      // if (!approx_eq(ui, U)) {
      //   Igor::Error("NGHOST={}: ({}, {}): Expected ui={} but got {}", NGHOST, i, j, U, ui);
      //   any_failed = true;
      // }
    });

    // Left: V
    grid.foreach_range(0, 1, 0, u.y.ny(), [=, &any_failed](Index i, Index j) {
      constexpr double V = 1.0;
      const auto vi      = (u.y(0, j) + u.y(-1, j)) / 2.0;
      if (!approx_eq(vi, V)) {
        Igor::Error("NGHOST={}: ({}, {}): Expected u.y={} but got {}", NGHOST, i, j, V, vi);
        any_failed = true;
      }
    });
  }

  return !any_failed;
}

auto main() -> int {
  bool any_failed = false;
  if (!test_dirichlet()) {
    Igor::Error("Dirichlet boundart conditions failed.");
    any_failed = true;
  }

  return any_failed ? 1 : 0;
}
