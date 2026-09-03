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

template <typename Float, Layout LAYOUT>
void print(const FaceVector<Float, LAYOUT> u) {
  std::println("u.x =");
  for (Index j = -u.x.nghost(); j < u.x.ny() + u.x.nghost(); ++j) {
    for (Index i = -u.x.nghost(); i < u.x.nx() + u.x.nghost(); ++i) {
      if ((i < 0 || i >= u.x.nx()) || (j < 0 || j >= u.x.ny())) {
        std::print("\033[94m{: .2f}\033[0m\t", u.x(i, j));
      } else {
        std::print("\033[32m{: .2f}\033[0m\t", u.x(i, j));
      }
    }
    std::cout << '\n';
  }
  std::cout << '\n';

  std::println("u.y =");
  for (Index j = -u.y.nghost(); j < u.y.ny() + u.y.nghost(); ++j) {
    for (Index i = -u.y.nghost(); i < u.y.nx() + u.y.nghost(); ++i) {
      if ((i < 0 || i >= u.y.nx()) || (j < 0 || j >= u.y.ny())) {
        std::print("\033[94m{: .2f}\033[0m\t", u.y(i, j));
      } else {
        std::print("\033[32m{: .2f}\033[0m\t", u.y(i, j));
      }
    }
    std::cout << '\n';
  }
}

auto test_dirichlet() -> bool {
  bool any_failed        = false;

  constexpr Index NGHOST = 3;
  Grid<double> grid(0.0, 1.0, 2, 0.0, 1.0, 2, NGHOST);
  auto u = grid.alloc_face_vector();
  fill(u.x, std::numeric_limits<double>::quiet_NaN());
  fill(u.y, std::numeric_limits<double>::quiet_NaN());
  grid.foreach_face_i<Dimension::X>(FOREACH_FUNC { u.x(i, j) = i + u.x.nx() * j; });
  grid.foreach_face_i<Dimension::Y>(FOREACH_FUNC { u.y(i, j) = i + u.y.nx() * j; });
  print(u);
  std::cout << "\n============================================================\n";

  Dirichlet2::apply_left_align<double>(grid, u.x, true, 0.0, 1.0);
  Dirichlet2::apply_right_align<double>(grid, u.x, true, 0.0, 1.0);
  Dirichlet2::apply_bottom_offset<double>(grid, u.x, false, 0.0, 1.0);
  Dirichlet2::apply_top_offset<double>(grid, u.x, false, 0.0, 1.0);

  Dirichlet2::apply_left_offset<double>(grid, u.y, false, 0.0, 1.0);
  Dirichlet2::apply_right_offset<double>(grid, u.y, false, 0.0, 1.0);
  Dirichlet2::apply_bottom_align<double>(grid, u.y, true, 0.0, 1.0);
  Dirichlet2::apply_top_align<double>(grid, u.y, true, 0.0, 1.0);
  print(u);

  return !any_failed;
}

auto test_neumann() -> bool {
  bool any_failed        = false;

  constexpr Index NGHOST = 2;
  Grid<double> grid(0.0, 1.0, 2, 0.0, 1.0, 2, NGHOST);
  auto u = grid.alloc_face_vector();
  fill(u.x, std::numeric_limits<double>::quiet_NaN());
  fill(u.y, std::numeric_limits<double>::quiet_NaN());
  grid.foreach_face_i<Dimension::X>(FOREACH_FUNC { u.x(i, j) = i + u.x.nx() * j; });
  grid.foreach_face_i<Dimension::Y>(FOREACH_FUNC { u.y(i, j) = i + u.y.nx() * j; });
  print(u);
  std::cout << "\n============================================================\n";

  Neumann2::apply_left_align<double>(grid, u.x);
  Neumann2::apply_right_align<double>(grid, u.x);
  Neumann2::apply_bottom_offset<double>(grid, u.x);
  Neumann2::apply_top_offset<double>(grid, u.x);

  Neumann2::apply_left_offset<double>(grid, u.y);
  Neumann2::apply_right_offset<double>(grid, u.y);
  Neumann2::apply_bottom_align<double>(grid, u.y);
  Neumann2::apply_top_align<double>(grid, u.y);
  print(u);

  return !any_failed;
}

auto main() -> int {
  bool any_failed = false;
  // if (!test_dirichlet()) {
  //   Igor::Error("Dirichlet boundart conditions failed.");
  //   any_failed = true;
  // }

  if (!test_neumann()) {
    Igor::Error("Neumann boundart conditions failed.");
    any_failed = true;
  }

  return any_failed ? 1 : 0;
}
