#include <numbers>

#include <Igor/Defer.hpp>
#include <Igor/Logging.hpp>

#include "Grid.hpp"
#include "VTKWriter.hpp"

using Float           = double;
constexpr Index NX    = 128;
constexpr Index NY    = 128;

constexpr Float x_min = 0.0;
constexpr Float x_max = 1.0;
constexpr Float y_min = 0.0;
constexpr Float y_max = 1.0;

constexpr Float rho   = 1.0;
constexpr Float mu    = 1.0;

auto F(Float t) -> Float { return std::exp(-2.0 * mu / rho * t); }
auto u_analytical(Float x, Float y, Float t) -> Float {
  constexpr auto pi = std::numbers::pi_v<Float>;
  return std::sin(2.0 * pi * x) * std::cos(2.0 * pi * y) * F(t);
}
auto v_analytical(Float x, Float y, Float t) -> Float {
  constexpr auto pi = std::numbers::pi_v<Float>;
  return -std::cos(2.0 * pi * x) * std::sin(2.0 * pi * y) * F(t);
}

auto main() -> int {
  Grid grid(x_min, x_max, NX, y_min, y_max, NY);

  auto p = grid.alloc_scalar();
  Igor::Defer free_p([&] { grid.free(p); });

  auto u = grid.alloc_vector();
  Igor::Defer free_u([&] { grid.free(u); });

  const Float dx = (x_max - x_min) / NX;
  const Float dy = (y_max - y_min) / NY;
  for (Index i = 0; i < NX; ++i) {
    for (Index j = 0; j < NY; ++j) {
      const Float x = x_min + (i + 0.5) * dx;
      const Float y = y_min + (j + 0.5) * dy;
      u.x(i, j)     = u_analytical(x, y, 0.0);
      u.y(i, j)     = v_analytical(x, y, 0.0);
      p(i, j)       = x < 0.5 ? 1.0 : 0.0;
    }
  }

  VTKWriter writer("output/", grid);
  writer.add_field("p", p);
  writer.add_field("u", u);
  if (!writer.write(0.0)) { return 1; }

  Igor::Info("Ok.");
}
