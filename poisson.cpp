#include <charconv>
#include <numbers>

#include <poisfft.h>

#include <Igor/Defer.hpp>
#include <Igor/Logging.hpp>

#include "Grid.hpp"
#include "IO.hpp"
#include "VTKWriter.hpp"

using Float = double;

auto F(Float t) -> Float { return std::exp(-2.0 * t); }
auto u_analytical(Float x, Float y, Float t) -> Float {
  constexpr auto pi = std::numbers::pi_v<Float>;
  return std::sin(2.0 * pi * x) * std::cos(2.0 * pi * y) * F(t) * x;
}
auto v_analytical(Float x, Float y, Float t) -> Float {
  constexpr auto pi = std::numbers::pi_v<Float>;
  return -std::cos(2.0 * pi * x) * std::sin(2.0 * pi * y) * F(t);
}

template <typename Float, Layout LAYOUT>
void calc_div(const Grid<Float, LAYOUT>& grid,
              FaceVector<Float, LAYOUT> u,
              Scalar<Float, LAYOUT> div) {
  grid.foreach_i(FOREACH_FUNC {
    div(i, j) = (u.right(i, j) - u.left(i, j)) / grid.dx() +  //
                (u.top(i, j) - u.bottom(i, j)) / grid.dy();
  });
}

template <typename Float, Layout LAYOUT>
void interpolate(const Grid<Float, LAYOUT>& grid,
                 FaceVector<Float, LAYOUT> u,
                 Vector<Float, LAYOUT> ui) {
  grid.foreach_i(FOREACH_FUNC {
    ui.x(i, j) = (u.right(i, j) + u.left(i, j)) / 2.0;
    ui.y(i, j) = (u.top(i, j) + u.bottom(i, j)) / 2.0;
  });
}

template <typename Float, Layout LAYOUT>
void periodic_bconds(FaceVector<Float, LAYOUT> u) {
  IGOR_ASSERT(u.x.nghost() == 1, "Expected exactly one ghost cell but got {}", u.x.nghost());
  for (Index i = 0; i < u.x.nx(); ++i) {
    u.x(i, -1)           = u.x(i, u.x.ny() - 2);
    u.x(i, 0)            = u.x(i, u.x.ny() - 1);
    u.x(i, u.x.ny() - 1) = u.x(i, 0);
    u.x(i, u.x.ny())     = u.x(i, 1);
  }
  for (Index j = 0; j < u.x.ny(); ++j) {
    u.x(-1, j)           = u.x(u.x.nx() - 2, j);
    u.x(0, j)            = u.x(u.x.nx() - 1, j);
    u.x(u.x.nx() - 1, j) = u.x(0, j);
    u.x(u.x.nx(), j)     = u.x(1, j);
  }

  for (Index i = 0; i < u.y.nx(); ++i) {
    u.y(i, -1)           = u.y(i, u.y.ny() - 2);
    u.y(i, 0)            = u.y(i, u.y.ny() - 1);
    u.y(i, u.y.ny() - 1) = u.y(i, 0);
    u.y(i, u.y.ny())     = u.y(i, 1);
  }
  for (Index j = 0; j < u.y.ny(); ++j) {
    u.y(-1, j)           = u.y(u.y.nx() - 2, j);
    u.y(0, j)            = u.y(u.y.nx() - 1, j);
    u.y(u.y.nx() - 1, j) = u.y(0, j);
    u.y(u.y.nx(), j)     = u.y(1, j);
  }
}

auto main(int argc, char** argv) -> int {
  const auto usage_str = Igor::detail::format("Usage: {} <grid size N>", argv[0]);
  if (argc < 2) {
    Igor::Error("{}", usage_str);
    return 1;
  }

  Index N = 0;
  if (std::from_chars(argv[1], argv[1] + std::strlen(argv[1]), N).ec != std::errc{} || N <= 0) {
    Igor::Error("{}", usage_str);
    Igor::Error("  Invalid grid size N `{}`", argv[1]);
    return 1;
  }
  const auto output_dir = get_output_directory() + "/" + std::to_string(N);
  if (!init_output_directory(output_dir)) { return 1; }

  const Float L = 1.0;
  Grid<Float> grid(0.0, L, N, 0.0, L, N);

  auto p = grid.alloc_scalar();
  IGOR_DEFER(grid.free(p););

  auto u = grid.alloc_face_vector();
  IGOR_DEFER(grid.free(u););

  auto ui = grid.alloc_vector();
  IGOR_DEFER(grid.free(ui););

  auto div = grid.alloc_scalar();
  IGOR_DEFER(grid.free(div););

  auto rhs = grid.alloc_scalar();
  IGOR_DEFER(grid.free(rhs););

  VTKWriter writer(output_dir, grid);
  writer.add_field("p", p);
  writer.add_field("ui", ui);
  writer.add_field("div", div);
  writer.add_field("rhs", rhs);

  grid.foreach_i(FOREACH_FUNC {
    u.x(i, j) = u_analytical(grid.x(i), grid.ym(j), 0.0);
    u.y(i, j) = v_analytical(grid.xm(i), grid.y(j), 0.0);
  });
  periodic_bconds(u);
  interpolate(grid, u, ui);

  const Float dt = 1e-1;
  calc_div(grid, u, div);
  grid.foreach_i(FOREACH_FUNC { rhs(i, j) = div(i, j) / dt; });

  writer.write(0.0);

  const std::array<int, 2> ns   = {N, N};
  const std::array<Float, 2> Ls = {L, L};
  const std::array<int, 4> BCs  = {
      PoisFFT::NEUMANN, PoisFFT::NEUMANN, PoisFFT::NEUMANN, PoisFFT::NEUMANN};
  PoisFFT::Solver<2, Float> solver(ns.data(), Ls.data(), BCs.data(), PoisFFT::FINITE_DIFFERENCE_2);

  const std::array<int, 4> ngs = {grid.nghost(), grid.nghost(), grid.nghost(), grid.nghost()};
  solver.execute(p.data(), rhs.data(), ngs.data(), ngs.data());

  grid.foreach_face_i<Dimension::X>(FOREACH_FUNC {
    const auto dpdx  = (p(i, j) - p(i - 1, j)) / grid.dx();
    u.x(i, j)       -= dt * dpdx;
  });
  grid.foreach_face_i<Dimension::Y>(FOREACH_FUNC {
    const auto dpdy  = (p(i, j) - p(i, j - 1)) / grid.dy();
    u.y(i, j)       -= dt * dpdy;
  });
  periodic_bconds(u);
  interpolate(grid, u, ui);

  calc_div(grid, u, div);

  writer.write(1.0);
}
