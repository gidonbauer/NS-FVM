#include <charconv>
#include <random>

#include <poisfft.h>

#include <Igor/Defer.hpp>
#include <Igor/Logging.hpp>

#include "Grid.hpp"
#include "IO.hpp"
#include "VTKWriter.hpp"

using Float = double;

template <typename Float, Layout LAYOUT>
void calc_div(const Grid<Float, LAYOUT>& grid,
              FaceVector<Float, LAYOUT> uf,
              Scalar<Float, LAYOUT> div) {
  grid.foreach_i(FOREACH_FUNC {
    div(i, j) = (uf.right(i, j) - uf.left(i, j)) / grid.dx() +  //
                (uf.top(i, j) - uf.bottom(i, j)) / grid.dy();
  });
}

template <typename Float, Layout LAYOUT>
void interpolate(const Grid<Float, LAYOUT>& grid,
                 const Vector<Float, LAYOUT> u,
                 FaceVector<Float, LAYOUT> uf) {
  grid.template foreach_face_i<Dimension::X>(
      FOREACH_FUNC { uf.x(i, j) = (u.x(i, j) + u.x(i - 1, j)) / 2.0; });
  grid.template foreach_face_i<Dimension::Y>(
      FOREACH_FUNC { uf.y(i, j) = (u.y(i, j) + u.y(i, j - 1)) / 2.0; });
}

template <typename Float, Layout LAYOUT>
void interpolate(const Grid<Float, LAYOUT>& grid,
                 const FaceVector<Float, LAYOUT> uf,
                 Vector<Float, LAYOUT> u) {
  grid.foreach_i(FOREACH_FUNC {
    u.x(i, j) = (uf.right(i, j) + uf.left(i, j)) / 2.0;
    u.y(i, j) = (uf.top(i, j) + uf.bottom(i, j)) / 2.0;
  });
}

template <typename Float, Layout LAYOUT>
void dirichlet_bconds(Vector<Float, LAYOUT> u) {
  IGOR_ASSERT(u.x.nghost() == 1, "Expected exactly one ghost cell but got {}", u.x.nghost());
  for (Index i = 0; i < u.x.nx(); ++i) {
    u.x(i, -1)       = 0.0;
    u.x(i, u.x.ny()) = 0.0;
  }
  for (Index j = 0; j < u.x.ny(); ++j) {
    u.x(-1, j)       = 0.0;
    u.x(u.x.nx(), j) = 0.0;
  }

  for (Index i = 0; i < u.y.nx(); ++i) {
    u.y(i, -1)       = 0.0;
    u.y(i, u.y.ny()) = 0.0;
  }
  for (Index j = 0; j < u.y.ny(); ++j) {
    u.y(-1, j)       = 0.0;
    u.y(u.y.nx(), j) = 0.0;
  }
}

template <typename Float, Layout LAYOUT>
void dirichlet_bconds(FaceVector<Float, LAYOUT> uf) {
  IGOR_ASSERT(uf.x.nghost() == 1, "Expected exactly one ghost cell but got {}", uf.x.nghost());
  for (Index i = 0; i < uf.x.nx(); ++i) {
    uf.x(i, -1)            = 0.0;
    uf.x(i, 0)             = 0.0;
    uf.x(i, uf.x.ny() - 1) = 0.0;
    uf.x(i, uf.x.ny())     = 0.0;
  }
  for (Index j = 0; j < uf.x.ny(); ++j) {
    uf.x(-1, j)            = 0.0;
    uf.x(0, j)             = 0.0;
    uf.x(uf.x.nx() - 1, j) = 0.0;
    uf.x(uf.x.nx(), j)     = 0.0;
  }

  for (Index i = 0; i < uf.y.nx(); ++i) {
    uf.y(i, -1)            = 0.0;
    uf.y(i, 0)             = 0.0;
    uf.y(i, uf.y.ny() - 1) = 0.0;
    uf.y(i, uf.y.ny())     = 0.0;
  }
  for (Index j = 0; j < uf.y.ny(); ++j) {
    uf.y(-1, j)            = 0.0;
    uf.y(0, j)             = 0.0;
    uf.y(uf.y.nx() - 1, j) = 0.0;
    uf.y(uf.y.nx(), j)     = 0.0;
  }
}

template <typename Float, Layout LAYOUT>
constexpr auto L1_norm(const Grid<Float, LAYOUT>& grid, Scalar<Float, LAYOUT> s) -> Float {
  Float res = 0.0;
  grid.foreach_i([=, &res](Index i, Index j) { res += std::abs(s(i, j)) * grid.dv(); });
  return res;
}

auto main(int argc, char** argv) -> int {
  std::mt19937 rng(std::random_device{}());
  std::normal_distribution<Float> dist(0.0, 1e-3);

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

  auto uf = grid.alloc_face_vector();
  IGOR_DEFER(grid.free(uf););

  auto u = grid.alloc_vector();
  IGOR_DEFER(grid.free(u););

  auto div = grid.alloc_scalar();
  IGOR_DEFER(grid.free(div););

  auto rhs = grid.alloc_scalar();
  IGOR_DEFER(grid.free(rhs););

  VTKWriter writer(output_dir, grid);
  writer.add_field("p", p);
  writer.add_field("u", u);
  writer.add_field("div", div);
  writer.add_field("rhs", rhs);

  grid.foreach_face_i<Dimension::X>([=, &rng, &dist](Index i, Index j) { uf.x(i, j) = dist(rng); });
  grid.foreach_face_i<Dimension::Y>([=, &rng, &dist](Index i, Index j) { uf.y(i, j) = dist(rng); });
  dirichlet_bconds(uf);
  interpolate(grid, uf, u);

  const Float dt = 1e-1;
  calc_div(grid, uf, div);
  Igor::Info("||div|| = {}", L1_norm(grid, div));
  grid.foreach_i(FOREACH_FUNC { rhs(i, j) = div(i, j) / dt; });

  writer.write(0.0);

  const std::array<int, 2> ns   = {N + 2 * grid.nghost(), N + 2 * grid.nghost()};
  const std::array<Float, 2> Ls = {L + 2.0 * grid.dx(), L + 2.0 * grid.dy()};
  const std::array<int, 4> BCs  = {
      PoisFFT::NEUMANN, PoisFFT::NEUMANN, PoisFFT::NEUMANN, PoisFFT::NEUMANN};
  PoisFFT::Solver<2, Float> solver(ns.data(), Ls.data(), BCs.data(), PoisFFT::SPECTRAL);

  // const std::array<int, 4> ngs = {grid.nghost(), grid.nghost(), grid.nghost(), grid.nghost()};
  // solver.execute(p.data(), rhs.data(), ngs.data(), ngs.data());
  solver.execute(p.data(), rhs.data());

#if 0
  grid.foreach_face_i<Dimension::X>(FOREACH_FUNC {
    const auto dpdx  = (p(i, j) - p(i - 1, j)) / grid.dx();
    uf.x(i, j)      -= dt * dpdx;
  });
  grid.foreach_face_i<Dimension::Y>(FOREACH_FUNC {
    const auto dpdy  = (p(i, j) - p(i, j - 1)) / grid.dy();
    uf.y(i, j)      -= dt * dpdy;
  });
  dirichlet_bconds(uf);
  interpolate(grid, uf, u);
#else
  grid.foreach_i([=, &rng, &dist](Index i, Index j) {
    u.x(i, j) = dist(rng);
    u.y(i, j) = dist(rng);
  });
  dirichlet_bconds(u);
  interpolate(grid, u, uf);
  dirichlet_bconds(uf);
#endif

  calc_div(grid, uf, div);
  Igor::Info("||div|| = {}", L1_norm(grid, div));

  writer.write(1.0);
}
