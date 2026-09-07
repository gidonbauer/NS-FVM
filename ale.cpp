#include <charconv>

#include <poisfft.h>

#include <Igor/Defer.hpp>
#include <Igor/Logging.hpp>
#include <Igor/Math.hpp>
#include <Igor/Timer.hpp>

#include "BoundaryConditions.hpp"
#include "Common.hpp"
#include "Grid.hpp"
#include "IO.hpp"
#include "Mac.hpp"
#include "Monitor.hpp"
#include "VTKWriter.hpp"
#include "WENO5.hpp"

using Float = double;
struct Vec2 {
  Float x, y;
};

constexpr Float x_min    = 0.0;
constexpr Float x_max    = 1.0;
constexpr Float y_min    = 0.0;
constexpr Float y_max    = 1.0;

constexpr Vec2 U         = {.x = 1.0, .y = 1.0};
constexpr Vec2 W         = {.x = 0.5, .y = 0.5};

constexpr Float CFL      = 0.5;
constexpr Float tend     = 2.0;
constexpr Float dt_write = tend / 100.0;

// =================================================================================================
// = ALE Advection =================================================================================
// =================================================================================================
template <typename Float, Layout LAYOUT>
constexpr void advection_calc_flux(const Grid<Float, LAYOUT>& grid,
                                   const FaceVector<Float, LAYOUT> u,
                                   const Scalar<Float, LAYOUT> s,
                                   const VertexScalar<Float, LAYOUT> x,
                                   const VertexScalar<Float, LAYOUT> y,
                                   FaceVector<Float, LAYOUT> F) {
  static auto sL = grid.alloc_face_vector();
  static auto sR = grid.alloc_face_vector();
  weno_reconstruction(grid, s, sL, sR);

  grid.template foreach_face_i<Dimension::X>(FOREACH_FUNC {
    const auto si     = u.x(i, j) >= 0.0 ? sL.x(i, j) : sR.x(i, j);
    const auto ui     = u.x(i, j);
    const auto vi     = (u.y(i - 1, j) + u.y(i - 1, j + 1) + u.y(i, j) + u.y(i, j + 1)) / 4.0;
    const auto dxdeta = (x(i, j + 1) - x(i, j)) / grid.dy();
    const auto dydeta = (y(i, j + 1) - y(i, j)) / grid.dy();
    F.x(i, j)         = dydeta * (-si * ui + si * W.x) - dxdeta * (-si * vi + si * W.y);
  });

  grid.template foreach_face_i<Dimension::Y>(FOREACH_FUNC {
    const auto si    = u.y(i, j) >= 0.0 ? sL.y(i, j) : sR.y(i, j);
    const auto ui    = (u.x(i - 1, j) + u.x(i - 1, j + 1) + u.x(i, j) + u.x(i, j + 1)) / 4.0;
    const auto vi    = u.y(i, j);
    const auto dxdxi = (x(i + 1, j) - x(i, j)) / grid.dx();
    const auto dydxi = (y(i + 1, j) - y(i, j)) / grid.dx();
    F.y(i, j)        = -dydxi * (-si * ui + si * W.x) + dxdxi * (-si * vi + si * W.y);
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
constexpr void ale_calc_J_geom(const Grid<Float, LAYOUT>& grid,
                               const VertexScalar<Float, LAYOUT> x,
                               const VertexScalar<Float, LAYOUT> y,
                               Scalar<Float, LAYOUT> J) {
  grid.foreach_i(FOREACH_FUNC {
    const auto dxdxi =
        ((x(i + 1, j) + x(i + 1, j + 1)) - (x(i, j) + x(i, j + 1))) / (2.0 * grid.dx());
    const auto dxdeta =
        ((x(i, j + 1) + x(i + 1, j + 1)) - (x(i, j) + x(i + 1, j))) / (2.0 * grid.dy());

    const auto dydxi =
        ((y(i + 1, j) + y(i + 1, j + 1)) - (y(i, j) + y(i, j + 1))) / (2.0 * grid.dx());
    const auto dydeta =
        ((y(i, j + 1) + y(i + 1, j + 1)) - (y(i, j) + y(i + 1, j))) / (2.0 * grid.dy());

    J(i, j) = dxdxi * dydeta - dxdeta * dydxi;
  });
}

template <typename Float, Layout LAYOUT>
constexpr void ale_calc_J_flux(const Grid<Float, LAYOUT>& grid,
                               const VertexScalar<Float, LAYOUT> x,
                               const VertexScalar<Float, LAYOUT> y,
                               FaceVector<Float, LAYOUT> F) {
  grid.template foreach_face_i<Dimension::X>(FOREACH_FUNC {
    const auto dxdeta = (x(i, j + 1) - x(i, j)) / grid.dy();
    const auto dydeta = (y(i, j + 1) - y(i, j)) / grid.dy();
    F.x(i, j)         = dydeta * W.x - dxdeta * W.y;
  });

  grid.template foreach_face_i<Dimension::Y>(FOREACH_FUNC {
    const auto dxdxi = (x(i + 1, j) - x(i, j)) / grid.dx();
    const auto dydxi = (y(i + 1, j) - y(i, j)) / grid.dx();
    F.y(i, j)        = -dydxi * W.x + dxdxi * W.y;
  });
}

template <typename Float, Layout LAYOUT>
constexpr void ale_update_J(const Grid<Float, LAYOUT>& grid,
                            Float dt,
                            const FaceVector<Float, LAYOUT> F,
                            const Scalar<Float, LAYOUT> J_old,
                            Scalar<Float, LAYOUT> J) {
  grid.foreach_i(FOREACH_FUNC {
    J(i, j) = J_old(i, j) + dt * ((F.right(i, j) - F.left(i, j)) / grid.dx() +
                                  (F.top(i, j) - F.bottom(i, j)) / grid.dy());
  });
}
// =================================================================================================
// = ALE Advection =================================================================================
// =================================================================================================

template <typename Float, Layout LAYOUT>
constexpr void interpolate_vertex(const Grid<Float, LAYOUT>& grid,
                                  const VertexScalar<Float, LAYOUT> v,
                                  Scalar<Float, LAYOUT> s) {
  grid.foreach_i(
      FOREACH_FUNC { s(i, j) = (v(i, j) + v(i + 1, j) + v(i, j + 1) + v(i + 1, j + 1)) / 4.0; });
}

// =================================================================================================
auto main(int argc, char** argv) -> int {
  const auto usage_str = Igor::detail::format("Usage: {} <grid size>", argv[0]);
  if (argc < 2) {
    Igor::Error("{}", usage_str);
    return 1;
  }

  Index N = 0;
  if (std::from_chars(argv[1], argv[1] + std::strlen(argv[1]), N).ec != std::errc{} || N <= 0) {
    Igor::Error("{}", usage_str);
    Igor::Error("  Invalid grid size `{}`", argv[1]);
    return 1;
  }

  const auto output_dir = get_output_directory();
  if (!init_output_directory(output_dir)) { return 1; }

  Grid<Float> grid(x_min, x_max, N, y_min, y_max, N, 3);

  // -----------------------------------------------------------------------------------------------
  auto s_old     = grid.alloc_scalar();
  auto s         = grid.alloc_scalar();
  auto Fs        = grid.alloc_face_vector();

  auto u         = grid.alloc_face_vector();

  auto geo_x_old = grid.alloc_vertex_scalar();
  auto geo_x     = grid.alloc_vertex_scalar();
  auto geo_y_old = grid.alloc_vertex_scalar();
  auto geo_y     = grid.alloc_vertex_scalar();

  auto J_old     = grid.alloc_scalar();
  auto J         = grid.alloc_scalar();
  auto FJ        = grid.alloc_face_vector();

  auto xi        = grid.alloc_scalar();
  auto yi        = grid.alloc_scalar();

  Float dt       = 0.0;
  Float t        = 0.0;
  // -----------------------------------------------------------------------------------------------

  const BConds<Float> bconds{
      .left   = Periodic{},
      .right  = Periodic{},
      .bottom = Periodic{},
      .top    = Periodic{},
  };

  fill(u.x, U.x);
  fill(u.y, U.y);
  apply_velocity_bconds(grid, bconds, bconds, u, t);

  grid.foreach_i(FOREACH_FUNC {
    const auto x = grid.xm(i);
    const auto y = grid.ym(j);
    s(i, j)      = static_cast<Float>(Igor::sqr(x - 0.5) + Igor::sqr(y - 0.5) < Igor::sqr(0.1));
  });
  apply_bconds(grid, bconds, s, 0.0);

  grid.foreach_vertex_a(FOREACH_FUNC {
    geo_x(i, j) = grid.x(i);
    geo_y(i, j) = grid.y(j);
  });
  interpolate_vertex(grid, geo_x, xi);
  interpolate_vertex(grid, geo_y, yi);

  ale_calc_J_geom(grid, geo_x, geo_y, J);

  VTKWriter writer(output_dir, grid);
  writer.add_field("s", s);
  writer.add_field("x", xi);
  writer.add_field("y", yi);
  writer.add_field("J", J);
  if (!writer.write(t)) { return 1; }

  Stats s_stats = stats(grid, s);
  Stats J_stats = stats(grid, J);

  Monitor<Float> monitor(output_dir + "/monitor.log");
  monitor.add_variable(&t, "t");
  monitor.add_variable(&dt, "dt");
  monitor.add_variable(&s_stats.min, "min(s)");
  monitor.add_variable(&s_stats.max, "max(s)");
  monitor.add_variable(&J_stats.min, "min(J)");
  monitor.add_variable(&J_stats.max, "max(J)");
  monitor.write();

  IGOR_TIME_SCOPE("Solver")
  while (t < tend) {

    dt = std::min({
        adjust_dt(grid, u, 1.0, 0.0, CFL),
        dt_write,
        tend - t,
    });

    copy(geo_x, geo_x_old);
    copy(geo_y, geo_y_old);
    copy(J, J_old);
    copy(s, s_old);

    for (Index sub_iter = 0; sub_iter < 2; ++sub_iter) {
      const auto local_dt = sub_iter == 0 ? dt / 2.0 : dt;

      advection_calc_flux(grid, u, s, geo_x, geo_y, Fs);
      advection_update_s(grid, local_dt, Fs, s_old, s);
      apply_bconds(grid, bconds, s, 0.0);

      ale_calc_J_flux(grid, geo_x, geo_y, FJ);
      ale_update_J(grid, local_dt, FJ, J_old, J);

      grid.foreach_vertex_a(FOREACH_FUNC {
        geo_x(i, j) = geo_x_old(i, j) + local_dt * W.x;
        geo_y(i, j) = geo_y_old(i, j) + local_dt * W.y;
      });
    }

    s_stats = stats(grid, s);
    J_stats = stats(grid, J);
    interpolate_vertex(grid, geo_x, xi);
    interpolate_vertex(grid, geo_y, yi);
    t += dt;
    if (should_save(t, dt, dt_write, tend)) {
      if (!writer.write(t)) { return 1; }
    }
    monitor.write();
  }

  Igor::Info("Ok.");
}
