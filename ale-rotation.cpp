#include <charconv>
#include <numbers>

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
  Float x;
  Float y;
};

constexpr Float pi       = std::numbers::pi_v<Float>;

constexpr Float x_min    = -0.5;
constexpr Float x_max    = 0.5;
constexpr Float y_min    = -0.5;
constexpr Float y_max    = 0.5;

constexpr Float CFL      = 0.5;
constexpr Float tend     = 1.0;
constexpr Float dt_write = tend / 100.0;
constexpr Float dt_max   = 1e-2;

// =================================================================================================
// = ALE Advection =================================================================================
// =================================================================================================
template <typename Float, Layout LAYOUT>
constexpr void advection_calc_flux(const Grid<Float, LAYOUT>& grid,
                                   const FaceVector<Float, LAYOUT> u,
                                   const Scalar<Float, LAYOUT> s,
                                   const FaceVector<Float, LAYOUT> w,
                                   const VertexScalar<Float, LAYOUT> x,
                                   const VertexScalar<Float, LAYOUT> y,
                                   FaceVector<Float, LAYOUT> F) {
  static auto sL = grid.alloc_face_vector();
  static auto sR = grid.alloc_face_vector();
  weno_reconstruction(grid, s, sL, sR);

  grid.template foreach_face_i<Dimension::X>(FOREACH_FUNC {
    const auto ux     = u.x(i, j);
    const auto uy     = (u.y(i - 1, j) + u.y(i - 1, j + 1) + u.y(i, j) + u.y(i, j + 1)) / 4.0;

    const auto wx     = w.x(i, j);
    const auto wy     = (w.y(i - 1, j) + w.y(i - 1, j + 1) + w.y(i, j) + w.y(i, j + 1)) / 4.0;

    const auto dxdeta = (x(i, j + 1) - x(i, j)) / grid.dy();
    const auto dydeta = (y(i, j + 1) - y(i, j)) / grid.dy();

    const auto si     = dydeta * (ux - wx) - dxdeta * (uy - wy) >= 0.0 ? sL.x(i, j) : sR.x(i, j);

    F.x(i, j)         = dydeta * (-si * ux + si * wx) - dxdeta * (-si * uy + si * wy);
  });

  grid.template foreach_face_i<Dimension::Y>(FOREACH_FUNC {
    const auto ux    = (u.x(i, j) + u.x(i, j - 1) + u.x(i + 1, j) + u.x(i + 1, j - 1)) / 4.0;
    const auto uy    = u.y(i, j);

    const auto wx    = (w.x(i, j) + w.x(i, j - 1) + w.x(i + 1, j) + w.x(i + 1, j - 1)) / 4.0;
    const auto wy    = w.y(i, j);

    const auto dxdxi = (x(i + 1, j) - x(i, j)) / grid.dx();
    const auto dydxi = (y(i + 1, j) - y(i, j)) / grid.dx();

    const auto si    = -dydxi * (ux - wx) + dxdxi * (uy - wy) >= 0.0 ? sL.y(i, j) : sR.y(i, j);

    F.y(i, j)        = -dydxi * (-si * ux + si * wx) + dxdxi * (-si * uy + si * wy);
  });
}

template <typename Float, Layout LAYOUT>
constexpr void advection_update_s(const Grid<Float, LAYOUT>& grid,
                                  Float dt,
                                  const FaceVector<Float, LAYOUT> F,
                                  const Scalar<Float, LAYOUT> s_old,
                                  const Scalar<Float, LAYOUT> J_old,
                                  const Scalar<Float, LAYOUT> J,
                                  Scalar<Float, LAYOUT> s) {
  grid.foreach_i(FOREACH_FUNC {
    s(i, j) = (J_old(i, j) * s_old(i, j) + dt * ((F.right(i, j) - F.left(i, j)) / grid.dx() +
                                                 (F.top(i, j) - F.bottom(i, j)) / grid.dy())) /
              J(i, j);
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
                               const FaceVector<Float, LAYOUT> w,
                               const VertexScalar<Float, LAYOUT> x,
                               const VertexScalar<Float, LAYOUT> y,
                               FaceVector<Float, LAYOUT> F) {
  grid.template foreach_face_i<Dimension::X>(FOREACH_FUNC {
    const auto wx     = w.x(i, j);
    const auto wy     = (w.y(i - 1, j) + w.y(i - 1, j + 1) + w.y(i, j) + w.y(i, j + 1)) / 4.0;
    const auto dxdeta = (x(i, j + 1) - x(i, j)) / grid.dy();
    const auto dydeta = (y(i, j + 1) - y(i, j)) / grid.dy();
    F.x(i, j)         = dydeta * wx - dxdeta * wy;
  });

  grid.template foreach_face_i<Dimension::Y>(FOREACH_FUNC {
    const auto wx    = (w.x(i, j - 1) + w.x(i, j) + w.x(i + 1, j - 1) + w.x(i + 1, j)) / 4.0;
    const auto wy    = w.y(i, j);
    const auto dxdxi = (x(i + 1, j) - x(i, j)) / grid.dx();
    const auto dydxi = (y(i + 1, j) - y(i, j)) / grid.dx();
    F.y(i, j)        = -dydxi * wx + dxdxi * wy;
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

constexpr auto cartesian2polar(Vec2 p) -> Vec2 {
  const auto r     = std::sqrt(Igor::sqr(p.x) + Igor::sqr(p.y));
  const auto theta = std::atan2(p.y, p.x);
  return {.x = theta, .y = r};
}

template <typename Float, Layout LAYOUT>
constexpr void calc_mesh_velocity(const Grid<Float, LAYOUT>& grid,
                                  const VertexScalar<Float, LAYOUT> x,
                                  const VertexScalar<Float, LAYOUT> y,
                                  FaceVector<Float, LAYOUT> w) {
  grid.template foreach_face_a<Dimension::X>(FOREACH_FUNC {
    const auto xi         = (x(i, j) + x(i, j + 1)) / 2.0;
    const auto yi         = (y(i, j) + y(i, j + 1)) / 2.0;
    const auto [theta, r] = cartesian2polar({.x = xi, .y = yi});

    const auto Utheta     = 2.0 * pi * r;
    w.x(i, j)             = -Utheta * std::sin(theta);
  });

  grid.template foreach_face_a<Dimension::Y>(FOREACH_FUNC {
    const auto xi         = (x(i, j) + x(i + 1, j)) / 2.0;
    const auto yi         = (y(i, j) + y(i + 1, j)) / 2.0;
    const auto [theta, r] = cartesian2polar({.x = xi, .y = yi});

    const auto Utheta     = 2.0 * pi * r;
    w.y(i, j)             = Utheta * std::cos(theta);
  });
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
  auto ui        = grid.alloc_vector();

  auto geo_x_old = grid.alloc_vertex_scalar();
  auto geo_x     = grid.alloc_vertex_scalar();
  auto geo_y_old = grid.alloc_vertex_scalar();
  auto geo_y     = grid.alloc_vertex_scalar();

  auto w         = grid.alloc_face_vector();
  auto wi        = grid.alloc_vector();

  auto J_old     = grid.alloc_scalar();
  auto J         = grid.alloc_scalar();
  auto FJ        = grid.alloc_face_vector();
  auto J_geom    = grid.alloc_scalar();

  auto xi        = grid.alloc_scalar();
  auto yi        = grid.alloc_scalar();

  Float dt       = 0.0;
  Float t        = 0.0;
  // -----------------------------------------------------------------------------------------------

  const BConds<Float> bconds{
      .left   = Neumann{},
      .right  = Neumann{},
      .bottom = Neumann{},
      .top    = Neumann{},
  };

  grid.foreach_vertex_a(FOREACH_FUNC {
    geo_x(i, j) = grid.x(i);
    geo_y(i, j) = grid.y(j);
  });
  interpolate_vertex(grid, geo_x, xi);
  interpolate_vertex(grid, geo_y, yi);

  ale_calc_J_geom(grid, geo_x, geo_y, J);
  // Make sure that these values are never used
  fill_ghost(grid, J, std::numeric_limits<Float>::quiet_NaN());
  ale_calc_J_geom(grid, geo_x, geo_y, J_geom);

#if STATIONARY
  fill(u, 0.0);
#else
  calc_mesh_velocity(grid, geo_x, geo_y, u);
  grid.foreach_face_i<Dimension::X>(FOREACH_FUNC { u.x(i, j) *= -1.0; });
  grid.foreach_face_i<Dimension::Y>(FOREACH_FUNC { u.y(i, j) *= -1.0; });
#endif
  apply_velocity_bconds(grid, bconds, bconds, u, t);
  interpolate(grid, u, ui);

  grid.foreach_i(FOREACH_FUNC {
    const auto x = grid.xm(i);
    const auto y = grid.ym(j);
    s(i, j)      = static_cast<Float>(Igor::sqr(x + 0.25) + Igor::sqr(y) < Igor::sqr(0.1));
  });
  apply_bconds(grid, bconds, s, 0.0);

  calc_mesh_velocity(grid, geo_x, geo_y, w);
  interpolate(grid, w, wi);

  const auto moving_dir     = output_dir + "moving";
  const auto stationary_dir = output_dir + "stationary";
  {
    std::error_code ec;
    std::filesystem::create_directories(moving_dir, ec);
    if (ec) {
      Igor::Warn("Could not create directory `{}`: {}", moving_dir, ec.message());
      return 1;
    }

    std::filesystem::create_directories(stationary_dir, ec);
    if (ec) {
      Igor::Warn("Could not create directory `{}`: {}", stationary_dir, ec.message());
      return 1;
    }
  }

  VTKWriter writer_moving_grid(moving_dir, grid);
  writer_moving_grid.add_field("s", s);
  writer_moving_grid.add_field("u", ui);
  writer_moving_grid.add_field("w", wi);
  writer_moving_grid.add_field("x", xi);
  writer_moving_grid.add_field("y", yi);
  writer_moving_grid.add_field("J", J);
  writer_moving_grid.add_field("J_geom", J_geom);
  if (!writer_moving_grid.write(geo_x, geo_y, t)) { return 1; }

  VTKWriter writer_stationary_grid(stationary_dir, grid);
  writer_stationary_grid.add_field("s", s);
  writer_stationary_grid.add_field("u", ui);
  writer_stationary_grid.add_field("w", wi);
  writer_stationary_grid.add_field("x", xi);
  writer_stationary_grid.add_field("y", yi);
  writer_stationary_grid.add_field("J", J);
  writer_stationary_grid.add_field("J_geom", J_geom);
  if (!writer_stationary_grid.write(t)) { return 1; }

  Stats s_stats      = stats(grid, s, J);
  const auto s0_sum  = s_stats.sum;
  Stats J_stats      = stats(grid, J, J);
  Stats J_geom_stats = stats(grid, J_geom, J);

  Monitor<Float> monitor(output_dir + "/monitor.log");
  monitor.add_variable(&t, "t");
  monitor.add_variable(&dt, "dt");
  monitor.add_variable(&s_stats.min, "min(s)");
  monitor.add_variable(&s_stats.max, "max(s)");
  monitor.add_variable(&s_stats.sum, "sum(s)");
  monitor.add_variable(&J_stats.min, "min(J)");
  monitor.add_variable(&J_stats.max, "max(J)");
  monitor.add_variable(&J_geom_stats.min, "min(J_geom)");
  monitor.add_variable(&J_geom_stats.max, "max(J_geom)");
  monitor.write();

  IGOR_TIME_SCOPE("Solver")
  while (t < tend) {

    dt = std::min({
        adjust_dt(grid, u, 1.0, 0.0, CFL),
        adjust_dt(grid, w, 1.0, 0.0, CFL),
        dt_max,
        dt_write,
        tend - t,
    });

    copy(geo_x, geo_x_old);
    copy(geo_y, geo_y_old);
    copy(J, J_old);
    copy(s, s_old);

    for (Index sub_iter = 0; sub_iter < 2; ++sub_iter) {
      const auto local_dt = sub_iter == 0 ? dt / 2.0 : dt;

      ale_calc_J_flux(grid, w, geo_x, geo_y, FJ);
      ale_update_J(grid, local_dt, FJ, J_old, J);
      // Make sure that these values are never used
      fill_ghost(grid, J, std::numeric_limits<Float>::quiet_NaN());

      advection_calc_flux(grid, u, s, w, geo_x, geo_y, Fs);
      advection_update_s(grid, local_dt, Fs, s_old, J_old, J, s);
      apply_bconds(grid, bconds, s, 0.0);

      // #define SPIRAL

      // This is a hack because I don't know how to properly handle the ghost cells otherwise
      grid.foreach_vertex_a(FOREACH_FUNC {
#ifndef SPIRAL
        const auto [theta, r] = cartesian2polar({.x = geo_x(i, j), .y = geo_y(i, j)});
#else
        const auto [theta, r] = cartesian2polar({.x = grid.x(i), .y = grid.y(j)});
#endif  // SPIRAL
        const auto Utheta = 2.0 * pi * r;
        const auto wx     = -Utheta * std::sin(theta);
        const auto wy     = Utheta * std::cos(theta);
        geo_x(i, j)       = geo_x_old(i, j) + local_dt * wx;
        geo_y(i, j)       = geo_y_old(i, j) + local_dt * wy;
      });
#ifndef SPIRAL
      calc_mesh_velocity(grid, geo_x, geo_y, w);
#endif

#ifndef STATIONARY
      calc_mesh_velocity(grid, geo_x, geo_y, u);
      grid.foreach_face_i<Dimension::X>(FOREACH_FUNC { u.x(i, j) *= -1.0; });
      grid.foreach_face_i<Dimension::Y>(FOREACH_FUNC { u.y(i, j) *= -1.0; });
      apply_velocity_bconds(grid, bconds, bconds, u);
#endif
    }

    ale_calc_J_geom(grid, geo_x, geo_y, J_geom);

    s_stats      = stats(grid, s, J);
    J_stats      = stats(grid, J, J);
    J_geom_stats = stats(grid, J_geom, J);
    interpolate(grid, w, wi);
    interpolate_vertex(grid, geo_x, xi);
    interpolate_vertex(grid, geo_y, yi);
    t += dt;
    if (should_save(t, dt, dt_write, tend)) {
      if (!writer_moving_grid.write(geo_x, geo_y, t)) { return 1; }
      if (!writer_stationary_grid.write(t)) { return 1; }
    }
    monitor.write();
  }

  Igor::Info("abserr(s) = {}", std::abs(s_stats.sum - s0_sum));
  Igor::Info("relerr(s) = {}", std::abs((s_stats.sum - s0_sum) / s0_sum));

#ifndef SPIRAL
  Igor::Info("abserr(J) = {}", std::abs(J_stats.max - 1.0));
#endif
  Igor::Info("abserr(J, J_geom) = {}", std::abs(J_stats.max - J_geom_stats.max));

  Igor::Info("Ok.");
}
