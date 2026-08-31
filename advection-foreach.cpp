#include <charconv>

#include <Igor/Defer.hpp>
#include <Igor/Logging.hpp>
#include <Igor/Timer.hpp>

#include "Grid.hpp"
#include "IO.hpp"
#include "Quadrature.hpp"
#include "VTKWriter.hpp"

using Float                  = double;

constexpr Float x_min        = 0.0;
constexpr Float x_max        = 1.0;
constexpr Float y_min        = 0.0;
constexpr Float y_max        = 1.0;

constexpr Float tend         = 0.5;
constexpr Float dt_write     = tend / 100.0;

constexpr Index NUM_SUB_ITER = 5;

constexpr Float a            = 0.5;

constexpr auto u_analytical(Float x, Float /*y*/, Float t) -> Float {
  return x < 0.5 + a * t ? 1.0 : 0.0;
}

constexpr auto left_biased_reconstruction(Float um2, Float um1, Float u0, Float up1, Float up2)
    -> Float {
  constexpr Float eps  = 1e-6;
  constexpr Float dL_0 = 1.0 / 10.0;
  constexpr Float dL_1 = 6.0 / 10.0;
  constexpr Float dL_2 = 3.0 / 10.0;

  const auto q0        = (2.0 * um2 - 7.0 * um1 + 11.0 * u0) / 6.0;
  const auto q1        = (-1.0 * um1 + 5.0 * u0 + 2.0 * up1) / 6.0;
  const auto q2        = (2.0 * u0 + 5.0 * up1 - 1.0 * up2) / 6.0;

  const auto b0        = (13.0 / 12.0) * Igor::sqr(um2 - 2.0 * um1 + u0) +
                         (1.0 / 4.0) * Igor::sqr(um2 - 4.0 * um1 + 3.0 * u0);
  const auto b1        = (13.0 / 12.0) * Igor::sqr(um1 - 2.0 * u0 + up1) +  //
                         (1.0 / 4.0) * Igor::sqr(um1 - up1);
  const auto b2        = (13.0 / 12.0) * Igor::sqr(u0 - 2.0 * up1 + up2) +
                         (1.0 / 4.0) * Igor::sqr(3.0 * u0 - 4.0 * up1 + up2);

  const auto a0        = dL_0 / Igor::sqr(b0 + eps);
  const auto a1        = dL_1 / Igor::sqr(b1 + eps);
  const auto a2        = dL_2 / Igor::sqr(b2 + eps);
  const auto asum      = a0 + a1 + a2;

  return (a0 * q0 + a1 * q1 + a2 * q2) / asum;
}

constexpr auto right_biased_reconstruction(Float um1, Float u0, Float up1, Float up2, Float up3)
    -> Float {
  constexpr Float eps  = 1e-6;
  constexpr Float dR_0 = 3.0 / 10.0;
  constexpr Float dR_1 = 6.0 / 10.0;
  constexpr Float dR_2 = 1.0 / 10.0;

  const Float q0       = (-1.0 * um1 + 5.0 * u0 + 2.0 * up1) / 6.0;
  const Float q1       = (2.0 * u0 + 5.0 * up1 - 1.0 * up2) / 6.0;
  const Float q2       = (11.0 * up1 - 7.0 * up2 + 2.0 * up3) / 6.0;

  const auto b0        = (13.0 / 12.0) * Igor::sqr(um1 - 2.0 * u0 + up1) +
                         (1.0 / 4.0) * Igor::sqr(um1 - 4.0 * u0 + 3.0 * up1);
  const auto b1        = (13.0 / 12.0) * Igor::sqr(u0 - 2.0 * up1 + up2) +  //
                         (1.0 / 4.0) * Igor::sqr(u0 - up2);
  const auto b2        = (13.0 / 12.0) * Igor::sqr(up1 - 2.0 * up2 + up3) +
                         (1.0 / 4.0) * Igor::sqr(3.0 * up1 - 4.0 * up2 + up3);

  const auto a0        = dR_0 / Igor::sqr(b0 + eps);
  const auto a1        = dR_1 / Igor::sqr(b1 + eps);
  const auto a2        = dR_2 / Igor::sqr(b2 + eps);
  const auto asum      = a0 + a1 + a2;

  return (a0 * q0 + a1 * q1 + a2 * q2) / asum;
}

template <typename Float, Layout LAYOUT>
constexpr void calc_reconstruction(const Grid<Float, LAYOUT>& grid,
                                   const Scalar<Float, LAYOUT> u,
                                   FaceVector<Float, LAYOUT> uL,
                                   FaceVector<Float, LAYOUT> uR) {
  grid.template foreach_face_i<Dimension::X>(FOREACH_FUNC {
    uL.x(i, j) =
        left_biased_reconstruction(u(i - 3, j), u(i - 2, j), u(i - 1, j), u(i, j), u(i + 1, j));
    uR.x(i, j) =
        right_biased_reconstruction(u(i - 2, j), u(i - 1, j), u(i, j), u(i + 1, j), u(i + 2, j));
  });

  grid.template foreach_face_i<Dimension::Y>(FOREACH_FUNC {
    uL.y(i, j) =
        left_biased_reconstruction(u(i, j - 3), u(i, j - 2), u(i, j - 1), u(i, j), u(i, j + 1));
    uR.y(i, j) =
        right_biased_reconstruction(u(i, j - 2), u(i, j - 1), u(i, j), u(i, j + 1), u(i, j + 2));
  });
}

template <typename Float, Layout LAYOUT>
constexpr void calc_flux(const Grid<Float, LAYOUT>& grid,
                         const FaceVector<Float, LAYOUT> uL,
                         const FaceVector<Float, LAYOUT> uR,
                         FaceVector<Float, LAYOUT> F) {
  grid.template foreach_face_i<Dimension::X>(
      FOREACH_FUNC { F.x(i, j) = a * (a >= 0.0 ? uL.x(i, j) : uR.x(i, j)); });
  grid.template foreach_face_i<Dimension::Y>(
      FOREACH_FUNC { F.y(i, j) = a * (a >= 0.0 ? uL.y(i, j) : uR.y(i, j)); });
}

template <typename Float, Layout LAYOUT>
constexpr void update_u(const Grid<Float, LAYOUT>& grid,
                        Float dt,
                        const FaceVector<Float, LAYOUT> F,
                        const Scalar<Float, LAYOUT> u_old,
                        Scalar<Float, LAYOUT> u) {
  grid.foreach_i(FOREACH_FUNC {
    u(i, j) = u_old(i, j) - dt * ((F.right(i, j) - F.left(i, j)) / grid.dx() +
                                  (F.top(i, j) - F.bottom(i, j)) / grid.dy());
  });
}

template <typename Float, Layout LAYOUT>
constexpr void boundary_conditions(Scalar<Float, LAYOUT> u) {
  // Homogeneous Neumann
  for (Index i = 0; i < u.nx(); ++i) {
    for (Index j = -u.nghost(); j < 0; ++j) {
      u(i, j) = u(i, -1 - j);
    }
    for (Index j = u.ny(); j < u.ny() + u.nghost(); ++j) {
      u(i, j) = u(i, 2 * u.ny() - 1 - j);
    }
  }

  for (Index j = 0; j < u.ny(); ++j) {
    for (Index i = -u.nghost(); i < 0; ++i) {
      u(i, j) = u(-1 - i, j);
    }
    for (Index i = u.nx(); i < u.nx() + u.nghost(); ++i) {
      u(i, j) = u(2 * u.ny() - 1 - i, j);
    }
  }
}

auto run(const std::string& output_base_dir, Index N) -> bool {
  const auto output_dir = output_base_dir + "/" + std::to_string(N);
  if (!init_output_directory(output_dir)) { return false; }

  Grid<Float, Layout::C> grid(x_min, x_max, N, y_min, y_max, N, 3);

  auto u_old = grid.alloc_scalar();
  auto u     = grid.alloc_scalar();
  auto uL    = grid.alloc_face_vector();
  auto uR    = grid.alloc_face_vector();
  auto F     = grid.alloc_face_vector();

  Float dt   = 0.0;
  Float t    = 0.0;

  grid.foreach_i(FOREACH_FUNC {
    const Float x = grid.x_min() + (i + 0.5) * grid.dx();
    const Float y = grid.y_min() + (j + 0.5) * grid.dy();
    u(i, j)       = u_analytical(x, y, t);
  });
  boundary_conditions(u);

  VTKWriter writer(output_dir, grid);
  writer.add_field("u", u);
  if (!writer.write(t)) { return false; }

  while (t < tend) {
    dt = 0.5 * std::min(grid.dx(), grid.dy()) / a;
    dt = std::min(dt, tend - t);

    copy(u, u_old);

#if 1
    for (Index sub_iter = 0; sub_iter < NUM_SUB_ITER; ++sub_iter) {
      grid.foreach_a(FOREACH_FUNC { u(i, j) = 0.5 * (u(i, j) + u_old(i, j)); });

      calc_reconstruction(grid, u, uL, uR);
      calc_flux(grid, uL, uR, F);
      update_u(grid, dt, F, u_old, u);
      boundary_conditions(u);
    }
#else
    calc_reconstruction(u, uL, uR);
    calc_flux(uL, uR, F);
    update_u(grid.dx(), grid.dy(), dt / 2.0, F, u_old, u);
    boundary_conditions(u);

    calc_reconstruction(u, uL, uR);
    calc_flux(uL, uR, F);
    update_u(grid.dx(), grid.dy(), dt, F, u_old, u);
    boundary_conditions(u);
#endif

    t += dt;
    if (should_save(t, dt, dt_write, tend) && !writer.write(t)) { return false; }
  }

  Float L1_error = 0.0;
  grid.foreach_i<Exec::SERIAL>([=, &L1_error](Index i, Index j) {
    const Float x0 = grid.x_min() + i * grid.dx();
    const Float x1 = x0 + grid.dx();
    const Float y0 = grid.y_min() + j * grid.dy();
    const Float y1 = y0 + grid.dy();

    const Float u_exp =
        quadrature([&](Float x, Float y) { return u_analytical(x, y, t); }, x0, x1, y0, y1) /
        grid.dv(i, j);
    L1_error += grid.dv(i, j) * std::abs(u(i, j) - u_exp);
  });
  Igor::Info("{} => {}", N, L1_error);

  return true;
}

auto main(int argc, char** argv) -> int {
  const auto usage_str = Igor::detail::format("Usage: {} <max. level>", argv[0]);
  if (argc < 2) {
    Igor::Error("{}", usage_str);
    return 1;
  }

  Index max_level = 0;
  if (std::from_chars(argv[1], argv[1] + std::strlen(argv[1]), max_level).ec != std::errc{} ||
      max_level <= 0) {
    Igor::Error("{}", usage_str);
    Igor::Error("  Invalid max. level `{}`", argv[1]);
    return 1;
  }

  const auto output_dir = get_output_directory();

  IGOR_TIME_SCOPE("Scaling") {
    // for (Index level = 2; level <= max_level; ++level) {
    //   const auto N = 1 << level;
    //   run(output_dir, N);
    // }
    run(output_dir, 1 << max_level);
  }

  Igor::Info("Ok.");
}
