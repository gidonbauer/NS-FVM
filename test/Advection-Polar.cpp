#include <charconv>
#include <numbers>

#include <Igor/Defer.hpp>
#include <Igor/Logging.hpp>
#include <Igor/Math.hpp>
#include <Igor/Timer.hpp>

#include "Advection-Diffusion.hpp"
#include "BoundaryConditions.hpp"
#include "Common.hpp"
#include "Grid.hpp"
#include "HDFWriter.hpp"
#include "IO.hpp"
#include "Mac.hpp"
#include "Monitor.hpp"

#include "Test-Common.hpp"

using Float               = double;
constexpr Float pi        = std::numbers::pi_v<Float>;

constexpr Float theta_min = 0.0;
constexpr Float theta_max = 2.0 * pi;
constexpr Float r_min     = 1.0;
constexpr Float r_max     = 2.0;

constexpr Float D         = 0.0;
constexpr Float U         = 2.0 * pi;

constexpr Float CFL       = 0.5;
constexpr Float tend      = 1.0;
constexpr Float dt_write  = tend / 100.0;

// =================================================================================================
constexpr auto s_analytical(Float theta, Float r) -> Float {
  // return static_cast<Float>(Igor::sqr(theta - pi) + Igor::sqr(r - 1.5) < Igor::sqr(0.25));
  return r * std::sin(theta);
}

// =================================================================================================
namespace Expected {
constexpr std::array ns  = {16, 32, 64, 128};
constexpr std::array L1s = {1.02105635e-01, 2.46850481e-02, 6.02777288e-03, 1.48930342e-03};
static_assert(ns.size() == L1s.size());
constexpr auto L1(Index n) { return interp_n2(ns, L1s, n); }
}  // namespace Expected
// =================================================================================================

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

  const std::string output_dir = "./test/output/Advection-Polar-" + std::to_string(N) + '/';
  if (!init_output_directory(output_dir)) { return 1; }

  Grid<Float> grid(theta_min, theta_max, N, r_min, r_max, N, 3, Coordinates::POLAR);

  auto u     = grid.alloc_face_vector();
  auto ui    = grid.alloc_vector();

  auto s_old = grid.alloc_scalar();
  auto s     = grid.alloc_scalar();
  auto FS    = grid.alloc_face_vector();

  Float dt   = 0.0;
  Float t    = 0.0;

  const BConds<Float> bconds{
      .left   = Periodic{},
      .right  = Periodic{},
      .bottom = Neumann{},
      .top    = Neumann{},
  };

  grid.foreach_i(FOREACH_FUNC { s(i, j) = s_analytical(grid.xm(i), grid.ym(j)); });
  apply_bconds(grid, bconds, s, t);

  grid.foreach_face_i<Dimension::X>(FOREACH_FUNC {
    const auto r = grid.ym(j);
    u.x(i, j)    = r * U;
  });
  apply_velocity_bconds(grid, bconds, bconds, u);
  interpolate(grid, u, ui);

  HDFWriter writer(output_dir, grid);
  writer.add_field("u", ui);
  writer.add_field("s", s);
  if (!writer.write(t)) { return 1; }

  Stats s_stats = stats(grid, s);

  Monitor<Float> monitor(output_dir + "/monitor.log");
  monitor.add_variable(&t, "t");
  monitor.add_variable(&dt, "dt");
  monitor.add_variable(&s_stats.min, "min(s)");
  monitor.add_variable(&s_stats.max, "max(s)");
  monitor.write();

  IGOR_TIME_SCOPE("Advection-Polar-" + std::to_string(N))
  while (t < tend) {

    dt = std::min({
        adjust_dt(grid, u, 1.0, 0.0, CFL),
        advection_adjust_dt(grid, D, CFL),
        tend - t,
    });

    copy(s, s_old);

    for (Index sub_iter = 0; sub_iter < 2; ++sub_iter) {
      const auto local_dt = sub_iter == 0 ? dt / 2.0 : dt;
      advection_calc_flux(grid, u, s, D, FS);
      advection_update_s(grid, local_dt, FS, s_old, s);
      apply_bconds(grid, bconds, s, t);
    }

    s_stats  = stats(grid, s);
    t       += dt;
    if (should_save(t, dt, dt_write, tend)) {
      if (!writer.write(t)) { return 1; }
    }
    monitor.write();
  }

  Float L1 = 0.0;
  grid.foreach_i<Exec::SERIAL>([=, &L1](Index i, Index j) {
    const auto s_exp  = s_analytical(grid.xm(i), grid.ym(j));
    L1               += std::abs(s_exp - s(i, j)) * grid.dv(i, j);
  });

  if (L1 > 1.1 * Expected::L1(N)) {
    Igor::Error("s error does not match expected value: expected {:.8e} but got {:.8e}",
                Expected::L1(N),
                L1);
    return 1;
  }
  return 0;
}
