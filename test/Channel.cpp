#include <charconv>

#include <Igor/Defer.hpp>
#include <Igor/Logging.hpp>
#include <Igor/Math.hpp>
#include <Igor/Timer.hpp>

#if !FFT_POISSON && !MG_POISSON || FFT_POISSON && MG_POISSON
#error "Exactly one of FFT_POISSON or MG_POISSON must be true"
#endif

#if FFT_POISSON
#include <poisfft.h>
#else
#include "MultigridPoisson.hpp"
#endif

#include "BoundaryConditions.hpp"
#include "Common.hpp"
#include "Grid.hpp"
#include "IO.hpp"
#include "Mac.hpp"
#include "Monitor.hpp"
#include "VTKWriter.hpp"

using Float              = double;

constexpr Index ratio    = 1;
constexpr Float y_min    = 0.0;
constexpr Float y_max    = 1.0;
constexpr Float x_min    = 0.0;
constexpr Float x_max    = ratio * y_max;

constexpr Float rho      = 1.0;
constexpr Float mu       = 1.0;
constexpr Float Uavg     = 1.0;

constexpr Float CFL      = 0.7;
constexpr Float tend     = 0.5;
constexpr Float dt_write = tend / 100.0;

// =================================================================================================
namespace Expected {

constexpr std::array ns  = {16, 32, 64, 128, 256};
constexpr std::array L1s = {
    2.11220054e-03, 5.32794280e-04, 1.33995883e-04, 3.35113462e-05, 8.37842234e-06};
static_assert(ns.size() == L1s.size());

template <std::size_t N>
constexpr auto interp_n2(const std::array<Index, N>& xs, const std::array<double, N>& es, Index n)
    -> double {
  static_assert(N >= 2, "need at least two samples");

  // Exact hits: return the tabulated value bit-for-bit.
  for (std::size_t i = 0; i < N; ++i) {
    if (xs[i] == n) { return es[i]; }
  }

  const auto coeff = [&](std::size_t i) {
    const auto xi = static_cast<double>(xs[i]);
    return es[i] * xi * xi;
  };
  const auto nd = static_cast<double>(n);

  // Extrapolation: hold C at the nearest end -> pure quadratic scaling.
  if (n < xs.front()) { return coeff(0) / (nd * nd); }
  if (n > xs.back()) { return coeff(N - 1) / (nd * nd); }

  // Bracket: xs[i] < n < xs[i + 1]
  std::size_t i = 0;
  while (xs[i + 1] < n) {
    ++i;
  }

  const auto x0 = static_cast<double>(xs[i]);
  const auto x1 = static_cast<double>(xs[i + 1]);
  const auto t  = (nd - x0) / (x1 - x0);

  return ((1.0 - t) * coeff(i) + t * coeff(i + 1)) / (nd * nd);
}
constexpr auto L1(Index n) { return interp_n2(ns, L1s, n); }

}  // namespace Expected

// =================================================================================================
constexpr auto u_analytical(Float y) -> Float {
  constexpr Float H = y_max - y_min;
  const Float s     = (y - y_min) / H;
  return Uavg * 6.0 * s * (1.0 - s);
}

// =================================================================================================
template <typename Float, Layout LAYOUT>
void correct_outflow(const Grid<Float, LAYOUT>& grid, FaceVector<Float, LAYOUT> u) {
  // Rescale the outflow so that it matches the inflow exactly (global continuity).
  Float Qin  = 0.0;
  Float Qout = 0.0;
  for (Index j = 0; j < u.x.ny(); ++j) {
    Qin  += u.x(0, j) * grid.dy();
    Qout += u.x(u.x.nx() - 1, j) * grid.dy();
  }
  const Float corr = (Qin - Qout) / (static_cast<Float>(u.x.ny()) * grid.dy());
  for (Index j = 0; j < u.x.ny(); ++j) {
    u.x(u.x.nx() - 1, j) += corr;
  }
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

#if FFT_POISSON
  const std::string output_dir = "./test/output/Channel-FFT-" + std::to_string(N) + '/';
#else
  const std::string output_dir = "./test/output/Channel-MG-" + std::to_string(N) + '/';
#endif
  if (!init_output_directory(output_dir)) { return 1; }

  Grid<Float, Layout::C> grid(x_min, x_max, ratio * N, y_min, y_max, N, 1);

  auto u_old = grid.alloc_face_vector();
  auto u     = grid.alloc_face_vector();

  auto FUX   = grid.alloc_scalar();
  auto FUY   = grid.alloc_vertex_scalar();
  auto FVX   = grid.alloc_vertex_scalar();
  auto FVY   = grid.alloc_scalar();

  auto ui    = grid.alloc_vector();
  auto p     = grid.alloc_scalar();
  auto dp    = grid.alloc_scalar();
  auto div   = grid.alloc_scalar();

  auto u_exp = grid.alloc_vector();
  grid.foreach_a(FOREACH_FUNC {
    u_exp.x(i, j) = u_analytical(grid.ym(j));
    u_exp.y(i, j) = 0.0;
  });

  Float dt = 0.0;
  Float t  = 0.0;

#if FFT_POISSON
  const std::array<int, 2> ns   = {grid.nx(), grid.ny()};
  const std::array<Float, 2> Ls = {grid.x_max() - grid.x_min(), grid.y_max() - grid.y_min()};
  const std::array<int, 4> BCs  = {
      PoisFFT::NEUMANN_STAG, PoisFFT::NEUMANN_STAG, PoisFFT::NEUMANN_STAG, PoisFFT::NEUMANN_STAG};
  PoisFFT::Solver<2, Float> solver(ns.data(), Ls.data(), BCs.data(), PoisFFT::FINITE_DIFFERENCE_2);
  const std::array<int, 2> ngs = {grid.nghost(), grid.nghost()};
#else
  MultigridSolver solver(grid);
  Float mg_res    = 0.0;
  Index mg_cycles = 0;
#endif

  const BConds<Float> u_bconds{
      .left   = Dirichlet<Float>{.val = [](Float y, Float /*t*/) { return u_analytical(y); }},
      .right  = Neumann{},
      .bottom = Dirichlet<Float>{.val = 0.0},
      .top    = Dirichlet<Float>{.val = 0.0},
  };

  const BConds<Float> v_bconds{
      .left   = Dirichlet<Float>{.val = 0.0},
      .right  = Neumann{},
      .bottom = Dirichlet<Float>{.val = 0.0},
      .top    = Dirichlet<Float>{.val = 0.0},
  };

  grid.foreach_face_i<Dimension::X>(FOREACH_FUNC { u.x(i, j) = 0.0; });
  grid.foreach_face_i<Dimension::Y>(FOREACH_FUNC { u.y(i, j) = 0.0; });
  apply_velocity_bconds(grid, u_bconds, v_bconds, u);
  interpolate(grid, u, ui);

  VTKWriter writer(output_dir, grid);
  writer.add_field("u", ui);
  writer.add_field("p", p);
  writer.add_field("div", div);
  writer.add_field("ua", u_exp);
  if (!writer.write(t)) { return 1; }

  Stats p_stats   = stats(grid, p);
  Stats u_stats   = stats(grid, u.x);
  Stats v_stats   = stats(grid, u.y);
  Stats div_stats = stats(grid, div);
  Float div_max   = std::max(std::abs(div_stats.min), std::abs(div_stats.max));

  Monitor<Float> monitor(output_dir + "/monitor.log");
  monitor.add_variable(&t, "t");
  monitor.add_variable(&dt, "dt");
  monitor.add_variable(&p_stats.max, "max(p)");
  monitor.add_variable(&u_stats.max, "max(u)");
  monitor.add_variable(&v_stats.max, "max(v)");
  monitor.add_variable(&div_max, "absmax(div)");
#if MG_POISSON
  monitor.add_variable(&mg_res, "residual(MG)");
  monitor.add_variable(&mg_cycles, "cycles(MG)");
#endif
  monitor.write();

  IGOR_TIME_SCOPE("Solver")
  while (t < tend) {
    // Time-step size: convective CFL and the (2D) explicit-diffusion limit dt <= h^2 / (4 nu).
    dt = adjust_dt(grid, u, rho, mu, CFL);
    dt = std::min({dt, dt_write, tend - t});

    copy(u, u_old);

    for (Index sub_iter = 0; sub_iter < 2; ++sub_iter) {
      const auto local_dt = sub_iter == 0 ? dt / 2.0 : dt;

      // 1) Predictor
      calc_flux(grid, u, p, rho, mu, FUX, FUY, FVX, FVY);
      update_u(grid, local_dt, FUX, FUY, FVX, FVY, u_old, u);
      apply_velocity_bconds(grid, u_bconds, v_bconds, u);
      correct_outflow(grid, u);

      // 2) Pressure correction
      calc_div(grid, u, div);
      grid.foreach_i(FOREACH_FUNC { div(i, j) *= rho / local_dt; });
#if FFT_POISSON
      solver.execute(dp.data(), div.data(), ngs.data(), ngs.data());
#else
      solver.solve(dp, div);
      mg_res    = solver.res();
      mg_cycles = solver.num_cycles();
#endif
      apply_neumann_bconds(grid, dp);

      // 3) Project
      correct_velocity(grid, dp, rho, local_dt, u, p);
    }

    interpolate(grid, u, ui);
    calc_div(grid, u, div);

    p_stats    = stats(grid, p);
    u_stats    = stats(grid, u.x);
    v_stats    = stats(grid, u.y);
    div_stats  = stats(grid, div);
    div_max    = std::max(std::abs(div_stats.min), std::abs(div_stats.max));

    t         += dt;
    if (should_save(t, dt, dt_write, tend)) {
      if (!writer.write(t)) { return 1; }
    }
    monitor.write();
  }

  Float L1 = 0.0;
  grid.foreach_range<Exec::SERIAL>(
      u.x.nx() / 2, u.x.nx() / 2 + 1, 0, u.x.ny(), [=, &L1](Index i, Index j) {
        const auto u_exp  = u_analytical(grid.ym(j));
        L1               += std::abs(u_exp - u.x(i, j)) * grid.dy();
      });
  if (L1 > 1.1 * Expected::L1(N)) {
    Igor::Error("u error does not match expected value: expected {:.8e} but got {:.8e}",
                Expected::L1(N),
                L1);
    return 1;
  }
}
