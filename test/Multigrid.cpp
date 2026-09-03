#include <charconv>
#include <cmath>
#include <numbers>

#include <poisfft.h>

#include <Igor/Logging.hpp>
#include <Igor/Math.hpp>
#include <Igor/Timer.hpp>

#include "Grid.hpp"
#include "IO.hpp"
#include "MultigridPoisson.hpp"
#include "VTKWriter.hpp"

using Float           = double;
constexpr auto pi     = std::numbers::pi_v<Float>;
constexpr Float x_min = 0.0;
constexpr Float x_max = 1.0;
constexpr Float y_min = 0.0;
constexpr Float y_max = 1.0;

// f(x)  = cos(2pi x) + cos(2pi y)
// f'(x) = -2pi sin(2pi x) - 2pi sin(2pi y)
// f"(x) = -4pi^2 cos(2pi x) - 4pi^2 cos(2pi y)
constexpr auto F(Float x, Float y) -> Float {
  return std::cos(2.0 * pi * x) + std::cos(2.0 * pi * y);
}
constexpr auto ddF(Float x, Float y) -> Float {
  return -4.0 * Igor::sqr(pi) * (std::cos(2.0 * pi * x) + std::cos(2.0 * pi * y));
}

template <typename Float, Layout LAYOUT>
constexpr auto L1error(const Grid<Float, LAYOUT>& grid,
                       const Scalar<Float, LAYOUT> f_true,
                       const Scalar<Float, LAYOUT> f_pred) -> Float {
  Float L1 = 0.0;
  grid.template foreach_i<Exec::SERIAL>(
      [=, &L1](Index i, Index j) { L1 += std::abs(f_true(i, j) - f_pred(i, j)) * grid.dv(i, j); });
  return L1;
}

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
  Grid<Float> grid(x_min, x_max, N, y_min, y_max, N, 1);

  auto rhs   = grid.alloc_scalar();
  auto f_exp = grid.alloc_scalar();
  auto f_fft = grid.alloc_scalar();
  auto f_gs  = grid.alloc_scalar();
  auto f_mg  = grid.alloc_scalar();

  grid.foreach_i(FOREACH_FUNC { rhs(i, j) = ddF(grid.xm(i), grid.ym(j)); });
  grid.foreach_i(FOREACH_FUNC { f_exp(i, j) = F(grid.xm(i), grid.ym(j)); });

  Float L1_fft = -1.0;
  {
    Igor::ScopeTimer timer("FFT-" + std::to_string(N));
    const std::array<int, 2> ns   = {grid.nx(), grid.ny()};
    const std::array<Float, 2> Ls = {grid.x_max() - grid.x_min(), grid.y_max() - grid.y_min()};
    const std::array<int, 4> BCs  = {
        PoisFFT::NEUMANN_STAG, PoisFFT::NEUMANN_STAG, PoisFFT::NEUMANN_STAG, PoisFFT::NEUMANN_STAG};
    const std::array<int, 2> ngs = {grid.nghost(), grid.nghost()};
    PoisFFT::Solver<2, Float> solver(
        ns.data(), Ls.data(), BCs.data(), PoisFFT::FINITE_DIFFERENCE_2);
    solver.execute(f_fft.data(), rhs.data(), ngs.data(), ngs.data());
    L1_fft = L1error(grid, f_exp, f_fft);
  }

  bool mg_converged = false;
  Float L1_mg       = -1.0;
  {
    Igor::ScopeTimer timer("Multigrid-" + std::to_string(N));
    MultigridSolver solver(grid);
    mg_converged = solver.solve(f_mg, rhs, 1e-6, 100'000);
    L1_mg        = L1error(grid, f_exp, f_mg);
    if (!mg_converged) { Igor::Error("Multigrid solver did not converge."); }
  }

  if (L1_mg > 1.1 * L1_fft) {
    Igor::Error("Multigrid solver did not converge to the same solution as the FFT solver: L1_fft "
                "= {:.12e}, L1_mg = {:.12e}",
                L1_fft,
                L1_mg);
    mg_converged = false;
  }

  // Igor::Debug("L1_fft = {:.8e}", L1_fft);
  // Igor::Debug("L1_mg  = {:.8e}", L1_mg);

  const auto output_dir = "./test/output/Multigrid-" + std::to_string(N);
  if (!init_output_directory(output_dir)) { return 1; }
  VTKWriter writer(output_dir, grid);
  writer.add_field("f_exp", f_exp);
  writer.add_field("f_fft", f_fft);
  writer.add_field("f_gs", f_gs);
  writer.add_field("f_mg", f_mg);
  if (!writer.write(0.0)) { return 1; }

  return mg_converged ? 0 : 1;
}
