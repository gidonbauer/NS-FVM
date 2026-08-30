#include <charconv>
#include <cmath>
#include <numbers>

#include <poisfft.h>

#include <Igor/Logging.hpp>
#include <Igor/Math.hpp>
#include <Igor/Timer.hpp>

#include "BoundaryConditions.hpp"
#include "Grid.hpp"
#include "IO.hpp"
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
constexpr auto poisson_gauss_seidel(const Grid<Float, LAYOUT>& grid,
                                    Scalar<Float, LAYOUT> x,
                                    const Scalar<Float, LAYOUT> rhs,
                                    Float tol      = 1e-6,
                                    Index max_iter = 100) -> bool {
  const Float idx2  = 1.0 / Igor::sqr(grid.dx());
  const Float idy2  = 1.0 / Igor::sqr(grid.dy());
  const Float idiag = 1.0 / (2.0 * (idx2 + idy2));

  bool converged    = false;
  for (Index iter = 0; iter <= max_iter; ++iter) {
    apply_neumann_bconds(grid, x);

    // Residual r = rhs - laplace(x) in max-norm
    Float res = 0.0;
    grid.template foreach_i<Exec::SERIAL>([=, &res](Index i, Index j) {
      const Float laplace_x = (x(i - 1, j) - Float{2} * x(i, j) + x(i + 1, j)) * idx2 +
                              (x(i, j - 1) - Float{2} * x(i, j) + x(i, j + 1)) * idy2;
      res                   = std::max(res, Igor::abs(rhs(i, j) - laplace_x));
    });
    if (res <= tol) {
      converged = true;
      break;
    }
    if (iter == max_iter) { break; }

    // Gauss-Seidel sweep
    grid.template foreach_i<Exec::SERIAL>([=](Index i, Index j) {
      x(i, j) =
          ((x(i - 1, j) + x(i + 1, j)) * idx2 + (x(i, j - 1) + x(i, j + 1)) * idy2 - rhs(i, j)) *
          idiag;
    });
  }

  // Make the solution mean-free
  Float mean = 0.0;
  grid.template foreach_i<Exec::SERIAL>([=, &mean](Index i, Index j) { mean += x(i, j); });
  mean /= static_cast<Float>(grid.nx() * grid.ny());
  grid.foreach_i(FOREACH_FUNC { x(i, j) -= mean; });

  return converged;
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

  auto f_exp = grid.alloc_scalar();
  auto f     = grid.alloc_scalar();
  auto rhs   = grid.alloc_scalar();
  auto f_gs  = grid.alloc_scalar();

  grid.foreach_i(FOREACH_FUNC { rhs(i, j) = ddF(grid.xm(i), grid.ym(j)); });
  grid.foreach_i(FOREACH_FUNC { f_exp(i, j) = F(grid.xm(i), grid.ym(j)); });

  const std::array<int, 2> ns   = {grid.nx(), grid.ny()};
  const std::array<Float, 2> Ls = {grid.x_max() - grid.x_min(), grid.y_max() - grid.y_min()};
  const std::array<int, 4> BCs  = {
      PoisFFT::NEUMANN_STAG, PoisFFT::NEUMANN_STAG, PoisFFT::NEUMANN_STAG, PoisFFT::NEUMANN_STAG};
  const std::array<int, 2> ngs = {grid.nghost(), grid.nghost()};

  {
    Igor::ScopeTimer timer("FINITE_DIFFERENCE_2");
    PoisFFT::Solver<2, Float> solver(
        ns.data(), Ls.data(), BCs.data(), PoisFFT::FINITE_DIFFERENCE_2);
    solver.execute(f.data(), rhs.data(), ngs.data(), ngs.data());
    Float L1 = 0.0;
    grid.foreach_i<Exec::SERIAL>(
        [=, &L1](Index i, Index j) { L1 += std::abs(f_exp(i, j) - f(i, j)) * grid.dv(); });
    Igor::Info("FINITE_DIFFERENCE_2");
    Igor::Info("  L1({}) = {:.8e}", N, L1);
  }

  {
    Igor::ScopeTimer timer("FINITE_DIFFERENCE_4");
    PoisFFT::Solver<2, Float> solver(ns.data(), Ls.data(), BCs.data(), 4);
    solver.execute(f.data(), rhs.data(), ngs.data(), ngs.data());
    Float L1 = 0.0;
    grid.foreach_i<Exec::SERIAL>(
        [=, &L1](Index i, Index j) { L1 += std::abs(f_exp(i, j) - f(i, j)) * grid.dv(); });
    Igor::Info("FINITE_DIFFERENCE_4");
    Igor::Info("  L1({}) = {:.8e}", N, L1);
  }

  {
    Igor::ScopeTimer timer("SPECTRAL");
    PoisFFT::Solver<2, Float> solver(ns.data(), Ls.data(), BCs.data(), PoisFFT::SPECTRAL);
    solver.execute(f.data(), rhs.data(), ngs.data(), ngs.data());
    Float L1 = 0.0;
    grid.foreach_i<Exec::SERIAL>(
        [=, &L1](Index i, Index j) { L1 += std::abs(f_exp(i, j) - f(i, j)) * grid.dv(); });
    Igor::Info("SPECTRAL");
    Igor::Info("  L1({}) = {:.8e}", N, L1);
  }

  bool gs_converged = false;
  {
    Igor::ScopeTimer timer("GAUSS_SEIDEL");
    gs_converged = poisson_gauss_seidel(grid, f_gs, rhs, 1e-6, 100'000);
    Float L1     = 0.0;
    grid.foreach_i<Exec::SERIAL>(
        [=, &L1](Index i, Index j) { L1 += std::abs(f_exp(i, j) - f_gs(i, j)) * grid.dv(); });
    Igor::Info("GAUSS_SEIDEL");
    if (!gs_converged) { Igor::Warn("  Gauss-Seidel did not converge."); }
    Igor::Info("  L1({}) = {:.8e}", N, L1);
  }

  const auto output_dir = get_output_directory("test/output");
  if (!init_output_directory(output_dir)) { return 1; }
  VTKWriter writer(output_dir, grid);
  writer.add_field("f", f);
  writer.add_field("f_exp", f_exp);
  writer.add_field("f_gs", f_gs);
  if (!writer.write(0.0)) { return 1; }

  return gs_converged ? 0 : 1;
}
