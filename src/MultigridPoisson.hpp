#pragma once

#include <Igor/Math.hpp>

#include "BoundaryConditions.hpp"
#include "Grid.hpp"

template <typename Float, Layout LAYOUT>
class MultigridSolver {
  using Grid   = Grid<Float, LAYOUT>;
  using Scalar = Scalar<Float, LAYOUT>;

  struct Level {
    Grid grid;
    Scalar sol;
    Scalar rhs;
    Scalar res;
  };
  std::vector<Level> m_levels;
  BConds<Float> m_bconds;
  Index m_num_iter_pre;
  Index m_num_iter_post;
  Index m_num_iter_coarse;
  Index m_max_iter_coarse;
  Index m_num_cycles                = 0;
  Float m_res                       = 0.0;

  static constexpr Float COARSE_TOL = 1e-3;

  // -----------------------------------------------------------------------------------------------
  constexpr void make_mean_free(const Grid& grid, Scalar s) const noexcept {
    Float mean = 0.0;
    Float vol  = 0.0;
    grid.template foreach_i<Exec::SERIAL>([=, &mean, &vol](Index i, Index j) {
      mean += s(i, j) * grid.dv(i, j);
      vol  += grid.dv(i, j);
    });
    mean /= vol;
    grid.foreach_i(FOREACH_FUNC { s(i, j) -= mean; });
  }

  // -----------------------------------------------------------------------------------------------
  constexpr auto residual(const Level& level) const noexcept -> Float {
    const Float inv_dx2 = 1.0 / Igor::sqr(level.grid.dx());
    const Float inv_dy  = 1.0 / level.grid.dy();
    const Float inv_dy2 = 1.0 / Igor::sqr(level.grid.dy());
    auto sol            = level.sol;
    auto rhs            = level.rhs;
    auto res            = level.res;

    // apply_neumann_bconds(level.grid, sol);
    apply_bconds(level.grid, m_bconds, sol, -1.0);
    std::atomic<Float> max_res = 0.0;
    switch (level.grid.coords()) {
      case Coordinates::CARTESIAN:
        level.grid.foreach_i([=, &max_res](Index i, Index j) {
          const Float L = (sol(i - 1, j) - 2.0 * sol(i, j) + sol(i + 1, j)) * inv_dx2 +
                          (sol(i, j - 1) - 2.0 * sol(i, j) + sol(i, j + 1)) * inv_dy2;
          res(i, j)     = rhs(i, j) - L;
          update_maximum_atomic(max_res, std::abs(res(i, j)));
        });
        break;
      case Coordinates::POLAR:
        level.grid.foreach_i([=, &max_res](Index i, Index j) {
          const Float dpdr     = (sol(i, j + 1) - sol(i, j - 1)) * 0.5 * inv_dy;
          const Float ddpdrr   = (sol(i, j - 1) - 2.0 * sol(i, j) + sol(i, j + 1)) * inv_dy2;
          const Float ddpdthth = (sol(i - 1, j) - 2.0 * sol(i, j) + sol(i + 1, j)) * inv_dx2;
          const Float r        = level.grid.ym(j);

          const Float L        = ddpdrr + dpdr / r + ddpdthth / Igor::sqr(r);
          res(i, j)            = rhs(i, j) - L;
          update_maximum_atomic(max_res, std::abs(res(i, j)));
        });
        break;
    }
    return static_cast<Float>(max_res);
  }

  // -----------------------------------------------------------------------------------------------
  constexpr void smooth_cartesian(const Level& level, Index num_iter) {
    const Float inv_dx2 = 1.0 / Igor::sqr(level.grid.dx());
    const Float inv_dy2 = 1.0 / Igor::sqr(level.grid.dy());
    const Float idiag   = 1.0 / (2.0 * (inv_dx2 + inv_dy2));
    auto sol            = level.sol;
    auto rhs            = level.rhs;

    // Do exactly num_iter iterations of the Gauss-Seidel algorithm
    for (Index iter = 0; iter < num_iter; ++iter) {
      // apply_neumann_bconds(level.grid, sol);
      apply_bconds(level.grid, m_bconds, sol, -1.0);
#ifndef NS_FVM_PARALLEL
      level.grid.template foreach_i<Exec::SERIAL>(FOREACH_FUNC {
        sol(i, j) = ((sol(i - 1, j) + sol(i + 1, j)) * inv_dx2 +  //
                     (sol(i, j - 1) + sol(i, j + 1)) * inv_dy2 -  //
                     rhs(i, j)) *
                    idiag;
      });
#else
      // Red-black Gauss-Seidel for parallel execution
      const auto nx_half = level.grid.nx() / 2;
      const auto nx_rem  = level.grid.nx() % 2;
      const auto nx      = level.grid.nx();
      // First pass
      level.grid.foreach_range(
          0, nx_half + nx_rem, 0, level.grid.ny(), FOREACH_FUNC {
            const Index nj = j;
            const Index ni = j % 2 == 0 ? 2 * i : 2 * i + 1;
            if (ni >= nx) { return; }
            sol(ni, nj) = ((sol(ni - 1, nj) + sol(ni + 1, nj)) * inv_dx2 +  //
                           (sol(ni, nj - 1) + sol(ni, nj + 1)) * inv_dy2 -  //
                           rhs(ni, nj)) *
                          idiag;
          });
      // Second pass
      level.grid.foreach_range(
          0, nx_half + nx_rem, 0, level.grid.ny(), FOREACH_FUNC {
            const Index nj = j;
            const Index ni = j % 2 == 1 ? 2 * i : 2 * i + 1;
            if (ni >= nx) { return; }
            sol(ni, nj) = ((sol(ni - 1, nj) + sol(ni + 1, nj)) * inv_dx2 +  //
                           (sol(ni, nj - 1) + sol(ni, nj + 1)) * inv_dy2 -  //
                           rhs(ni, nj)) *
                          idiag;
          });
#endif  // NS_FVM_PARALLEL
    }
  }

  constexpr void smooth_polar(const Level& level, Index num_iter) {
    const Float inv_dx2 = 1.0 / Igor::sqr(level.grid.dx());
    const Float inv_dy  = 1.0 / level.grid.dy();
    const Float inv_dy2 = 1.0 / Igor::sqr(level.grid.dy());
    auto sol            = level.sol;
    auto rhs            = level.rhs;

    // Do exactly num_iter iterations of the Gauss-Seidel algorithm
    for (Index iter = 0; iter < num_iter; ++iter) {
      // apply_neumann_bconds(level.grid, sol);
      apply_bconds(level.grid, m_bconds, sol, -1.0);
#ifndef NS_FVM_PARALLEL
      level.grid.template foreach_i<Exec::SERIAL>(FOREACH_FUNC {
        const Float r = level.grid.ym(j);
        sol(i, j)     = ((sol(i - 1, j) + sol(i + 1, j)) * inv_dx2 / Igor::sqr(r) +  //
                         (sol(i, j - 1) + sol(i, j + 1)) * inv_dy2 +                 //
                         (sol(i, j + 1) - sol(i, j - 1)) * 0.5 * inv_dy / r -        //
                         rhs(i, j)) /
                        (2.0 * inv_dx2 / Igor::sqr(r) + 2.0 * inv_dy2);
      });
#else
      // Red-black Gauss-Seidel for parallel execution
      const auto nx_half = level.grid.nx() / 2;
      const auto nx_rem  = level.grid.nx() % 2;
      const auto nx      = level.grid.nx();
      // First pass
      level.grid.foreach_range(
          0, nx_half + nx_rem, 0, level.grid.ny(), FOREACH_FUNC {
            const Index nj = j;
            const Index ni = j % 2 == 0 ? 2 * i : 2 * i + 1;
            if (ni >= nx) { return; }
            const Float r = level.grid.ym(nj);
            sol(ni, nj)   = ((sol(ni - 1, nj) + sol(ni + 1, nj)) * inv_dx2 / Igor::sqr(r) +  //
                             (sol(ni, nj - 1) + sol(ni, nj + 1)) * inv_dy2 +                 //
                             (sol(ni, nj + 1) - sol(ni, nj - 1)) * 0.5 * inv_dy / r -        //
                             rhs(ni, nj)) /
                            (2.0 * inv_dx2 / Igor::sqr(r) + 2.0 * inv_dy2);
          });
      // Second pass
      level.grid.foreach_range(
          0, nx_half + nx_rem, 0, level.grid.ny(), FOREACH_FUNC {
            const Index nj = j;
            const Index ni = j % 2 == 1 ? 2 * i : 2 * i + 1;
            if (ni >= nx) { return; }
            const Float r = level.grid.ym(nj);
            sol(ni, nj)   = ((sol(ni - 1, nj) + sol(ni + 1, nj)) * inv_dx2 / Igor::sqr(r) +  //
                             (sol(ni, nj - 1) + sol(ni, nj + 1)) * inv_dy2 +                 //
                             (sol(ni, nj + 1) - sol(ni, nj - 1)) * 0.5 * inv_dy / r -        //
                             rhs(ni, nj)) /
                            (2.0 * inv_dx2 / Igor::sqr(r) + 2.0 * inv_dy2);
          });
#endif  // NS_FVM_PARALLEL
    }
  }

  constexpr void smooth(const Level& level, Index num_iter) {
    switch (level.grid.coords()) {
      case Coordinates::CARTESIAN: return smooth_cartesian(level, num_iter);
      case Coordinates::POLAR:     return smooth_polar(level, num_iter);
    }
    Igor::Panic("Unreachable");
  }

  // -----------------------------------------------------------------------------------------------
  constexpr void restrict_residual(const Level& level, const Level& coarse) {
    IGOR_ASSERT(level.grid.nx() / 2 == coarse.grid.nx() && level.grid.ny() / 2 == coarse.grid.ny(),
                "Expected `coarse` to be the next coarser level but we skipped something.");
    auto res         = level.res;
    auto rhs         = coarse.rhs;
    const auto fgrid = level.grid;

    // Residual of `level` becomes the rhs of `coarse`. The average must be weighted by the cell
    // volumes, otherwise the restriction is not conservative on non-uniform volumes (polar).
    coarse.grid.foreach_i(FOREACH_FUNC {
      const Float wb = fgrid.dv(2 * i, 2 * j);
      const Float wt = fgrid.dv(2 * i, 2 * j + 1);
      rhs(i, j)      = (wb * (res(2 * i, 2 * j) + res(2 * i + 1, 2 * j)) +
                        wt * (res(2 * i, 2 * j + 1) + res(2 * i + 1, 2 * j + 1))) /
                       (2.0 * (wb + wt));
    });
    make_mean_free(coarse.grid, rhs);
  }

  // -----------------------------------------------------------------------------------------------
  constexpr void prolongate_and_correct(const Level& coarse, const Level& level) {
    IGOR_ASSERT(level.grid.nx() / 2 == coarse.grid.nx() && level.grid.ny() / 2 == coarse.grid.ny(),
                "Expected `coarse` to be the next coarser level but we skipped something.");
    auto lsol = level.sol;
    auto csol = coarse.sol;

    // apply_neumann_bconds(coarse.grid, coarse.sol);
    apply_bconds(coarse.grid, m_bconds, coarse.sol, -1.0);
    // Bilinear interpolation of the coarse correction onto the finer solution
    level.grid.foreach_i(FOREACH_FUNC {
      const Index ic  = i / 2;
      const Index jc  = j / 2;
      const Index di  = 2 * (i % 2) - 1;
      const Index dj  = 2 * (j % 2) - 1;
      lsol(i, j)     += (9.0 * csol(ic, jc) + 3.0 * (csol(ic + di, jc) + csol(ic, jc + dj)) +
                         csol(ic + di, jc + dj)) /
                        16.0;
    });
  }

  // -----------------------------------------------------------------------------------------------
  constexpr void vcycle(size_t l) {
    IGOR_ASSERT(l < m_levels.size(), "Level {} is out of bounds for {} levels", l, num_levels());
    Level& level = m_levels[l];

    // Check if we are at the coarsest level
    if (l + 1 == m_levels.size()) {
      for (Index iter = 0; iter < m_max_iter_coarse; iter += m_num_iter_coarse) {
        smooth(level, m_num_iter_coarse);
        if (residual(level) <= COARSE_TOL) { break; }
      }
      return;
    }

    Level& coarse = m_levels[l + 1];
    smooth(level, m_num_iter_pre);  // Do some iterations to remove high frequencies from residual
    residual(level);                // Calculate the residual to be used in correction equation

    restrict_residual(level, coarse);  // Interpolate the residual of `level` onto `rhs` of coarse
    fill(coarse.sol, 0.0);
    vcycle(l + 1);  // Solve correction equation for `coarse`
    make_mean_free(coarse.grid, coarse.sol);

    // Bilinear interpolation of the coarse correction onto level
    prolongate_and_correct(coarse, level);
    smooth(level, m_num_iter_post);
  }

 public:
  constexpr MultigridSolver(const Grid& grid,
                            BConds<Float> bconds  = {.left   = Neumann{},
                                                     .right  = Neumann{},
                                                     .bottom = Neumann{},
                                                     .top    = Neumann{}},
                            Index min_size        = 2,
                            Index num_iter_pre    = 2,
                            Index num_iter_post   = 2,
                            Index num_iter_coarse = 50,
                            Index max_iter_coarse = 5000)
      : m_bconds(std::move(bconds)),
        m_num_iter_pre(num_iter_pre),
        m_num_iter_post(num_iter_post),
        m_num_iter_coarse(num_iter_coarse),
        m_max_iter_coarse(max_iter_coarse) {
    IGOR_ASSERT(grid.nghost() >= 1, "Expected at least one ghost cell, but got {}", grid.nghost());
    IGOR_ASSERT(min_size >= 1, "Expected a positive minimum grid size, but got {}", min_size);

    Grid level_grid(grid.x_min(),
                    grid.x_max(),
                    grid.nx(),
                    grid.y_min(),
                    grid.y_max(),
                    grid.ny(),
                    grid.nghost(),
                    grid.coords());
    while (true) {
      m_levels.emplace_back(level_grid,
                            level_grid.alloc_scalar(),
                            level_grid.alloc_scalar(),
                            level_grid.alloc_scalar());

      if (level_grid.nx() % 2 != 0 || level_grid.ny() % 2 != 0 || level_grid.nx() / 2 < min_size ||
          level_grid.ny() / 2 < min_size) {
        break;
      }
      level_grid = Grid(level_grid.x_min(),
                        level_grid.x_max(),
                        level_grid.nx() / 2,
                        level_grid.y_min(),
                        level_grid.y_max(),
                        level_grid.ny() / 2,
                        level_grid.nghost(),
                        level_grid.coords());
    }
  }

  constexpr auto solve(Scalar sol, Scalar rhs, Float tol = 1e-4, Index max_iter = 100) -> bool {
    const Level& fine = m_levels[0];
    IGOR_ASSERT(sol.nx() == fine.sol.nx() && sol.ny() == fine.sol.ny() &&
                    sol.nghost() == fine.sol.nghost(),
                "Field `sol` does not match the grid the solver was constructed with.");
    IGOR_ASSERT(rhs.nx() == fine.rhs.nx() && rhs.ny() == fine.rhs.ny() &&
                    rhs.nghost() == fine.rhs.nghost(),
                "Field `rhs` does not match the grid the solver was constructed with.");

    copy(sol, fine.sol);
    copy(rhs, fine.rhs);
    make_mean_free(fine.grid, fine.rhs);

    bool converged = false;
    for (m_num_cycles = 0; true; ++m_num_cycles) {
      m_res = residual(fine);
      if (m_res <= tol) {
        converged = true;
        break;
      }
      if (m_num_cycles == max_iter) { break; }
      vcycle(0);
    }

    make_mean_free(fine.grid, fine.sol);
    copy(fine.sol, sol);

    return converged;
  }

  [[nodiscard]] constexpr auto num_levels() const noexcept -> Index {
    return static_cast<Index>(m_levels.size());
  }
  [[nodiscard]] constexpr auto num_cycles() const noexcept -> Index { return m_num_cycles; }
  [[nodiscard]] constexpr auto res() const noexcept -> Float { return m_res; }
};
