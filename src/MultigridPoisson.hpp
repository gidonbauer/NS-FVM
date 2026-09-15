#pragma once

#include <Igor/Math.hpp>

#include "BoundaryConditions.hpp"
#include "Grid.hpp"

template <typename Float, Layout LAYOUT>
class MultigridSolver {
  using Grid                               = Grid<Float, LAYOUT>;
  using Scalar                             = Scalar<Float, LAYOUT>;

  static constexpr Index MIN_NUM_ITER_POST = 2;
  static constexpr Index MAX_NUM_ITER_POST = 100;
  static constexpr Index MIN_NUM_ITER_PRE  = 1;
  static constexpr Index MAX_NUM_ITER_PRE  = 100;
  static constexpr Float OMEGA             = 1.2;

  struct Level {
    Grid grid;
    Scalar sol;
    Scalar rhs;
    Scalar res;

    // Restriction weights; depend on cell volume
    std::vector<Float> restrict_w_lo;
    std::vector<Float> restrict_w_hi;

    // = Relevant only for polar coordinates ===================================
    // Thomas coefficients for linear solver in r-direction in polar case
    std::vector<Float> tri_a;      // Sub-diagonal, zero at j = 0
    std::vector<Float> tri_cstar;  // Super-diagonal divided by the pivot
    std::vector<Float> tri_minv;   // Reciprocal of the pivot
    Float tri_bnd_lo = 0.0;  // Coupling to the ghost below, zero if it folds into the diagonal
    Float tri_bnd_hi = 0.0;  // Coupling to the ghost above, zero if it folds into the diagonal

    // Polar stencil coefficients
    std::vector<Float> tri_theta;  // (i, j)
    std::vector<Float> tri_lo;     // (i, j - 1)
    std::vector<Float> tri_hi;     // (i, j + 1)
    // = Relevant only for polar coordinates ===================================
  };
  std::vector<Level> m_levels;
  BConds<Float> m_bconds;
  Index m_num_iter_pre;
  Index m_num_iter_post;
  Index m_num_cycles = 0;
  Float m_res        = 0.0;

  // -----------------------------------------------------------------------------------------------
  static constexpr void precompute_restriction_weights(const Level& fine, Level& coarse) {
    const Index cny = coarse.grid.ny();
    coarse.restrict_w_lo.resize(static_cast<size_t>(cny));
    coarse.restrict_w_hi.resize(static_cast<size_t>(cny));

    for (Index jc = 0; jc < cny; ++jc) {
      const Float w_lo                              = fine.grid.dv(0, 2 * jc);
      const Float w_hi                              = fine.grid.dv(0, 2 * jc + 1);
      const Float norm                              = 1.0 / (2.0 * (w_lo + w_hi));

      coarse.restrict_w_lo[static_cast<size_t>(jc)] = w_lo * norm;
      coarse.restrict_w_hi[static_cast<size_t>(jc)] = w_hi * norm;
    }
  }

  // -----------------------------------------------------------------------------------------------
  // Precompute the coefficients for the Thomas algorithm in r-direction (y-direction);
  // Only for polar coordinates
  static constexpr void precompute_tridiag(Level& level, const BConds<Float>& bconds) {
    if (level.grid.coords() != Coordinates::POLAR) { return; }

    const Grid& grid    = level.grid;
    const Index n       = grid.ny();
    const Float inv_dx2 = 1.0 / Igor::sqr(grid.dx());
    const Float inv_dy  = 1.0 / grid.dy();
    const Float inv_dy2 = 1.0 / Igor::sqr(grid.dy());

    level.tri_a.resize(static_cast<size_t>(n));
    level.tri_cstar.resize(static_cast<size_t>(n));
    level.tri_minv.resize(static_cast<size_t>(n));

    level.tri_theta.resize(static_cast<size_t>(n));
    level.tri_lo.resize(static_cast<size_t>(n));
    level.tri_hi.resize(static_cast<size_t>(n));

    // Account for periodic boundary conditions
    const bool fold_lo = std::holds_alternative<Neumann>(bconds.bottom);
    const bool fold_hi = std::holds_alternative<Neumann>(bconds.top);

    Float cstar_prev   = 0.0;
    for (Index j = 0; j < n; ++j) {
      const Float r = grid.ym(j);
      const Float a = inv_dy2 - 0.5 * inv_dy / r;  // Coefficient of lower diagonal
      const Float c = inv_dy2 + 0.5 * inv_dy / r;  // Coefficient of upper diagonal
      Float b       = -2.0 * inv_dy2 - 2.0 * inv_dx2 / Igor::sqr(r);  // Coefficient of diagonal

      if (j == 0) {
        if (fold_lo) {
          b += a;  // Adjustment for Neumann boundary condition
        } else {
          level.tri_bnd_lo = a;
        }
      }
      if (j == n - 1) {
        if (fold_hi) {
          b += c;  // Adjustment for Neumann boundary condition
        } else {
          level.tri_bnd_hi = c;
        }
      }

      // The sub-diagonal is zero in first row
      const Float a_in = j == 0 ? 0.0 : a;
      // Pivot: b_j if j = 0, b_j - a_j*cstar_(j-1) otherwise
      const Float m                           = b - a_in * cstar_prev;

      level.tri_a[static_cast<size_t>(j)]     = a_in;
      level.tri_minv[static_cast<size_t>(j)]  = 1.0 / m;
      level.tri_cstar[static_cast<size_t>(j)] = c / m;
      level.tri_theta[static_cast<size_t>(j)] = inv_dx2 / Igor::sqr(r);
      level.tri_lo[static_cast<size_t>(j)]    = a;
      level.tri_hi[static_cast<size_t>(j)]    = c;
      cstar_prev                              = level.tri_cstar[static_cast<size_t>(j)];
    }
  }

  // -----------------------------------------------------------------------------------------------
  constexpr void make_mean_free(const Grid& grid, Scalar s) const noexcept {
    struct SumVol {
      Float sum, vol;
    };
    const SumVol sv = grid.transform_reduce_i(
        SumVol{.sum = 0.0, .vol = 0.0},
        FOREACH_FUNC { return SumVol{.sum = s(i, j) * grid.dv(i, j), .vol = grid.dv(i, j)}; },
        [](SumVol lhs, const SumVol& rhs) {
          lhs.sum += rhs.sum;
          lhs.vol += rhs.vol;
          return lhs;
        });
    const auto mean = sv.sum / sv.vol;
    grid.foreach_i(FOREACH_FUNC { s(i, j) -= mean; });
  }

  // -----------------------------------------------------------------------------------------------
  constexpr void residual(const Level& level) const noexcept {
    const Float inv_dx2 = 1.0 / Igor::sqr(level.grid.dx());
    const Float inv_dy2 = 1.0 / Igor::sqr(level.grid.dy());
    auto sol            = level.sol;
    auto rhs            = level.rhs;
    auto res            = level.res;

    apply_bconds(level.grid, m_bconds, sol, -1.0);
    switch (level.grid.coords()) {
      case Coordinates::CARTESIAN:
        level.grid.foreach_i(FOREACH_FUNC {
          const Float c = sol(i, j);
          const Float L = (sol(i - 1, j) - 2.0 * c + sol(i + 1, j)) * inv_dx2 +
                          (sol(i, j - 1) - 2.0 * c + sol(i, j + 1)) * inv_dy2;
          res(i, j)     = rhs(i, j) - L;
        });
        break;
      case Coordinates::POLAR:
        {
          // Get precomputed coefficients
          const Float* theta = level.tri_theta.data();  // inv_dx2/sq(r)
          const Float* c_lo  = level.tri_lo.data();     // inv_dy2 - 0.5*inv_dy/r
          const Float* c_hi  = level.tri_hi.data();     // inv_dy2 + 0.5*inv_dy/r

          level.grid.foreach_i(FOREACH_FUNC {
            const Float c = sol(i, j);
            const Float L =
                // d^2(sol)/dr^2 + 1/r*d(sol)/dr
                c_lo[j] * sol(i, j - 1) + c_hi[j] * sol(i, j + 1) + -2.0 * inv_dy2 * c +
                // 1/r^2*d^2(sol)/d(theta)^2
                theta[j] * (sol(i - 1, j) - 2.0 * c + sol(i + 1, j));
            res(i, j) = rhs(i, j) - L;
          });
        }
        break;
    }
  }

  // -----------------------------------------------------------------------------------------------
  constexpr auto max_res(const Level& level) const noexcept -> Float {
    const auto res = level.res;
    return level.grid.transform_reduce_i(
        0.0,
        FOREACH_FUNC { return std::abs(res(i, j)); },
        [](Float lhs, Float rhs) { return std::max(lhs, rhs); });
  }

  // -----------------------------------------------------------------------------------------------
  constexpr void smooth_cartesian(const Level& level, Index num_iter) {
    const Float inv_dx2 = 1.0 / Igor::sqr(level.grid.dx());
    const Float inv_dy2 = 1.0 / Igor::sqr(level.grid.dy());
    const Float idiag   = 1.0 / (2.0 * (inv_dx2 + inv_dy2));
    const auto sol      = level.sol;
    const auto rhs      = level.rhs;
    const Index nx      = level.grid.nx();
    const Index ny      = level.grid.ny();

    for (Index iter = 0; iter < num_iter; ++iter) {
      // Red-black Gauss-Seidel with over-relaxation
      for (Index parity = 0; parity < 2; ++parity) {
        apply_bconds(level.grid, m_bconds, sol, -1.0);

        level.grid.foreach_range(0, nx, 0, (ny + 1) / 2, [=](Index i, Index jj) {
          const Index j = 2 * jj + (i + parity) % 2;
          if (j >= ny) { return; }

          sol(i, j) = (1.0 - OMEGA) * sol(i, j) +
                      OMEGA *
                          ((sol(i - 1, j) + sol(i + 1, j)) * inv_dx2 +
                           (sol(i, j - 1) + sol(i, j + 1)) * inv_dy2 - rhs(i, j)) *
                          idiag;
        });
      }
    }
  }

  // -----------------------------------------------------------------------------------------------
  constexpr void smooth_polar(const Level& level, Index num_iter) {
    const auto sol     = level.sol;
    const auto rhs     = level.rhs;
    const Index nx     = level.grid.nx();
    const Index ny     = level.grid.ny();

    const Index sol_si = sol.stride_x();
    const Index sol_sj = sol.stride_y();
    const Index rhs_sj = rhs.stride_y();

    const Float* tri_a = level.tri_a.data();
    const Float* cstar = level.tri_cstar.data();
    const Float* minv  = level.tri_minv.data();
    const Float* theta = level.tri_theta.data();
    const Float bnd_lo = level.tri_bnd_lo;
    const Float bnd_hi = level.tri_bnd_hi;

    for (Index iter = 0; iter < num_iter; ++iter) {
      apply_bconds(level.grid, m_bconds, sol, -1.0);

      // Zebra line relaxation along r:
      // - Solve linear system for an entire row in r-direction (y-direction) -> Thomas algorithm
      // - Skip every second line for parallelization
      for (Index parity = 0; parity < 2; ++parity) {
        level.grid.foreach_range(0, (nx + 1 - parity) / 2, 0, 1, [=](Index ii, Index /*unused*/) {
          const Index i = 2 * ii + parity;

          // Use arrays for SIMD
          Float* sol_row       = sol.at(i, 0);      // Current row of sol; we solve for this
          const Float* row_lo  = sol_row - sol_si;  // Left/previous row
          const Float* row_hi  = sol_row + sol_si;  // Right/next row
          const Float* rhs_row = rhs.at(i, 0);      // Current row of rhs

          // - Thomas algorithm ----------------------------
          // NOLINTBEGIN
          Float prev = 0.0;
          for (Index j = 0; j < ny; ++j) {
            // RHS for tridiagonal system; theta-derivative (x-derivative) is pulled to the right
            Float d = rhs_row[j * rhs_sj] - theta[j] * (row_lo[j * sol_sj] + row_hi[j * sol_sj]);
            // Periodic boundary conditions
            if (j == 0) { d -= bnd_lo * sol_row[-sol_sj]; }
            if (j == ny - 1) { d -= bnd_hi * sol_row[ny * sol_sj]; }

            prev                = (d - tri_a[j] * prev) * minv[j];
            sol_row[j * sol_sj] = prev;
          }

          for (Index j = ny - 2; j >= 0; --j) {
            sol_row[j * sol_sj] -= cstar[j] * sol_row[(j + 1) * sol_sj];
          }
          // NOLINTEND
          // - Thomas algorithm ----------------------------
        });
      }
    }
  }

  // -----------------------------------------------------------------------------------------------
  constexpr void smooth(const Level& level, Index num_iter) {
    // clang-format off
    switch (level.grid.coords()) {
      case Coordinates::CARTESIAN: return smooth_cartesian(level, num_iter);  // Red-black GS
      case Coordinates::POLAR:     return smooth_polar(level, num_iter);      // Zebra-line Thomas
    }
    // clang-format on
    Igor::Panic("Unreachable");
  }

  // -----------------------------------------------------------------------------------------------
  constexpr void restrict_residual([[maybe_unused]] const Level& level,
                                   const Scalar fine_res,  // can be level.res or level.rhs
                                   const Level& coarse) {
    IGOR_ASSERT(level.grid.nx() / 2 == coarse.grid.nx() && level.grid.ny() / 2 == coarse.grid.ny(),
                "Expected `coarse` to be the next coarser level but we skipped something.");
    const auto res = fine_res;
    const auto rhs = coarse.rhs;

    // Residual of `level` becomes the rhs of `coarse`. Volume weighted average.
    const Float* w_lo = coarse.restrict_w_lo.data();
    const Float* w_hi = coarse.restrict_w_hi.data();

    coarse.grid.foreach_i(FOREACH_FUNC {
      rhs(i, j) = w_lo[j] * (res(2 * i, 2 * j) + res(2 * i + 1, 2 * j)) +
                  w_hi[j] * (res(2 * i, 2 * j + 1) + res(2 * i + 1, 2 * j + 1));
    });
  }

  // -----------------------------------------------------------------------------------------------
  constexpr void prolongate_and_correct(const Level& coarse, const Level& level) {
    IGOR_ASSERT(level.grid.nx() / 2 == coarse.grid.nx() && level.grid.ny() / 2 == coarse.grid.ny(),
                "Expected `coarse` to be the next coarser level but we skipped something.");
    const auto lsol = level.sol;
    const auto csol = coarse.sol;

    apply_bconds(coarse.grid, m_bconds, coarse.sol, -1.0);

    // Bilinear interpolation of the coarse correction onto the finer solution
    constexpr Float W = 1.0 / 16.0;
    coarse.grid.foreach_i(FOREACH_FUNC {
      const Float center = 9.0 * csol(i, j);
      const Float left   = 3.0 * csol(i - 1, j);
      const Float right  = 3.0 * csol(i + 1, j);
      const Float bottom = 3.0 * csol(i, j - 1);
      const Float top    = 3.0 * csol(i, j + 1);

      // clang-format off
      lsol(2 * i, 2 * j)         += (center + left  + bottom + csol(i - 1, j - 1)) * W;
      lsol(2 * i, 2 * j + 1)     += (center + left  + top    + csol(i - 1, j + 1)) * W;
      lsol(2 * i + 1, 2 * j)     += (center + right + bottom + csol(i + 1, j - 1)) * W;
      lsol(2 * i + 1, 2 * j + 1) += (center + right + top    + csol(i + 1, j + 1)) * W;
      // clang-format on
    });
  }

  // -----------------------------------------------------------------------------------------------
  constexpr void vcycle(size_t l, bool sol_is_zero) {
    IGOR_ASSERT(l < m_levels.size(), "Level {} is out of bounds for {} levels", l, num_levels());
    Level& level = m_levels[l];

    // Check if we are at the coarsest level
    if (l + 1 == m_levels.size()) {
      smooth(level, m_num_iter_post);
      return;
    }

    Level& coarse = m_levels[l + 1];
    if (m_num_iter_pre > 0) {
      smooth(level, m_num_iter_pre);  // Remove high frequencies from residual
      residual(level);                // Iterate changed -> recompute residual
      sol_is_zero = false;
    }

    // Interpolate the residual of `level` onto `rhs` of coarse
    restrict_residual(level, sol_is_zero ? level.rhs : level.res, coarse);
    fill(coarse.sol, 0.0);
    vcycle(l + 1, true);  // Solve correction equation for `coarse`

    // Bilinear interpolation of the coarse correction onto level
    prolongate_and_correct(coarse, level);
    smooth(level, m_num_iter_post);
  }

 public:
  // -----------------------------------------------------------------------------------------------
  constexpr MultigridSolver(const Grid& grid,
                            BConds<Float> bconds = {.left   = Neumann{},
                                                    .right  = Neumann{},
                                                    .bottom = Neumann{},
                                                    .top    = Neumann{}},
                            Index min_size       = 2,
                            Index num_iter_pre   = 0,
                            Index num_iter_post  = 4)
      : m_bconds(std::move(bconds)),
        m_num_iter_pre(num_iter_pre),
        m_num_iter_post(num_iter_post) {
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
      precompute_tridiag(m_levels.back(), m_bconds);
      if (m_levels.size() > 1) {
        precompute_restriction_weights(m_levels[m_levels.size() - 2], m_levels.back());
      }

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

  // -----------------------------------------------------------------------------------------------
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

    residual(fine);
    m_res            = max_res(fine);
    Float res_before = m_res;
    for (m_num_cycles = 0; true; ++m_num_cycles) {
      if (m_res <= tol) {
        converged = true;
        break;
      }
      if (m_num_cycles == max_iter) { break; }

      vcycle(0, false);
      residual(fine);
      m_res = max_res(fine);

      // Dynamically adapt the number of smoothing iterations; adapted from Basilisk.
      // Pre- and post-smoothing move together.
      if (m_res > tol) {
        const Float reduction = res_before / m_res;
        if (reduction < 1.2) {
          if (m_num_iter_post < MAX_NUM_ITER_POST) { m_num_iter_post += 1; }
          if (m_num_iter_pre > 0 && m_num_iter_pre < MAX_NUM_ITER_PRE) { m_num_iter_pre += 1; }
        } else if (reduction > 10.0) {
          if (m_num_iter_post > MIN_NUM_ITER_POST) { m_num_iter_post -= 1; }
          if (m_num_iter_pre > MIN_NUM_ITER_PRE) { m_num_iter_pre -= 1; }
        }
      }
      res_before = m_res;
    }

    make_mean_free(fine.grid, fine.sol);
    copy(fine.sol, sol);

    return converged;
  }

  // -----------------------------------------------------------------------------------------------
  [[nodiscard]] constexpr auto num_levels() const noexcept -> Index {
    return static_cast<Index>(m_levels.size());
  }
  [[nodiscard]] constexpr auto num_cycles() const noexcept -> Index { return m_num_cycles; }
  [[nodiscard]] constexpr auto res() const noexcept -> Float { return m_res; }

  [[nodiscard]] constexpr auto num_iter_pre() noexcept -> Index& { return m_num_iter_pre; }
  [[nodiscard]] constexpr auto num_iter_pre() const noexcept -> Index { return m_num_iter_pre; }

  [[nodiscard]] constexpr auto num_iter_post() noexcept -> Index& { return m_num_iter_post; }
  [[nodiscard]] constexpr auto num_iter_post() const noexcept -> Index { return m_num_iter_post; }
};
