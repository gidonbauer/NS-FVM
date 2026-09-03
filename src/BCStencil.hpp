#pragma once

#include <algorithm>
#include <array>
#include <cmath>

#include <Igor/Logging.hpp>

#include "Grid.hpp"

// Highest polynomial degree supported by the boundary extrapolation.
inline constexpr Index MAX_BC_ORDER = 5;
// Highest number of ghost layers a single boundary stencil can fill.
inline constexpr Index MAX_BC_NGHOST = 8;

// Position of the sample points of a field relative to the boundary, along the boundary normal.
//   ALIGNED   : a sample point lies exactly on the boundary (u.x on the left/right walls, u.y on
//               the bottom/top walls).
//   STAGGERED : the boundary falls halfway between the last interior and the first ghost sample
//               (u.y on the left/right walls, u.x on the bottom/top walls, cell-centered scalars).
enum class Stagger { ALIGNED, STAGGERED };

enum class BCKind { DIRICHLET, NEUMANN };

// -------------------------------------------------------------------------------------------------
template <typename Float>
struct BCStencil {
  // Ghost layer m (m == 0 being the layer closest to the wall) is reconstructed as
  //     f_ghost(m) = w[m][0] * datum + sum_k w[m][1 + k] * f_sample(k)
  // with `f_sample(k)` the k-th interior sample away from the wall and `datum` the boundary datum
  // scaled as documented in `make_bc_stencil`.
  std::array<std::array<Float, MAX_BC_ORDER + 1>, MAX_BC_NGHOST> w{};
  Index nghost = 0;
  Index order  = 0;
};

// -------------------------------------------------------------------------------------------------
// Build the ghost-fill stencil for one boundary.
//
// Everything happens in the scaled coordinate xi = s / h, where s is the distance from the wall,
// h the mesh spacing along the normal and xi grows towards the interior. The interior samples sit
// at xi_k = k + 1/2 (STAGGERED) or xi_k = k + 1 (ALIGNED, the sample at xi = 0 is the wall itself)
// and the ghosts at xi = -(m + 1/2) resp. xi = -(m + 1).
//
// The ghost values are those of the unique polynomial p of degree `order` with
//     p(0)     = datum                      (DIRICHLET)
//     p'(0)    = datum                      (NEUMANN)
//     p(xi_k)  = f_sample(k)                for k = 0 ... order - 1
//
// Because xi points into the domain, the NEUMANN datum is the *inward* normal derivative scaled by
// h, i.e. the caller has to pass `-h * df/dn` with n the outward normal.
//
// `order == 1` reproduces the familiar two-point formulas: `2 * g - f_sample(0)` for a STAGGERED
// Dirichlet boundary, `2 * g - f_sample(0)` (odd reflection) for an ALIGNED one, and the mirror for
// the *first* Neumann ghost layer. Note that for Neumann the deeper layers are *not* mirrored -- a
// degree-1 polynomial with a prescribed derivative at the wall extends the field linearly, which
// beyond the first layer is worse than reflection. Use `mirror_ghosts_x` / `mirror_ghosts_y` for
// that; see the comment there.
template <typename Float>
[[nodiscard]] constexpr auto
make_bc_stencil(BCKind kind, Stagger stagger, Index nghost, Index order) noexcept
    -> BCStencil<Float> {
  IGOR_ASSERT(order >= 1 && order <= MAX_BC_ORDER,
              "Boundary extrapolation order must be in [1, {}], but is {}",
              MAX_BC_ORDER,
              order);
  IGOR_ASSERT(nghost >= 1 && nghost <= MAX_BC_NGHOST,
              "Number of ghost cells must be in [1, {}], but is {}",
              MAX_BC_NGHOST,
              nghost);

  constexpr Index MAXN = MAX_BC_ORDER + 1;
  const Index n        = order + 1;

  const auto idx = [](Index r, Index c) -> size_t {
    return static_cast<size_t>(r) * static_cast<size_t>(MAXN) + static_cast<size_t>(c);
  };
  const auto xi_sample = [stagger](Index k) -> Float {
    return stagger == Stagger::STAGGERED ? static_cast<Float>(k) + static_cast<Float>(0.5)
                                         : static_cast<Float>(k + 1);
  };
  const auto xi_ghost = [stagger](Index m) -> Float {
    return stagger == Stagger::STAGGERED ? -(static_cast<Float>(m) + static_cast<Float>(0.5))
                                         : -static_cast<Float>(m + 1);
  };

  // A * a = b for the monomial coefficients a of p; b = (datum, f_sample(0), ..., f_sample(d-1)).
  std::array<Float, MAXN * MAXN> A{};
  std::array<Float, MAXN * MAXN> Ainv{};

  // Row 0: the boundary condition.
  A[idx(0, kind == BCKind::DIRICHLET ? 0 : 1)] = 1.0;
  // Rows 1 ... order: interpolation of the interior samples.
  for (Index k = 0; k < order; ++k) {
    Float p = 1.0;
    for (Index c = 0; c < n; ++c) {
      A[idx(k + 1, c)]  = p;
      p                *= xi_sample(k);
    }
  }

  // Gauss-Jordan inversion with partial pivoting; n <= MAX_BC_ORDER + 1, so this is cheap.
  for (Index i = 0; i < n; ++i) {
    Ainv[idx(i, i)] = 1.0;
  }
  for (Index c = 0; c < n; ++c) {
    Index piv = c;
    for (Index r = c + 1; r < n; ++r) {
      if (std::abs(A[idx(r, c)]) > std::abs(A[idx(piv, c)])) { piv = r; }
    }
    IGOR_ASSERT(std::abs(A[idx(piv, c)]) > 0.0,
                "Singular boundary stencil system for order {}.",
                order);
    if (piv != c) {
      for (Index k = 0; k < n; ++k) {
        std::swap(A[idx(c, k)], A[idx(piv, k)]);
        std::swap(Ainv[idx(c, k)], Ainv[idx(piv, k)]);
      }
    }
    const Float inv_pivot = 1.0 / A[idx(c, c)];
    for (Index k = 0; k < n; ++k) {
      A[idx(c, k)]    *= inv_pivot;
      Ainv[idx(c, k)] *= inv_pivot;
    }
    for (Index r = 0; r < n; ++r) {
      if (r == c || A[idx(r, c)] == 0.0) { continue; }
      const Float f = A[idx(r, c)];
      for (Index k = 0; k < n; ++k) {
        A[idx(r, k)]    -= f * A[idx(c, k)];
        Ainv[idx(r, k)] -= f * Ainv[idx(c, k)];
      }
    }
  }

  // p(xi_g) = sum_c xi_g^c * a_c = sum_c xi_g^c * sum_j Ainv(c, j) * b_j
  BCStencil<Float> stencil{.w = {}, .nghost = nghost, .order = order};
  for (Index m = 0; m < nghost; ++m) {
    Float p = 1.0;
    for (Index c = 0; c < n; ++c) {
      for (Index j = 0; j < n; ++j) {
        stencil.w[static_cast<size_t>(m)][static_cast<size_t>(j)] += p * Ainv[idx(c, j)];
      }
      p *= xi_ghost(m);
    }
  }
  return stencil;
}

// -------------------------------------------------------------------------------------------------
// Index helpers along the boundary normal. `n` is the number of interior samples of the field along
// that direction (`f.nx()` or `f.ny()`), `high_side` selects the right/top instead of the
// left/bottom wall.

// Index of the k-th interior sample away from the wall; k == 0 is the closest one.
[[nodiscard]] constexpr auto
bc_sample_idx(Index n, Stagger stagger, bool high_side, Index k) noexcept -> Index {
  const Index off = stagger == Stagger::ALIGNED ? 1 : 0;
  return high_side ? n - 1 - off - k : off + k;
}

// Index of the m-th ghost layer; m == 0 is the one closest to the wall.
[[nodiscard]] constexpr auto bc_ghost_idx(Index n, bool high_side, Index m) noexcept -> Index {
  return high_side ? n + m : -1 - m;
}

// Index of the sample sitting exactly on the wall; only meaningful for Stagger::ALIGNED.
[[nodiscard]] constexpr auto bc_wall_idx(Index n, bool high_side) noexcept -> Index {
  return high_side ? n - 1 : 0;
}

// Index the m-th ghost layer wraps onto for a periodic boundary. For ALIGNED fields the samples on
// the two walls are the same physical location, hence the period is one sample shorter.
[[nodiscard]] constexpr auto
bc_periodic_idx(Index n, Stagger stagger, bool high_side, Index m) noexcept -> Index {
  const Index period = stagger == Stagger::ALIGNED ? n - 1 : n;
  const Index ghost  = bc_ghost_idx(n, high_side, m);
  return high_side ? ghost - period : ghost + period;
}

// -------------------------------------------------------------------------------------------------
// Assert that the stencil only reads interior samples.
template <typename Float, Layout LAYOUT>
constexpr void
check_bc_stencil_fits(const Scalar<Float, LAYOUT>& f, Stagger stagger, Index n, Index order) {
  const Index off = stagger == Stagger::ALIGNED ? 1 : 0;
  IGOR_ASSERT(off + order <= n,
              "Boundary extrapolation of order {} needs {} interior samples but only {} are "
              "available along this direction.",
              order,
              off + order,
              n);
  IGOR_ASSERT(f.nghost() <= MAX_BC_NGHOST,
              "Number of ghost cells must be at most {}, but is {}",
              MAX_BC_NGHOST,
              f.nghost());
}

// -------------------------------------------------------------------------------------------------
// Fill the ghost layers of `f` on a boundary whose normal is the x-axis, for the rows
// [jlo, jhi). `datum(j)` returns the boundary datum (see `make_bc_stencil`) and `post(v)`
// post-processes the reconstructed ghost value (used by the clipped outflow condition).
template <typename Float, Layout LAYOUT, typename DATUM, typename POST>
constexpr void fill_ghosts_x(const Grid<Float, LAYOUT>& grid,
                             Scalar<Float, LAYOUT> f,
                             Stagger stagger,
                             bool high_side,
                             const BCStencil<Float>& stencil,
                             DATUM datum,
                             POST post,
                             Index jlo,
                             Index jhi) noexcept {
  check_bc_stencil_fits(f, stagger, f.nx(), stencil.order);
  grid.foreach_range(0, 1, jlo, jhi, [=](Index, Index j) {
    const Float b = datum(j);
    for (Index m = 0; m < stencil.nghost; ++m) {
      Float v = stencil.w[static_cast<size_t>(m)][0] * b;
      for (Index k = 0; k < stencil.order; ++k) {
        v += stencil.w[static_cast<size_t>(m)][static_cast<size_t>(k + 1)] *
             f(bc_sample_idx(f.nx(), stagger, high_side, k), j);
      }
      f(bc_ghost_idx(f.nx(), high_side, m), j) = post(v);
    }
  });
}

// Same as `fill_ghosts_x`, but for a boundary whose normal is the y-axis.
template <typename Float, Layout LAYOUT, typename DATUM, typename POST>
constexpr void fill_ghosts_y(const Grid<Float, LAYOUT>& grid,
                             Scalar<Float, LAYOUT> f,
                             Stagger stagger,
                             bool high_side,
                             const BCStencil<Float>& stencil,
                             DATUM datum,
                             POST post,
                             Index ilo,
                             Index ihi) noexcept {
  check_bc_stencil_fits(f, stagger, f.ny(), stencil.order);
  grid.foreach_range(ilo, ihi, 0, 1, [=](Index i, Index) {
    const Float b = datum(i);
    for (Index m = 0; m < stencil.nghost; ++m) {
      Float v = stencil.w[static_cast<size_t>(m)][0] * b;
      for (Index k = 0; k < stencil.order; ++k) {
        v += stencil.w[static_cast<size_t>(m)][static_cast<size_t>(k + 1)] *
             f(i, bc_sample_idx(f.ny(), stagger, high_side, k));
      }
      f(i, bc_ghost_idx(f.ny(), high_side, m)) = post(v);
    }
  });
}

// -------------------------------------------------------------------------------------------------
// Distance of the k-th interior sample from the wall, in units of the mesh spacing.
[[nodiscard]] constexpr auto bc_sample_dist(Stagger stagger, Index k) noexcept -> double {
  return stagger == Stagger::STAGGERED ? static_cast<double>(k) + 0.5 : static_cast<double>(k + 1);
}

// Even reflection of the interior samples about the wall, plus the linear ramp that carries a
// non-zero outward normal derivative `dfdn`:
//     f_ghost(m) = f_sample(m) + 2 * h * dist(m) * dfdn
// Ghost layer m and interior sample m are mirror partners, so this is a pure copy for the
// homogeneous case: `f(-1) = f(0)` for a STAGGERED and `f(-1) = f(1)` for an ALIGNED field.
//
// Reflection is not a member of the `make_bc_stencil` family: it is exact whenever the odd part of
// the field about the wall is linear, hence third order accurate for a general smooth field, while
// its weights stay bounded by one. That makes it both more accurate and far more robust than the
// order-1 polynomial fit, whose zero-derivative constraint degenerates into a constant extension
// beyond the first ghost layer.
template <typename Float, Layout LAYOUT, typename POST>
constexpr void mirror_ghosts_x(const Grid<Float, LAYOUT>& grid,
                               Scalar<Float, LAYOUT> f,
                               Stagger stagger,
                               bool high_side,
                               Float dfdn,
                               POST post,
                               Index jlo,
                               Index jhi) noexcept {
  check_bc_stencil_fits(f, stagger, f.nx(), f.nghost());
  const Float h = grid.dx();
  grid.foreach_range(0, 1, jlo, jhi, [=](Index, Index j) {
    for (Index m = 0; m < f.nghost(); ++m) {
      const Float ramp = 2.0 * h * static_cast<Float>(bc_sample_dist(stagger, m)) * dfdn;
      f(bc_ghost_idx(f.nx(), high_side, m), j) =
          post(f(bc_sample_idx(f.nx(), stagger, high_side, m), j) + ramp);
    }
  });
}

// Same as `mirror_ghosts_x`, but for a boundary whose normal is the y-axis.
template <typename Float, Layout LAYOUT, typename POST>
constexpr void mirror_ghosts_y(const Grid<Float, LAYOUT>& grid,
                               Scalar<Float, LAYOUT> f,
                               Stagger stagger,
                               bool high_side,
                               Float dfdn,
                               POST post,
                               Index ilo,
                               Index ihi) noexcept {
  check_bc_stencil_fits(f, stagger, f.ny(), f.nghost());
  const Float h = grid.dy();
  grid.foreach_range(ilo, ihi, 0, 1, [=](Index i, Index) {
    for (Index m = 0; m < f.nghost(); ++m) {
      const Float ramp = 2.0 * h * static_cast<Float>(bc_sample_dist(stagger, m)) * dfdn;
      f(i, bc_ghost_idx(f.ny(), high_side, m)) =
          post(f(i, bc_sample_idx(f.ny(), stagger, high_side, m)) + ramp);
    }
  });
}

// -------------------------------------------------------------------------------------------------
// Post-processing functors for `fill_ghosts_x` / `fill_ghosts_y`.
struct BCNoPost {
  template <typename Float>
  constexpr auto operator()(Float v) const noexcept -> Float {
    return v;
  }
};

// Suppress inflow through a low-side (left/bottom) outflow boundary.
struct BCClipInflowLow {
  template <typename Float>
  constexpr auto operator()(Float v) const noexcept -> Float {
    return std::min(v, static_cast<Float>(0));
  }
};

// Suppress inflow through a high-side (right/top) outflow boundary.
struct BCClipInflowHigh {
  template <typename Float>
  constexpr auto operator()(Float v) const noexcept -> Float {
    return std::max(v, static_cast<Float>(0));
  }
};

inline constexpr BCNoPost bc_no_post{};
inline constexpr BCClipInflowLow bc_clip_inflow_low{};
inline constexpr BCClipInflowHigh bc_clip_inflow_high{};
