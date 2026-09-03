#pragma once

#include <variant>

#include "BCStencil.hpp"
#include "Grid.hpp"

// The four boundary conditions below fill an arbitrary number of ghost layers. All of them share
// the same convention:
//
//   * `apply_left` / `apply_right` only touch the rows a component actually has interior data in,
//     whereas `apply_bottom` / `apply_top` sweep the *extended* i-range including the ghost
//     columns. `apply_velocity_bconds` runs the x-boundaries first, so the ghost corners end up as
//     the tensor product of the two one-dimensional fills and nothing ever reads uninitialized
//     data.
//
//   * A component that has a sample point exactly on the wall (u.x on left/right, u.y on
//     bottom/top) is Stagger::ALIGNED there, the other one is Stagger::STAGGERED. See
//     `make_bc_stencil` for how the ghost values are reconstructed.
//
//   * `Dirichlet::order` is the degree of the polynomial used for the ghost reconstruction.
//     `order == 1` is linear extrapolation, i.e. `2 * g - u_0`. Raise it to match a higher order
//     interior discretization (e.g. `order = 3` alongside WENO-5), keeping in mind that the
//     extrapolation weights grow quickly with the order and that the ghost *corners*, being filled
//     by two successive extrapolations, amplify that growth.
//
//   * `Neumann` defaults to an even reflection about the wall rather than to a polynomial fit; see
//     `NeumannFill` below.

// -------------------------------------------------------------------------------------------------
template <typename Float>
struct Dirichlet {
  std::variant<Float, Float (*)(Float, Float)> U{}, V{};
  Index order = 1;

  // = LEFT ========
  template <Layout LAYOUT>
  constexpr void apply_left(const Grid<Float, LAYOUT>& grid,
                            FaceVector<Float, LAYOUT> uf,
                            Float t) const noexcept {
    IGOR_ASSERT(grid.nghost() >= 1, "Expected at least one ghost cell, but got {}", grid.nghost());

    const auto U_bc = [=, this](Index j) -> Float {
      return std::holds_alternative<Float>(U) ? std::get<0>(U) : std::get<1>(U)(grid.ym(j), t);
    };
    const auto V_bc = [=, this](Index j) -> Float {
      return std::holds_alternative<Float>(V) ? std::get<0>(V) : std::get<1>(V)(grid.y(j), t);
    };

    // u.x lies on the wall: prescribe it there, then extrapolate into the ghosts.
    grid.foreach_range(0, 1, 0, uf.x.ny(), [=](Index, Index j) {
      uf.x(bc_wall_idx(uf.x.nx(), false), j) = U_bc(j);
    });
    const auto sx =
        make_bc_stencil<Float>(BCKind::DIRICHLET, Stagger::ALIGNED, uf.x.nghost(), order);
    fill_ghosts_x(grid, uf.x, Stagger::ALIGNED, false, sx, U_bc, bc_no_post, 0, uf.x.ny());

    // u.y straddles the wall: the value is only enforced through the reconstruction.
    const auto sy =
        make_bc_stencil<Float>(BCKind::DIRICHLET, Stagger::STAGGERED, uf.y.nghost(), order);
    fill_ghosts_x(grid, uf.y, Stagger::STAGGERED, false, sy, V_bc, bc_no_post, 0, uf.y.ny());
  }

  // = RIGHT =======
  template <Layout LAYOUT>
  constexpr void apply_right(const Grid<Float, LAYOUT>& grid,
                             FaceVector<Float, LAYOUT> uf,
                             Float t) const noexcept {
    IGOR_ASSERT(grid.nghost() >= 1, "Expected at least one ghost cell, but got {}", grid.nghost());

    const auto U_bc = [=, this](Index j) -> Float {
      return std::holds_alternative<Float>(U) ? std::get<0>(U) : std::get<1>(U)(grid.ym(j), t);
    };
    const auto V_bc = [=, this](Index j) -> Float {
      return std::holds_alternative<Float>(V) ? std::get<0>(V) : std::get<1>(V)(grid.y(j), t);
    };

    grid.foreach_range(0, 1, 0, uf.x.ny(), [=](Index, Index j) {
      uf.x(bc_wall_idx(uf.x.nx(), true), j) = U_bc(j);
    });
    const auto sx =
        make_bc_stencil<Float>(BCKind::DIRICHLET, Stagger::ALIGNED, uf.x.nghost(), order);
    fill_ghosts_x(grid, uf.x, Stagger::ALIGNED, true, sx, U_bc, bc_no_post, 0, uf.x.ny());

    const auto sy =
        make_bc_stencil<Float>(BCKind::DIRICHLET, Stagger::STAGGERED, uf.y.nghost(), order);
    fill_ghosts_x(grid, uf.y, Stagger::STAGGERED, true, sy, V_bc, bc_no_post, 0, uf.y.ny());
  }

  // = BOTTOM ======
  template <Layout LAYOUT>
  constexpr void apply_bottom(const Grid<Float, LAYOUT>& grid,
                              FaceVector<Float, LAYOUT> uf,
                              Float t) const noexcept {
    IGOR_ASSERT(grid.nghost() >= 1, "Expected at least one ghost cell, but got {}", grid.nghost());

    const auto U_bc = [=, this](Index i) -> Float {
      return std::holds_alternative<Float>(U) ? std::get<0>(U) : std::get<1>(U)(grid.x(i), t);
    };
    const auto V_bc = [=, this](Index i) -> Float {
      return std::holds_alternative<Float>(V) ? std::get<0>(V) : std::get<1>(V)(grid.xm(i), t);
    };

    // u.x straddles the wall.
    const auto sx =
        make_bc_stencil<Float>(BCKind::DIRICHLET, Stagger::STAGGERED, uf.x.nghost(), order);
    fill_ghosts_y(grid,
                  uf.x,
                  Stagger::STAGGERED,
                  false,
                  sx,
                  U_bc,
                  bc_no_post,
                  -uf.x.nghost(),
                  uf.x.nx() + uf.x.nghost());

    // u.y lies on the wall.
    grid.foreach_range(-uf.y.nghost(), uf.y.nx() + uf.y.nghost(), 0, 1, [=](Index i, Index) {
      uf.y(i, bc_wall_idx(uf.y.ny(), false)) = V_bc(i);
    });
    const auto sy =
        make_bc_stencil<Float>(BCKind::DIRICHLET, Stagger::ALIGNED, uf.y.nghost(), order);
    fill_ghosts_y(grid,
                  uf.y,
                  Stagger::ALIGNED,
                  false,
                  sy,
                  V_bc,
                  bc_no_post,
                  -uf.y.nghost(),
                  uf.y.nx() + uf.y.nghost());
  }

  // = TOP =========
  template <Layout LAYOUT>
  constexpr void
  apply_top(const Grid<Float, LAYOUT>& grid, FaceVector<Float, LAYOUT> uf, Float t) const noexcept {
    IGOR_ASSERT(grid.nghost() >= 1, "Expected at least one ghost cell, but got {}", grid.nghost());

    const auto U_bc = [=, this](Index i) -> Float {
      return std::holds_alternative<Float>(U) ? std::get<0>(U) : std::get<1>(U)(grid.x(i), t);
    };
    const auto V_bc = [=, this](Index i) -> Float {
      return std::holds_alternative<Float>(V) ? std::get<0>(V) : std::get<1>(V)(grid.xm(i), t);
    };

    const auto sx =
        make_bc_stencil<Float>(BCKind::DIRICHLET, Stagger::STAGGERED, uf.x.nghost(), order);
    fill_ghosts_y(grid,
                  uf.x,
                  Stagger::STAGGERED,
                  true,
                  sx,
                  U_bc,
                  bc_no_post,
                  -uf.x.nghost(),
                  uf.x.nx() + uf.x.nghost());

    grid.foreach_range(-uf.y.nghost(), uf.y.nx() + uf.y.nghost(), 0, 1, [=](Index i, Index) {
      uf.y(i, bc_wall_idx(uf.y.ny(), true)) = V_bc(i);
    });
    const auto sy =
        make_bc_stencil<Float>(BCKind::DIRICHLET, Stagger::ALIGNED, uf.y.nghost(), order);
    fill_ghosts_y(grid,
                  uf.y,
                  Stagger::ALIGNED,
                  true,
                  sy,
                  V_bc,
                  bc_no_post,
                  -uf.y.nghost(),
                  uf.y.nx() + uf.y.nghost());
  }
};

// -------------------------------------------------------------------------------------------------
// How the ghosts of a Neumann boundary are built.
//   MIRROR      : even reflection about the wall plus the linear ramp carrying a non-zero normal
//                 derivative. Third order accurate with weights bounded by one -- the robust
//                 default. For the component lying on the wall this means `u.x(-1, j) = u.x(1, j)`,
//                 *not* `u.x(-1, j) = u.x(0, j)`, which would centre the zero-gradient condition
//                 half a cell off the boundary and is only first order accurate.
//   EXTRAPOLATE : the constrained polynomial fit of degree `order`. Needed to go beyond third
//                 order, at the price of rapidly growing extrapolation weights.
enum class NeumannFill { MIRROR, EXTRAPOLATE };

// -------------------------------------------------------------------------------------------------
// Prescribes the *outward* normal derivative of both velocity components.
struct Neumann {
  // Clip the normal velocity in the ghosts so that the boundary cannot draw fluid back in.
  bool clipped     = false;
  double dUdn      = 0.0;
  double dVdn      = 0.0;
  NeumannFill fill = NeumannFill::MIRROR;
  // Only used with NeumannFill::EXTRAPOLATE.
  Index order = 1;

  // = LEFT ========
  template <typename Float, Layout LAYOUT>
  constexpr void apply_left(const Grid<Float, LAYOUT>& grid,
                            FaceVector<Float, LAYOUT> uf,
                            Float /*t*/) const noexcept {
    IGOR_ASSERT(grid.nghost() >= 1, "Expected at least one ghost cell, but got {}", grid.nghost());

    apply_x(grid, uf, false, 0, uf.x.ny(), 0, uf.y.ny());
  }

  // = RIGHT =======
  template <typename Float, Layout LAYOUT>
  constexpr void apply_right(const Grid<Float, LAYOUT>& grid,
                             FaceVector<Float, LAYOUT> uf,
                             Float /*t*/) const noexcept {
    IGOR_ASSERT(grid.nghost() >= 1, "Expected at least one ghost cell, but got {}", grid.nghost());

    apply_x(grid, uf, true, 0, uf.x.ny(), 0, uf.y.ny());
  }

  // = BOTTOM ======
  template <typename Float, Layout LAYOUT>
  constexpr void apply_bottom(const Grid<Float, LAYOUT>& grid,
                              FaceVector<Float, LAYOUT> uf,
                              Float /*t*/) const noexcept {
    IGOR_ASSERT(grid.nghost() >= 1, "Expected at least one ghost cell, but got {}", grid.nghost());

    apply_y(grid,
            uf,
            false,
            -uf.x.nghost(),
            uf.x.nx() + uf.x.nghost(),
            -uf.y.nghost(),
            uf.y.nx() + uf.y.nghost());
  }

  // = TOP =========
  template <typename Float, Layout LAYOUT>
  constexpr void apply_top(const Grid<Float, LAYOUT>& grid,
                           FaceVector<Float, LAYOUT> uf,
                           Float /*t*/) const noexcept {
    IGOR_ASSERT(grid.nghost() >= 1, "Expected at least one ghost cell, but got {}", grid.nghost());

    apply_y(grid,
            uf,
            true,
            -uf.x.nghost(),
            uf.x.nx() + uf.x.nghost(),
            -uf.y.nghost(),
            uf.y.nx() + uf.y.nghost());
  }

 private:
  // On the left/right walls u.x is ALIGNED and u.y is STAGGERED; the clipping acts on the normal
  // component, i.e. on u.x.
  template <typename Float, Layout LAYOUT>
  constexpr void apply_x(const Grid<Float, LAYOUT>& grid,
                         FaceVector<Float, LAYOUT> uf,
                         bool high_side,
                         Index jlo_x,
                         Index jhi_x,
                         Index jlo_y,
                         Index jhi_y) const noexcept {
    const auto du = static_cast<Float>(dUdn);
    const auto dv = static_cast<Float>(dVdn);
    if (clipped && high_side) {
      fill_x(grid, uf.x, Stagger::ALIGNED, high_side, du, bc_clip_inflow_high, jlo_x, jhi_x);
    } else if (clipped) {
      fill_x(grid, uf.x, Stagger::ALIGNED, high_side, du, bc_clip_inflow_low, jlo_x, jhi_x);
    } else {
      fill_x(grid, uf.x, Stagger::ALIGNED, high_side, du, bc_no_post, jlo_x, jhi_x);
    }
    fill_x(grid, uf.y, Stagger::STAGGERED, high_side, dv, bc_no_post, jlo_y, jhi_y);
  }

  // On the bottom/top walls the roles are swapped: u.y is ALIGNED and carries the clipping.
  template <typename Float, Layout LAYOUT>
  constexpr void apply_y(const Grid<Float, LAYOUT>& grid,
                         FaceVector<Float, LAYOUT> uf,
                         bool high_side,
                         Index ilo_x,
                         Index ihi_x,
                         Index ilo_y,
                         Index ihi_y) const noexcept {
    const auto du = static_cast<Float>(dUdn);
    const auto dv = static_cast<Float>(dVdn);
    fill_y(grid, uf.x, Stagger::STAGGERED, high_side, du, bc_no_post, ilo_x, ihi_x);
    if (clipped && high_side) {
      fill_y(grid, uf.y, Stagger::ALIGNED, high_side, dv, bc_clip_inflow_high, ilo_y, ihi_y);
    } else if (clipped) {
      fill_y(grid, uf.y, Stagger::ALIGNED, high_side, dv, bc_clip_inflow_low, ilo_y, ihi_y);
    } else {
      fill_y(grid, uf.y, Stagger::ALIGNED, high_side, dv, bc_no_post, ilo_y, ihi_y);
    }
  }

  template <typename Float, Layout LAYOUT, typename POST>
  constexpr void fill_x(const Grid<Float, LAYOUT>& grid,
                        Scalar<Float, LAYOUT> f,
                        Stagger stagger,
                        bool high_side,
                        Float dfdn,
                        POST post,
                        Index jlo,
                        Index jhi) const noexcept {
    if (fill == NeumannFill::MIRROR) {
      mirror_ghosts_x(grid, f, stagger, high_side, dfdn, post, jlo, jhi);
    } else {
      // `make_bc_stencil` wants the inward normal derivative scaled by the mesh spacing.
      const Float b = -grid.dx() * dfdn;
      const auto s  = make_bc_stencil<Float>(BCKind::NEUMANN, stagger, f.nghost(), order);
      fill_ghosts_x(
          grid, f, stagger, high_side, s, [=](Index) -> Float { return b; }, post, jlo, jhi);
    }
  }

  template <typename Float, Layout LAYOUT, typename POST>
  constexpr void fill_y(const Grid<Float, LAYOUT>& grid,
                        Scalar<Float, LAYOUT> f,
                        Stagger stagger,
                        bool high_side,
                        Float dfdn,
                        POST post,
                        Index ilo,
                        Index ihi) const noexcept {
    if (fill == NeumannFill::MIRROR) {
      mirror_ghosts_y(grid, f, stagger, high_side, dfdn, post, ilo, ihi);
    } else {
      const Float b = -grid.dy() * dfdn;
      const auto s  = make_bc_stencil<Float>(BCKind::NEUMANN, stagger, f.nghost(), order);
      fill_ghosts_y(
          grid, f, stagger, high_side, s, [=](Index) -> Float { return b; }, post, ilo, ihi);
    }
  }
};

// -------------------------------------------------------------------------------------------------
struct Periodic {
  // = LEFT ========
  template <typename Float, Layout LAYOUT>
  constexpr void apply_left(const Grid<Float, LAYOUT>& grid,
                            FaceVector<Float, LAYOUT> uf,
                            Float /*t*/) const noexcept {
    IGOR_ASSERT(grid.nghost() >= 1, "Expected at least one ghost cell, but got {}", grid.nghost());
    grid.foreach_range(0, 1, 0, uf.x.ny(), [=](Index, Index j) {
      for (Index m = 0; m < uf.x.nghost(); ++m) {
        uf.x(bc_ghost_idx(uf.x.nx(), false, m), j) =
            uf.x(bc_periodic_idx(uf.x.nx(), Stagger::ALIGNED, false, m), j);
      }
    });
    grid.foreach_range(0, 1, 0, uf.y.ny(), [=](Index, Index j) {
      for (Index m = 0; m < uf.y.nghost(); ++m) {
        uf.y(bc_ghost_idx(uf.y.nx(), false, m), j) =
            uf.y(bc_periodic_idx(uf.y.nx(), Stagger::STAGGERED, false, m), j);
      }
    });
  }

  // = RIGHT =======
  template <typename Float, Layout LAYOUT>
  constexpr void apply_right(const Grid<Float, LAYOUT>& grid,
                             FaceVector<Float, LAYOUT> uf,
                             Float /*t*/) const noexcept {
    IGOR_ASSERT(grid.nghost() >= 1, "Expected at least one ghost cell, but got {}", grid.nghost());
    grid.foreach_range(0, 1, 0, uf.x.ny(), [=](Index, Index j) {
      for (Index m = 0; m < uf.x.nghost(); ++m) {
        uf.x(bc_ghost_idx(uf.x.nx(), true, m), j) =
            uf.x(bc_periodic_idx(uf.x.nx(), Stagger::ALIGNED, true, m), j);
      }
    });
    grid.foreach_range(0, 1, 0, uf.y.ny(), [=](Index, Index j) {
      for (Index m = 0; m < uf.y.nghost(); ++m) {
        uf.y(bc_ghost_idx(uf.y.nx(), true, m), j) =
            uf.y(bc_periodic_idx(uf.y.nx(), Stagger::STAGGERED, true, m), j);
      }
    });
  }

  // = BOTTOM ======
  template <typename Float, Layout LAYOUT>
  constexpr void apply_bottom(const Grid<Float, LAYOUT>& grid,
                              FaceVector<Float, LAYOUT> uf,
                              Float /*t*/) const noexcept {
    IGOR_ASSERT(grid.nghost() >= 1, "Expected at least one ghost cell, but got {}", grid.nghost());
    grid.foreach_range(
        -uf.x.nghost(), uf.x.nx() + uf.x.nghost(), 0, 1, [=](Index i, Index) {
          for (Index m = 0; m < uf.x.nghost(); ++m) {
            uf.x(i, bc_ghost_idx(uf.x.ny(), false, m)) =
                uf.x(i, bc_periodic_idx(uf.x.ny(), Stagger::STAGGERED, false, m));
          }
        });
    grid.foreach_range(
        -uf.y.nghost(), uf.y.nx() + uf.y.nghost(), 0, 1, [=](Index i, Index) {
          for (Index m = 0; m < uf.y.nghost(); ++m) {
            uf.y(i, bc_ghost_idx(uf.y.ny(), false, m)) =
                uf.y(i, bc_periodic_idx(uf.y.ny(), Stagger::ALIGNED, false, m));
          }
        });
  }

  // = TOP =========
  template <typename Float, Layout LAYOUT>
  constexpr void apply_top(const Grid<Float, LAYOUT>& grid,
                           FaceVector<Float, LAYOUT> uf,
                           Float /*t*/) const noexcept {
    IGOR_ASSERT(grid.nghost() >= 1, "Expected at least one ghost cell, but got {}", grid.nghost());
    grid.foreach_range(
        -uf.x.nghost(), uf.x.nx() + uf.x.nghost(), 0, 1, [=](Index i, Index) {
          for (Index m = 0; m < uf.x.nghost(); ++m) {
            uf.x(i, bc_ghost_idx(uf.x.ny(), true, m)) =
                uf.x(i, bc_periodic_idx(uf.x.ny(), Stagger::STAGGERED, true, m));
          }
        });
    grid.foreach_range(
        -uf.y.nghost(), uf.y.nx() + uf.y.nghost(), 0, 1, [=](Index i, Index) {
          for (Index m = 0; m < uf.y.nghost(); ++m) {
            uf.y(i, bc_ghost_idx(uf.y.ny(), true, m)) =
                uf.y(i, bc_periodic_idx(uf.y.ny(), Stagger::ALIGNED, true, m));
          }
        });
  }
};

// -------------------------------------------------------------------------------------------------
template <typename Float>
using BCond_t = std::variant<Dirichlet<Float>, Neumann, Periodic>;

template <typename Float>
struct VelocityBConds {
  BCond_t<Float> left;
  BCond_t<Float> right;
  BCond_t<Float> bottom;
  BCond_t<Float> top;
};

// -------------------------------------------------------------------------------------------------
template <typename Float, Layout LAYOUT>
void apply_velocity_bconds(const Grid<Float, LAYOUT>& grid,
                           const VelocityBConds<Float>& bconds,
                           FaceVector<Float, LAYOUT> uf,
                           Float t = -1.0) {
  // The order matters: left/right fill the interior rows only, bottom/top then sweep the extended
  // i-range and thereby fill the ghost corners as a tensor product of the two 1D reconstructions.
  std::visit([&](auto&& bcond) { bcond.apply_left(grid, uf, t); }, bconds.left);
  std::visit([&](auto&& bcond) { bcond.apply_right(grid, uf, t); }, bconds.right);
  std::visit([&](auto&& bcond) { bcond.apply_bottom(grid, uf, t); }, bconds.bottom);
  std::visit([&](auto&& bcond) { bcond.apply_top(grid, uf, t); }, bconds.top);
}

// -------------------------------------------------------------------------------------------------
// Homogeneous Neumann on all four sides of a cell-centered scalar, i.e. a mirror of the interior
// values about the wall.
template <typename Float, Layout LAYOUT>
constexpr void apply_neumann_bconds(const Grid<Float, LAYOUT>& grid,
                                    Scalar<Float, LAYOUT> field) noexcept {
  IGOR_ASSERT(field.nghost() > 0, "Expected at least one ghost cell, but got {}", field.nghost());

  const Float zero = 0.0;
  mirror_ghosts_x(grid, field, Stagger::STAGGERED, false, zero, bc_no_post, 0, field.ny());
  mirror_ghosts_x(grid, field, Stagger::STAGGERED, true, zero, bc_no_post, 0, field.ny());
  mirror_ghosts_y(grid,
                  field,
                  Stagger::STAGGERED,
                  false,
                  zero,
                  bc_no_post,
                  -field.nghost(),
                  field.nx() + field.nghost());
  mirror_ghosts_y(grid,
                  field,
                  Stagger::STAGGERED,
                  true,
                  zero,
                  bc_no_post,
                  -field.nghost(),
                  field.nx() + field.nghost());
}

// -------------------------------------------------------------------------------------------------
// How a Dirichlet condition is realized on a cell-centered scalar.
//   CONSTANT    : write the boundary value into every ghost cell. This is what an upwinded inflow
//                 flux wants (the flux across the boundary face then carries exactly `value`) but
//                 it only represents the boundary value to first order for gradients and
//                 reconstructions.
//   EXTRAPOLATE : reconstruct the ghosts so that the interpolant hits `value` at the wall. Second
//                 order (or `order + 1`), but the inflow flux no longer carries `value` exactly.
enum class DirichletFill { CONSTANT, EXTRAPOLATE };

template <typename Float, Layout LAYOUT>
constexpr void apply_dirichlet_bconds(const Grid<Float, LAYOUT>& grid,
                                      Scalar<Float, LAYOUT> field,
                                      Float value,
                                      DirichletFill fill = DirichletFill::CONSTANT,
                                      Index order        = 1) noexcept {
  IGOR_ASSERT(field.nghost() > 0, "Expected at least one ghost cell, but got {}", field.nghost());

  if (fill == DirichletFill::CONSTANT) {
    grid.foreach_range(
        0, 1, -field.nghost(), field.ny() + field.nghost(), [=](Index, Index j) {
          for (Index i = -field.nghost(); i < 0; ++i) {
            field(i, j) = value;
          }
          for (Index i = field.nx(); i < field.nx() + field.nghost(); ++i) {
            field(i, j) = value;
          }
        });
    grid.foreach_range(
        -field.nghost(), field.nx() + field.nghost(), 0, 1, [=](Index i, Index) {
          for (Index j = -field.nghost(); j < 0; ++j) {
            field(i, j) = value;
          }
          for (Index j = field.ny(); j < field.ny() + field.nghost(); ++j) {
            field(i, j) = value;
          }
        });
    return;
  }

  const auto s =
      make_bc_stencil<Float>(BCKind::DIRICHLET, Stagger::STAGGERED, field.nghost(), order);
  const auto bc = [=](Index) -> Float { return value; };

  fill_ghosts_x(grid, field, Stagger::STAGGERED, false, s, bc, bc_no_post, 0, field.ny());
  fill_ghosts_x(grid, field, Stagger::STAGGERED, true, s, bc, bc_no_post, 0, field.ny());
  fill_ghosts_y(grid,
                field,
                Stagger::STAGGERED,
                false,
                s,
                bc,
                bc_no_post,
                -field.nghost(),
                field.nx() + field.nghost());
  fill_ghosts_y(grid,
                field,
                Stagger::STAGGERED,
                true,
                s,
                bc,
                bc_no_post,
                -field.nghost(),
                field.nx() + field.nghost());
}
