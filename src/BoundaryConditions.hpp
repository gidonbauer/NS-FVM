#pragma once

#include <variant>

#include "Grid.hpp"

// -------------------------------------------------------------------------------------------------
template <typename Float>
struct Dirichlet {
  std::variant<Float, Float (*)(Float, Float)> U{}, V{};

  // = LEFT ========
  template <Layout LAYOUT>
  constexpr void apply_left(const Grid<Float, LAYOUT>& grid,
                            FaceVector<Float, LAYOUT> uf,
                            Float t) const noexcept {
    IGOR_ASSERT(grid.nghost() >= 1, "Expected at least one ghost cell, but got {}", grid.nghost());

    grid.foreach_range(
        -uf.x.nghost(), 1, -uf.x.nghost(), uf.x.ny() + uf.x.nghost(), [=, this](Index i, Index j) {
          const Float U_bc =
              std::holds_alternative<Float>(U) ? std::get<0>(U) : std::get<1>(U)(grid.ym(j), t);
          uf.x(i, j) = U_bc;
        });
    grid.foreach_range(
        -uf.y.nghost(), 0, -uf.y.nghost(), uf.y.ny() + uf.y.nghost(), [=, this](Index i, Index j) {
          const Float V_bc =
              std::holds_alternative<Float>(V) ? std::get<0>(V) : std::get<1>(V)(grid.y(j), t);
          if (i == -1) {
            uf.y(i, j) = 2.0 * V_bc - uf.y(0, j);
          } else {
            uf.y(i, j) = V_bc;
          }
        });
  }

  // = RIGHT =======
  template <Layout LAYOUT>
  constexpr void apply_right(const Grid<Float, LAYOUT>& grid,
                             FaceVector<Float, LAYOUT> uf,
                             Float t) const noexcept {
    IGOR_ASSERT(grid.nghost() == 1, "Expected exactly one ghost cell, but got {}", grid.nghost());
    grid.foreach_range(0, 1, -uf.x.nghost(), uf.x.ny() + uf.x.nghost(), [=, this](Index, Index j) {
      const Float U_bc =
          std::holds_alternative<Float>(U) ? std::get<0>(U) : std::get<1>(U)(grid.ym(j), t);
      uf.x(grid.nx(), j)     = U_bc;
      uf.x(grid.nx() + 1, j) = U_bc;
    });
    grid.foreach_range(0, 1, -uf.y.nghost(), uf.y.ny() + uf.y.nghost(), [=, this](Index, Index j) {
      const Float V_bc =
          std::holds_alternative<Float>(V) ? std::get<0>(V) : std::get<1>(V)(grid.y(j), t);
      uf.y(grid.nx(), j) = 2.0 * V_bc - uf.y(grid.nx() - 1, j);
    });
  }

  // = BOTTOM ======
  template <Layout LAYOUT>
  constexpr void apply_bottom(const Grid<Float, LAYOUT>& grid,
                              FaceVector<Float, LAYOUT> uf,
                              Float t) const noexcept {
    IGOR_ASSERT(grid.nghost() == 1, "Expected exactly one ghost cell, but got {}", grid.nghost());
    grid.foreach_range(-uf.x.nghost(), uf.x.nx() + uf.x.nghost(), 0, 1, [=, this](Index i, Index) {
      const Float U_bc =
          std::holds_alternative<Float>(U) ? std::get<0>(U) : std::get<1>(U)(grid.x(i), t);
      uf.x(i, -1) = 2.0 * U_bc - uf.x(i, 0);
    });
    grid.foreach_range(-uf.y.nghost(), uf.y.nx() + uf.y.nghost(), 0, 1, [=, this](Index i, Index) {
      const Float V_bc =
          std::holds_alternative<Float>(V) ? std::get<0>(V) : std::get<1>(V)(grid.xm(i), t);
      uf.y(i, -1) = V_bc;
      uf.y(i, 0)  = V_bc;
    });
  }

  // = TOP =========
  template <Layout LAYOUT>
  constexpr void
  apply_top(const Grid<Float, LAYOUT>& grid, FaceVector<Float, LAYOUT> uf, Float t) const noexcept {
    IGOR_ASSERT(grid.nghost() == 1, "Expected exactly one ghost cell, but got {}", grid.nghost());
    grid.foreach_range(-uf.x.nghost(), uf.x.nx() + uf.x.nghost(), 0, 1, [=, this](Index i, Index) {
      const Float U_bc =
          std::holds_alternative<Float>(U) ? std::get<0>(U) : std::get<1>(U)(grid.x(i), t);
      uf.x(i, grid.ny()) = 2.0 * U_bc - uf.x(i, grid.ny() - 1);
    });
    grid.foreach_range(-uf.y.nghost(), uf.y.nx() + uf.y.nghost(), 0, 1, [=, this](Index i, Index) {
      const Float V_bc =
          std::holds_alternative<Float>(V) ? std::get<0>(V) : std::get<1>(V)(grid.xm(i), t);
      uf.y(i, grid.ny())     = V_bc;
      uf.y(i, grid.ny() + 1) = V_bc;
    });
  }
};

// -------------------------------------------------------------------------------------------------
struct Neumann {
  bool clipped = false;

  template <typename Float, Layout LAYOUT>
  constexpr void apply_left(const Grid<Float, LAYOUT>& grid,
                            FaceVector<Float, LAYOUT> uf,
                            Float /*t*/) const noexcept {
    IGOR_ASSERT(grid.nghost() == 1, "Expected exactly one ghost cell, but got {}", grid.nghost());
    if (clipped) {
      grid.foreach_range(0, 1, -uf.x.nghost(), uf.x.ny() + uf.x.nghost(), [=](Index, Index j) {
        uf.x(-1, j) = std::min(uf.x(0, j), 0.0);
      });
    } else {
      grid.foreach_range(0, 1, -uf.x.nghost(), uf.x.ny() + uf.x.nghost(), [=](Index, Index j) {
        uf.x(-1, j) = uf.x(0, j);
      });
    }
    grid.foreach_range(0, 1, -uf.y.nghost(), uf.y.ny() + uf.y.nghost(), [=](Index, Index j) {
      uf.y(-1, j) = uf.y(0, j);
    });
  }

  template <typename Float, Layout LAYOUT>
  constexpr void apply_right(const Grid<Float, LAYOUT>& grid,
                             FaceVector<Float, LAYOUT> uf,
                             Float /*t*/) const noexcept {
    IGOR_ASSERT(grid.nghost() == 1, "Expected exactly one ghost cell, but got {}", grid.nghost());
    if (clipped) {
      grid.foreach_range(0, 1, -uf.x.nghost(), uf.x.ny() + uf.x.nghost(), [=](Index, Index j) {
        uf.x(grid.nx() + 1, j) = std::max(uf.x(grid.nx(), j), 0.0);
      });
    } else {
      grid.foreach_range(0, 1, -uf.x.nghost(), uf.x.ny() + uf.x.nghost(), [=](Index, Index j) {
        uf.x(grid.nx() + 1, j) = uf.x(grid.nx(), j);
      });
    }
    grid.foreach_range(0, 1, -uf.y.nghost(), uf.y.ny() + uf.y.nghost(), [=](Index, Index j) {
      uf.y(grid.nx(), j) = uf.y(grid.nx() - 1, j);
    });
  }

  template <typename Float, Layout LAYOUT>
  constexpr void apply_bottom(const Grid<Float, LAYOUT>& grid,
                              FaceVector<Float, LAYOUT> uf,
                              Float /*t*/) const noexcept {
    IGOR_ASSERT(grid.nghost() == 1, "Expected exactly one ghost cell, but got {}", grid.nghost());
    grid.foreach_range(-uf.x.nghost(), uf.x.nx() + uf.x.nghost(), 0, 1, [=](Index i, Index) {
      uf.x(i, -1) = uf.x(i, 0);
    });
    if (clipped) {
      grid.foreach_range(-uf.y.nghost(), uf.y.nx() + uf.y.nghost(), 0, 1, [=](Index i, Index) {
        uf.y(i, -1) = std::min(uf.y(i, 0), 0.0);
      });
    } else {
      grid.foreach_range(-uf.y.nghost(), uf.y.nx() + uf.y.nghost(), 0, 1, [=](Index i, Index) {
        uf.y(i, -1) = uf.y(i, 0);
      });
    }
  }

  template <typename Float, Layout LAYOUT>
  constexpr void apply_top(const Grid<Float, LAYOUT>& grid,
                           FaceVector<Float, LAYOUT> uf,
                           Float /*t*/) const noexcept {
    IGOR_ASSERT(grid.nghost() == 1, "Expected exactly one ghost cell, but got {}", grid.nghost());
    grid.foreach_range(-uf.x.nghost(), uf.x.nx() + uf.x.nghost(), 0, 1, [=](Index i, Index) {
      uf.x(i, grid.ny()) = uf.x(i, grid.ny() - 1);
    });
    if (clipped) {
      grid.foreach_range(-uf.y.nghost(), uf.y.nx() + uf.y.nghost(), 0, 1, [=](Index i, Index) {
        uf.y(i, grid.ny() + 1) = std::max(uf.y(i, grid.ny()), 0.0);
      });
    } else {
      grid.foreach_range(-uf.y.nghost(), uf.y.nx() + uf.y.nghost(), 0, 1, [=](Index i, Index) {
        uf.y(i, grid.ny() + 1) = uf.y(i, grid.ny());
      });
    }
  }
};

// -------------------------------------------------------------------------------------------------
struct Periodic {
  template <typename Float, Layout LAYOUT>
  constexpr void apply_left(const Grid<Float, LAYOUT>& grid,
                            FaceVector<Float, LAYOUT> uf,
                            Float /*t*/) const noexcept {
    IGOR_ASSERT(grid.nghost() == 1, "Expected exactly one ghost cell, but got {}", grid.nghost());
    grid.foreach_range(0, 1, -uf.x.nghost(), uf.x.ny() + uf.x.nghost(), [=](Index, Index j) {
      uf.x(-1, j) = uf.x(grid.nx() - 1, j);
    });
    grid.foreach_range(0, 1, -uf.y.nghost(), uf.y.ny() + uf.y.nghost(), [=](Index, Index j) {
      uf.y(-1, j) = uf.y(grid.nx() - 1, j);
    });
  }

  template <typename Float, Layout LAYOUT>
  constexpr void apply_right(const Grid<Float, LAYOUT>& grid,
                             FaceVector<Float, LAYOUT> uf,
                             Float /*t*/) const noexcept {
    IGOR_ASSERT(grid.nghost() == 1, "Expected exactly one ghost cell, but got {}", grid.nghost());
    grid.foreach_range(0, 1, -uf.x.nghost(), uf.x.ny() + uf.x.nghost(), [=](Index, Index j) {
      uf.x(grid.nx() + 1, j) = uf.x(1, j);
    });
    grid.foreach_range(0, 1, -uf.y.nghost(), uf.y.ny() + uf.y.nghost(), [=](Index, Index j) {
      uf.y(grid.nx(), j) = uf.y(0, j);
    });
  }

  template <typename Float, Layout LAYOUT>
  constexpr void apply_bottom(const Grid<Float, LAYOUT>& grid,
                              FaceVector<Float, LAYOUT> uf,
                              Float /*t*/) const noexcept {
    IGOR_ASSERT(grid.nghost() == 1, "Expected exactly one ghost cell, but got {}", grid.nghost());
    grid.foreach_range(-uf.x.nghost(), uf.x.nx() + uf.x.nghost(), 0, 1, [=](Index i, Index) {
      uf.x(i, -1) = uf.x(i, grid.ny() - 1);
    });
    grid.foreach_range(-uf.y.nghost(), uf.y.nx() + uf.y.nghost(), 0, 1, [=](Index i, Index) {
      uf.y(i, -1) = uf.y(i, grid.ny() - 1);
    });
  }

  template <typename Float, Layout LAYOUT>
  constexpr void apply_top(const Grid<Float, LAYOUT>& grid,
                           FaceVector<Float, LAYOUT> uf,
                           Float /*t*/) const noexcept {
    IGOR_ASSERT(grid.nghost() == 1, "Expected exactly one ghost cell, but got {}", grid.nghost());
    grid.foreach_range(-uf.x.nghost(), uf.x.nx() + uf.x.nghost(), 0, 1, [=](Index i, Index) {
      uf.x(i, grid.ny()) = uf.x(i, 0);
    });
    grid.foreach_range(-uf.y.nghost(), uf.y.nx() + uf.y.nghost(), 0, 1, [=](Index i, Index) {
      uf.y(i, grid.ny() + 1) = uf.y(i, 1);
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
  std::visit([&](auto&& bcond) { bcond.apply_left(grid, uf, t); }, bconds.left);
  std::visit([&](auto&& bcond) { bcond.apply_right(grid, uf, t); }, bconds.right);
  std::visit([&](auto&& bcond) { bcond.apply_bottom(grid, uf, t); }, bconds.bottom);
  std::visit([&](auto&& bcond) { bcond.apply_top(grid, uf, t); }, bconds.top);
}

// -------------------------------------------------------------------------------------------------
template <typename Float, Layout LAYOUT>
constexpr void apply_neumann_bconds(const Grid<Float, LAYOUT>& grid,
                                    Scalar<Float, LAYOUT> field) noexcept {
  IGOR_ASSERT(field.nghost() > 0, "Expected at least one ghost cell, but got {}", field.nghost());

  grid.foreach_range(0, 1, -field.nghost(), field.ny() + field.nghost(), [=](Index, Index j) {
    // LEFT
    for (Index i = -field.nghost(); i < 0; ++i) {
      field(i, j) = field(0, j);
    }
    // RIGHT
    for (Index i = field.nx(); i < field.nx() + field.nghost(); ++i) {
      field(i, j) = field(field.nx() - 1, j);
    }
  });

  grid.foreach_range(-field.nghost(), field.nx() + field.nghost(), 0, 1, [=](Index i, Index) {
    // BOTTOM
    for (Index j = -field.nghost(); j < 0; ++j) {
      field(i, j) = field(i, 0);
    }
    // TOP
    for (Index j = field.ny(); j < field.ny() + field.nghost(); ++j) {
      field(i, j) = field(i, field.ny() - 1);
    }
  });
}

// -------------------------------------------------------------------------------------------------
template <typename Float, Layout LAYOUT>
constexpr void apply_dirichlet_bconds(const Grid<Float, LAYOUT>& grid,
                                      Scalar<Float, LAYOUT> field,
                                      Float value) noexcept {
  IGOR_ASSERT(field.nghost() > 0, "Expected at least one ghost cell, but got {}", field.nghost());

  grid.foreach_range(0, 1, -field.nghost(), field.ny() + field.nghost(), [=](Index, Index j) {
    // LEFT
    for (Index i = -field.nghost(); i < 0; ++i) {
      field(i, j) = value;
    }
    // RIGHT
    for (Index i = field.nx(); i < field.nx() + field.nghost(); ++i) {
      field(i, j) = value;
    }
  });

  grid.foreach_range(-field.nghost(), field.nx() + field.nghost(), 0, 1, [=](Index i, Index) {
    // BOTTOM
    for (Index j = -field.nghost(); j < 0; ++j) {
      field(i, j) = value;
    }
    // TOP
    for (Index j = field.ny(); j < field.ny() + field.nghost(); ++j) {
      field(i, j) = value;
    }
  });
}
