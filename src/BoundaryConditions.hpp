#pragma once

#include <variant>

#include "Grid.hpp"

// =================================================================================================
template <typename Float>
struct Dirichlet {
  std::variant<Float, Float (*)(Float, Float)> val;

  // = LEFT ========================================================================================
  template <Layout LAYOUT>
  static constexpr void
  apply_left_offset(const Grid<Float, LAYOUT>& grid,
                    Scalar<Float, LAYOUT> s,
                    bool use_ym,
                    Float t,
                    std::variant<Float, Float (*)(Float, Float)> value) noexcept {
    grid.foreach_range(
        0, 1, 0, s.ny(), FOREACH_FUNC {
          // Linear extrapolation
          const auto s0 = s(0, j);
          const auto y  = use_ym ? grid.ym(j) : grid.y(j);
          const auto v =
              std::holds_alternative<Float>(value) ? std::get<0>(value) : std::get<1>(value)(y, t);
          for (i = -s.nghost(); i < 0; ++i) {
            s(i, j) = s0 - 2.0 * (v - s0) * i;
          }
        });
  }

  template <Layout LAYOUT>
  static constexpr void
  apply_left_align(const Grid<Float, LAYOUT>& grid,
                   Scalar<Float, LAYOUT> s,
                   bool use_ym,
                   Float t,
                   std::variant<Float, Float (*)(Float, Float)> value) noexcept {
    grid.foreach_range(
        0, 1, 0, s.ny(), FOREACH_FUNC {
          // Linear extrapolation
          const auto s1 = s(1, j);
          const auto y  = use_ym ? grid.ym(j) : grid.y(j);
          const auto v =
              std::holds_alternative<Float>(value) ? std::get<0>(value) : std::get<1>(value)(y, t);
          for (i = -s.nghost(); i <= 0; ++i) {
            s(i, j) = s1 - (v - s1) * (i - 1);
          }
        });
  }

  // = RIGHT =======================================================================================
  template <Layout LAYOUT>
  static constexpr void
  apply_right_offset(const Grid<Float, LAYOUT>& grid,
                     Scalar<Float, LAYOUT> s,
                     bool use_ym,
                     Float t,
                     std::variant<Float, Float (*)(Float, Float)> value) noexcept {
    grid.foreach_range(
        0, 1, 0, s.ny(), FOREACH_FUNC {
          // Linear extrapolation
          const auto sN = s(s.nx() - 1, j);
          const auto y  = use_ym ? grid.ym(j) : grid.y(j);
          const auto v =
              std::holds_alternative<Float>(value) ? std::get<0>(value) : std::get<1>(value)(y, t);
          for (i = s.nx(); i < s.nx() + s.nghost(); ++i) {
            s(i, j) = sN + 2.0 * (v - sN) * (i - s.nx() + 1);
          }
        });
  }

  template <Layout LAYOUT>
  static constexpr void
  apply_right_align(const Grid<Float, LAYOUT>& grid,
                    Scalar<Float, LAYOUT> s,
                    bool use_ym,
                    Float t,
                    std::variant<Float, Float (*)(Float, Float)> value) noexcept {
    grid.foreach_range(
        0, 1, 0, s.ny(), FOREACH_FUNC {
          // Linear extrapolation
          const auto sN = s(s.nx() - 2, j);
          const auto y  = use_ym ? grid.ym(j) : grid.y(j);
          const auto v =
              std::holds_alternative<Float>(value) ? std::get<0>(value) : std::get<1>(value)(y, t);
          for (i = s.nx() - 1; i < s.nx() + s.nghost(); ++i) {
            s(i, j) = sN + (v - sN) * (i - s.nx() + 2);
          }
        });
  }

  // = BOTTOM ======================================================================================
  template <Layout LAYOUT>
  static constexpr void
  apply_bottom_offset(const Grid<Float, LAYOUT>& grid,
                      Scalar<Float, LAYOUT> s,
                      bool use_xm,
                      Float t,
                      std::variant<Float, Float (*)(Float, Float)> value) noexcept {
    grid.foreach_range(
        0, s.nx(), 0, 1, FOREACH_FUNC {
          // Linear extrapolation
          const auto s0 = s(i, 0);
          const auto x  = use_xm ? grid.xm(i) : grid.x(i);
          const auto v =
              std::holds_alternative<Float>(value) ? std::get<0>(value) : std::get<1>(value)(x, t);
          for (j = -s.nghost(); j < 0; ++j) {
            s(i, j) = s0 - 2.0 * (v - s0) * j;
          }
        });
  }

  template <Layout LAYOUT>
  static constexpr void
  apply_bottom_align(const Grid<Float, LAYOUT>& grid,
                     Scalar<Float, LAYOUT> s,
                     bool use_xm,
                     Float t,
                     std::variant<Float, Float (*)(Float, Float)> value) noexcept {
    grid.foreach_range(
        0, s.nx(), 0, 1, FOREACH_FUNC {
          // Linear extrapolation
          const auto s1 = s(i, 1);
          const auto x  = use_xm ? grid.xm(i) : grid.x(i);
          const auto v =
              std::holds_alternative<Float>(value) ? std::get<0>(value) : std::get<1>(value)(x, t);
          for (j = -s.nghost(); j <= 0; ++j) {
            s(i, j) = s1 - (v - s1) * (j - 1);
          }
        });
  }

  // = TOP =========================================================================================
  template <Layout LAYOUT>
  static constexpr void
  apply_top_offset(const Grid<Float, LAYOUT>& grid,
                   Scalar<Float, LAYOUT> s,
                   bool use_xm,
                   Float t,
                   std::variant<Float, Float (*)(Float, Float)> value) noexcept {
    grid.foreach_range(
        0, s.nx(), 0, 1, FOREACH_FUNC {
          // Linear extrapolation
          const auto sN = s(i, s.ny() - 1);
          const auto x  = use_xm ? grid.xm(i) : grid.x(i);
          const auto v =
              std::holds_alternative<Float>(value) ? std::get<0>(value) : std::get<1>(value)(x, t);
          for (j = s.ny(); j < s.ny() + s.nghost(); ++j) {
            s(i, j) = sN + 2.0 * (v - sN) * (j - s.ny() + 1);
          }
        });
  }

  template <Layout LAYOUT>
  static constexpr void
  apply_top_align(const Grid<Float, LAYOUT>& grid,
                  Scalar<Float, LAYOUT> s,
                  bool use_xm,
                  Float t,
                  std::variant<Float, Float (*)(Float, Float)> value) noexcept {
    grid.foreach_range(
        0, s.nx(), 0, 1, FOREACH_FUNC {
          // Linear extrapolation
          const auto sN = s(i, s.ny() - 2);
          const auto x  = use_xm ? grid.xm(j) : grid.x(j);
          const auto v =
              std::holds_alternative<Float>(value) ? std::get<0>(value) : std::get<1>(value)(x, t);
          for (j = s.ny() - 1; j < s.ny() + s.nghost(); ++j) {
            s(i, j) = sN + (v - sN) * (j - s.ny() + 2);
          }
        });
  }
};

// =================================================================================================
struct Neumann {

  // = LEFT ========================================================================================
  template <typename Float, Layout LAYOUT>
  static constexpr void apply_left_offset(const Grid<Float, LAYOUT>& grid,
                                          Scalar<Float, LAYOUT> s) noexcept {
    grid.foreach_range(
        0, 1, 0, s.ny(), FOREACH_FUNC {
          for (i = -s.nghost(); i < 0; ++i) {
            s(i, j) = s(-i - 1, j);
          }
        });
  }

  template <typename Float, Layout LAYOUT>
  static constexpr void apply_left_align(const Grid<Float, LAYOUT>& grid,
                                         Scalar<Float, LAYOUT> s) noexcept {
    return apply_left_offset(grid, s);  // The same, maybe not in the future
  }

  // = RIGHT =======================================================================================
  template <typename Float, Layout LAYOUT>
  static constexpr void apply_right_offset(const Grid<Float, LAYOUT>& grid,
                                           Scalar<Float, LAYOUT> s) noexcept {
    grid.foreach_range(
        0, 1, 0, s.ny(), FOREACH_FUNC {
          for (i = s.nx(); i < s.nx() + s.nghost(); ++i) {
            s(i, j) = s(2 * s.nx() - i - 1, j);
          }
        });
  }

  template <typename Float, Layout LAYOUT>
  static constexpr void apply_right_align(const Grid<Float, LAYOUT>& grid,
                                          Scalar<Float, LAYOUT> s) noexcept {
    return apply_right_offset(grid, s);  // The same, maybe not in the future
  }

  // = BOTTOM ======================================================================================
  template <typename Float, Layout LAYOUT>
  static constexpr void apply_bottom_offset(const Grid<Float, LAYOUT>& grid,
                                            Scalar<Float, LAYOUT> s) noexcept {
    grid.foreach_range(
        0, s.nx(), 0, 1, FOREACH_FUNC {
          for (j = -s.nghost(); j < 0; ++j) {
            s(i, j) = s(i, -j - 1);
          }
        });
  }

  template <typename Float, Layout LAYOUT>
  static constexpr void apply_bottom_align(const Grid<Float, LAYOUT>& grid,
                                           Scalar<Float, LAYOUT> s) noexcept {
    return apply_bottom_offset(grid, s);
  }

  // = TOP =========================================================================================
  template <typename Float, Layout LAYOUT>
  static constexpr void apply_top_offset(const Grid<Float, LAYOUT>& grid,
                                         Scalar<Float, LAYOUT> s) noexcept {
    grid.foreach_range(
        0, s.nx(), 0, 1, FOREACH_FUNC {
          for (j = s.ny(); j < s.ny() + s.nghost(); ++j) {
            s(i, j) = s(i, 2 * s.ny() - j - 1);
          }
        });
  }

  template <typename Float, Layout LAYOUT>
  static constexpr void apply_top_align(const Grid<Float, LAYOUT>& grid,
                                        Scalar<Float, LAYOUT> s) noexcept {
    return apply_top_offset(grid, s);
  }
};

// =================================================================================================
struct Periodic {

  // = LEFT ========================================================================================
  template <typename Float, Layout LAYOUT>
  static constexpr void apply_left_offset(const Grid<Float, LAYOUT>& grid,
                                          Scalar<Float, LAYOUT> s) noexcept {
    grid.foreach_range(
        0, 1, 0, s.ny(), FOREACH_FUNC {
          for (i = -s.nghost(); i < 0; ++i) {
            s(i, j) = s(s.nx() + i, j);
          }
        });
  }

  template <typename Float, Layout LAYOUT>
  static constexpr void apply_left_align(const Grid<Float, LAYOUT>& grid,
                                         Scalar<Float, LAYOUT> s) noexcept {
    grid.foreach_range(
        0, 1, 0, s.ny(), FOREACH_FUNC {
          for (i = -s.nghost(); i < 0; ++i) {
            s(i, j) = s(s.nx() + i - 1, j);
          }
        });
  }

  // = RIGHT =======================================================================================
  template <typename Float, Layout LAYOUT>
  static constexpr void apply_right_offset(const Grid<Float, LAYOUT>& grid,
                                           Scalar<Float, LAYOUT> s) noexcept {
    grid.foreach_range(
        0, 1, 0, s.ny(), FOREACH_FUNC {
          for (i = s.nx(); i < s.nx() + s.nghost(); ++i) {
            s(i, j) = s(i - s.nx(), j);
          }
        });
  }

  template <typename Float, Layout LAYOUT>
  static constexpr void apply_right_align(const Grid<Float, LAYOUT>& grid,
                                          Scalar<Float, LAYOUT> s) noexcept {
    grid.foreach_range(
        0, 1, 0, s.ny(), FOREACH_FUNC {
          for (i = s.nx(); i < s.nx() + s.nghost(); ++i) {
            s(i, j) = s(i - s.nx() + 1, j);
          }
        });
  }

  // = BOTTOM ======================================================================================
  template <typename Float, Layout LAYOUT>
  static constexpr void apply_bottom_offset(const Grid<Float, LAYOUT>& grid,
                                            Scalar<Float, LAYOUT> s) noexcept {
    grid.foreach_range(
        0, s.nx(), 0, 1, FOREACH_FUNC {
          for (j = -s.nghost(); j < 0; ++j) {
            s(i, j) = s(i, s.ny() + j);
          }
        });
  }

  template <typename Float, Layout LAYOUT>
  static constexpr void apply_bottom_align(const Grid<Float, LAYOUT>& grid,
                                           Scalar<Float, LAYOUT> s) noexcept {
    grid.foreach_range(
        0, s.nx(), 0, 1, FOREACH_FUNC {
          for (j = -s.nghost(); j < 0; ++j) {
            s(i, j) = s(i, s.ny() + j - 1);
          }
        });
  }

  // = TOP =========================================================================================
  template <typename Float, Layout LAYOUT>
  static constexpr void apply_top_offset(const Grid<Float, LAYOUT>& grid,
                                         Scalar<Float, LAYOUT> s) noexcept {
    grid.foreach_range(
        0, s.nx(), 0, 1, FOREACH_FUNC {
          for (j = s.ny(); j < s.ny() + s.nghost(); ++j) {
            s(i, j) = s(i, s.ny() - j);
          }
        });
  }

  template <typename Float, Layout LAYOUT>
  static constexpr void apply_top_align(const Grid<Float, LAYOUT>& grid,
                                        Scalar<Float, LAYOUT> s) noexcept {
    grid.foreach_range(
        0, s.nx(), 0, 1, FOREACH_FUNC {
          for (j = s.ny(); j < s.ny() + s.nghost(); ++j) {
            s(i, j) = s(i, s.ny() - j + 1);
          }
        });
  }
};

// -------------------------------------------------------------------------------------------------
template <typename Float>
using BCond_t = std::variant<Dirichlet<Float>, Neumann, Periodic>;

template <typename Float>
struct BConds {
  BCond_t<Float> left;
  BCond_t<Float> right;
  BCond_t<Float> bottom;
  BCond_t<Float> top;
};

// -------------------------------------------------------------------------------------------------
template <typename Float, Layout LAYOUT>
constexpr void apply_bconds(const Grid<Float, LAYOUT>& grid,
                            const BConds<Float> bconds,
                            Scalar<Float, LAYOUT> s,
                            Float t) noexcept {
#define APPLY(side)                                                                                \
  do {                                                                                             \
    if (std::holds_alternative<Dirichlet<Float>>(bconds.side)) {                                   \
      Dirichlet<Float>::apply_##side##_offset(grid, s, true, t, std::get<0>(bconds.side).val);     \
    } else if (std::holds_alternative<Neumann>(bconds.side)) {                                     \
      Neumann::apply_##side##_offset(grid, s);                                                     \
    } else if (std::holds_alternative<Periodic>(bconds.side)) {                                    \
      Periodic::apply_##side##_offset(grid, s);                                                    \
    } else {                                                                                       \
      Igor::Panic("Unreachable");                                                                  \
    }                                                                                              \
  } while (false)

  APPLY(left);
  APPLY(right);
  APPLY(bottom);
  APPLY(top);

#undef APPLY
}

// -------------------------------------------------------------------------------------------------
template <typename Float, Layout LAYOUT>
void apply_velocity_bconds(const Grid<Float, LAYOUT>& grid,
                           const BConds<Float>& u_bconds,
                           const BConds<Float>& v_bconds,
                           FaceVector<Float, LAYOUT> u,
                           Float t = -1.0) {
#define APPLY_U_X(side)                                                                            \
  do {                                                                                             \
    if (std::holds_alternative<Dirichlet<Float>>(u_bconds.side)) {                                 \
      Dirichlet<Float>::apply_##side##_align(grid, u.x, true, t, std::get<0>(u_bconds.side).val);  \
    } else if (std::holds_alternative<Neumann>(u_bconds.side)) {                                   \
      Neumann::apply_##side##_align(grid, u.x);                                                    \
    } else if (std::holds_alternative<Periodic>(u_bconds.side)) {                                  \
      Periodic::apply_##side##_align(grid, u.x);                                                   \
    } else {                                                                                       \
      Igor::Panic("Unreachable");                                                                  \
    }                                                                                              \
  } while (false)

#define APPLY_V_X(side)                                                                            \
  do {                                                                                             \
    if (std::holds_alternative<Dirichlet<Float>>(v_bconds.side)) {                                 \
      Dirichlet<Float>::apply_##side##_offset(                                                     \
          grid, u.y, false, t, std::get<0>(v_bconds.side).val);                                    \
    } else if (std::holds_alternative<Neumann>(v_bconds.side)) {                                   \
      Neumann::apply_##side##_offset(grid, u.y);                                                   \
    } else if (std::holds_alternative<Periodic>(v_bconds.side)) {                                  \
      Periodic::apply_##side##_offset(grid, u.y);                                                  \
    } else {                                                                                       \
      Igor::Panic("Unreachable");                                                                  \
    }                                                                                              \
  } while (false)

#define APPLY_U_Y(side)                                                                            \
  do {                                                                                             \
    if (std::holds_alternative<Dirichlet<Float>>(u_bconds.side)) {                                 \
      Dirichlet<Float>::apply_##side##_offset(                                                     \
          grid, u.x, false, t, std::get<0>(u_bconds.side).val);                                    \
    } else if (std::holds_alternative<Neumann>(u_bconds.side)) {                                   \
      Neumann::apply_##side##_offset(grid, u.x);                                                   \
    } else if (std::holds_alternative<Periodic>(u_bconds.side)) {                                  \
      Periodic::apply_##side##_offset(grid, u.x);                                                  \
    } else {                                                                                       \
      Igor::Panic("Unreachable");                                                                  \
    }                                                                                              \
  } while (false)

#define APPLY_V_Y(side)                                                                            \
  do {                                                                                             \
    if (std::holds_alternative<Dirichlet<Float>>(v_bconds.side)) {                                 \
      Dirichlet<Float>::apply_##side##_align(grid, u.y, true, t, std::get<0>(v_bconds.side).val);  \
    } else if (std::holds_alternative<Neumann>(v_bconds.side)) {                                   \
      Neumann::apply_##side##_align(grid, u.y);                                                    \
    } else if (std::holds_alternative<Periodic>(v_bconds.side)) {                                  \
      Periodic::apply_##side##_align(grid, u.y);                                                   \
    } else {                                                                                       \
      Igor::Panic("Unreachable");                                                                  \
    }                                                                                              \
  } while (false)

  APPLY_U_X(left);
  APPLY_U_X(right);
  APPLY_U_Y(bottom);
  APPLY_U_Y(top);

  APPLY_V_X(left);
  APPLY_V_X(right);
  APPLY_V_Y(bottom);
  APPLY_V_Y(top);

#undef APPLY_U_X
#undef APPLY_U_Y
#undef APPLY_V_X
#undef APPLY_V_Y
}

// -------------------------------------------------------------------------------------------------
template <typename Float, Layout LAYOUT>
constexpr void apply_dirichlet_bconds(const Grid<Float, LAYOUT>& grid,
                                      Scalar<Float, LAYOUT> field,
                                      Float value) noexcept {
  IGOR_ASSERT(field.nghost() > 0, "Expected at least one ghost cell, but got {}", field.nghost());
  Dirichlet<Float>::apply_left_offset(grid, field, true, -1.0, value);
  Dirichlet<Float>::apply_right_offset(grid, field, true, -1.0, value);
  Dirichlet<Float>::apply_bottom_offset(grid, field, true, -1.0, value);
  Dirichlet<Float>::apply_top_offset(grid, field, true, -1.0, value);
}

// -------------------------------------------------------------------------------------------------
template <typename Float, Layout LAYOUT>
constexpr void apply_neumann_bconds(const Grid<Float, LAYOUT>& grid,
                                    Scalar<Float, LAYOUT> field) noexcept {
  IGOR_ASSERT(field.nghost() > 0, "Expected at least one ghost cell, but got {}", field.nghost());
  Neumann::apply_left_offset(grid, field);
  Neumann::apply_right_offset(grid, field);
  Neumann::apply_bottom_offset(grid, field);
  Neumann::apply_top_offset(grid, field);
}

// -------------------------------------------------------------------------------------------------
template <typename Float, Layout LAYOUT>
constexpr void apply_periodic_bconds(const Grid<Float, LAYOUT>& grid,
                                     Scalar<Float, LAYOUT> field) noexcept {
  IGOR_ASSERT(field.nghost() > 0, "Expected at least one ghost cell, but got {}", field.nghost());
  Periodic::apply_left_offset(grid, field);
  Periodic::apply_right_offset(grid, field);
  Periodic::apply_bottom_offset(grid, field);
  Periodic::apply_top_offset(grid, field);
}
