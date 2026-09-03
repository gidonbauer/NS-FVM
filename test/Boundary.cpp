// Tests for the ghost-layer boundary conditions.
//
// Four groups of checks:
//   1. Polynomial exactness    -- a field that is a polynomial of degree `order` along the normal
//                                 must be reproduced *exactly* in every ghost cell, corners
//                                 included. This pins the extrapolation weights.
//   2. Discrete boundary data  -- the operators the solver actually evaluates at the wall (the
//                                 face value itself, the two-point tangential average, the
//                                 normal difference) must return the prescribed boundary data.
//   3. Structural invariants   -- Neumann is a mirror, Periodic is a wrap, and the first ghost
//                                 layer does not depend on how many ghost layers exist.
//   4. Order of accuracy       -- for a non-polynomial field the ghost error must decay like
//                                 h^(order+1).

#include <array>
#include <cmath>
#include <limits>
#include <numbers>

#include "BoundaryConditions.hpp"

using Float = double;

constexpr Float pi = std::numbers::pi_v<Float>;

constexpr Float X0 = 0.0;
constexpr Float X1 = 1.0;
constexpr Float Y0 = 0.0;
constexpr Float Y1 = 1.0;

// = Manufactured solutions ========================================================================
// `Dirichlet` takes plain function pointers, which cannot capture, so the degree of the
// manufactured polynomial is passed through a file-scope variable.
Index g_deg = 1;

constexpr std::array<Float, MAX_BC_ORDER + 1> PC{0.7, -1.3, 0.9, 0.4, -0.25, 0.6};
constexpr std::array<Float, MAX_BC_ORDER + 1> QC{-0.4, 1.1, -0.6, 0.8, 0.35, -0.5};

auto poly(const std::array<Float, MAX_BC_ORDER + 1>& c, Float t) -> Float {
  Float v = 0.0;
  Float p = 1.0;
  for (Index k = 0; k <= g_deg; ++k) {
    v += c[static_cast<size_t>(k)] * p;
    p *= t;
  }
  return v;
}

// Separable, so that the tensor-product fill of the ghost corners is exact as well.
auto U_poly(Float x, Float y) -> Float { return poly(PC, x) * poly(QC, y); }
auto V_poly(Float x, Float y) -> Float { return poly(QC, x) * poly(PC, y); }

// Smooth but not polynomial; used for the accuracy study and the operator checks.
auto U_trig(Float x, Float y) -> Float { return std::sin(2.0 * pi * x) * std::cos(2.0 * pi * y); }
auto V_trig(Float x, Float y) -> Float { return -std::cos(2.0 * pi * x) * std::sin(2.0 * pi * y); }

// Zero normal derivative on all four walls of [0, 1] x [0, 1].
auto U_neum(Float x, Float y) -> Float { return std::cos(2.0 * pi * x) * std::cos(2.0 * pi * y); }
auto V_neum(Float x, Float y) -> Float {
  return 0.5 * std::cos(2.0 * pi * x) * std::cos(2.0 * pi * y);
}

// Boundary traces of the manufactured fields, as function pointers.
auto U_poly_left(Float y, Float /*t*/) -> Float { return U_poly(X0, y); }
auto U_poly_right(Float y, Float /*t*/) -> Float { return U_poly(X1, y); }
auto U_poly_bottom(Float x, Float /*t*/) -> Float { return U_poly(x, Y0); }
auto U_poly_top(Float x, Float /*t*/) -> Float { return U_poly(x, Y1); }
auto V_poly_left(Float y, Float /*t*/) -> Float { return V_poly(X0, y); }
auto V_poly_right(Float y, Float /*t*/) -> Float { return V_poly(X1, y); }
auto V_poly_bottom(Float x, Float /*t*/) -> Float { return V_poly(x, Y0); }
auto V_poly_top(Float x, Float /*t*/) -> Float { return V_poly(x, Y1); }

auto U_trig_left(Float y, Float /*t*/) -> Float { return U_trig(X0, y); }
auto U_trig_right(Float y, Float /*t*/) -> Float { return U_trig(X1, y); }
auto U_trig_bottom(Float x, Float /*t*/) -> Float { return U_trig(x, Y0); }
auto U_trig_top(Float x, Float /*t*/) -> Float { return U_trig(x, Y1); }
auto V_trig_left(Float y, Float /*t*/) -> Float { return V_trig(X0, y); }
auto V_trig_right(Float y, Float /*t*/) -> Float { return V_trig(X1, y); }
auto V_trig_bottom(Float x, Float /*t*/) -> Float { return V_trig(x, Y0); }
auto V_trig_top(Float x, Float /*t*/) -> Float { return V_trig(x, Y1); }

// = Test helpers ==================================================================================
struct Checker {
  Index failures = 0;
  Index reported = 0;
  Float worst    = 0.0;

  void expect(Float got, Float want, Float tol, const char* what, Index i, Index j) {
    const Float err = std::abs(got - want);
    worst           = std::max(worst, err);
    if (!(err <= tol)) {
      failures += 1;
      if (reported < 5) {
        reported += 1;
        Igor::Error("{} at ({}, {}): expected {:.17e} but got {:.17e} (err {:.3e})",
                    what,
                    i,
                    j,
                    want,
                    got,
                    err);
      }
    }
  }
};

// Position of the sample points of the two face-velocity components.
auto x_of_ux(const Grid<Float>& grid, Index i) -> Float { return grid.x(i); }
auto y_of_ux(const Grid<Float>& grid, Index j) -> Float { return grid.ym(j); }
auto x_of_uy(const Grid<Float>& grid, Index i) -> Float { return grid.xm(i); }
auto y_of_uy(const Grid<Float>& grid, Index j) -> Float { return grid.y(j); }

void init_interior(const Grid<Float>& grid,
                   FaceVector<Float> u,
                   Float (*U)(Float, Float),
                   Float (*V)(Float, Float)) {
  grid.foreach_face_i<Dimension::X>(
      FOREACH_FUNC { u.x(i, j) = U(x_of_ux(grid, i), y_of_ux(grid, j)); });
  grid.foreach_face_i<Dimension::Y>(
      FOREACH_FUNC { u.y(i, j) = V(x_of_uy(grid, i), y_of_uy(grid, j)); });
}

// Largest deviation of the whole padded array (interior included) from the analytic field.
auto max_error(const Grid<Float>& grid,
               FaceVector<Float> u,
               Float (*U)(Float, Float),
               Float (*V)(Float, Float)) -> Float {
  Float e = 0.0;
  grid.foreach_range<Exec::SERIAL>(
      -u.x.nghost(),
      u.x.nx() + u.x.nghost(),
      -u.x.nghost(),
      u.x.ny() + u.x.nghost(),
      [&](Index i, Index j) {
        e = std::max(e, std::abs(u.x(i, j) - U(x_of_ux(grid, i), y_of_ux(grid, j))));
      });
  grid.foreach_range<Exec::SERIAL>(
      -u.y.nghost(),
      u.y.nx() + u.y.nghost(),
      -u.y.nghost(),
      u.y.ny() + u.y.nghost(),
      [&](Index i, Index j) {
        e = std::max(e, std::abs(u.y(i, j) - V(x_of_uy(grid, i), y_of_uy(grid, j))));
      });
  return e;
}

auto dirichlet_poly_bconds(Index order) -> VelocityBConds<Float> {
  return {
      .left   = Dirichlet<Float>{.U = U_poly_left, .V = V_poly_left, .order = order},
      .right  = Dirichlet<Float>{.U = U_poly_right, .V = V_poly_right, .order = order},
      .bottom = Dirichlet<Float>{.U = U_poly_bottom, .V = V_poly_bottom, .order = order},
      .top    = Dirichlet<Float>{.U = U_poly_top, .V = V_poly_top, .order = order},
  };
}

auto dirichlet_trig_bconds(Index order) -> VelocityBConds<Float> {
  return {
      .left   = Dirichlet<Float>{.U = U_trig_left, .V = V_trig_left, .order = order},
      .right  = Dirichlet<Float>{.U = U_trig_right, .V = V_trig_right, .order = order},
      .bottom = Dirichlet<Float>{.U = U_trig_bottom, .V = V_trig_bottom, .order = order},
      .top    = Dirichlet<Float>{.U = U_trig_top, .V = V_trig_top, .order = order},
  };
}

// = 1. Polynomial exactness =======================================================================
// A degree-`order` polynomial must be reconstructed exactly in every ghost cell, for every number
// of ghost layers. Because the field is separable this also covers the ghost corners, which are
// filled as a tensor product of the two one-dimensional reconstructions.
auto test_polynomial_exactness() -> bool {
  Checker c;
  constexpr Index N   = 16;
  constexpr Float tol = 1e-10;

  for (Index order = 1; order <= 3; ++order) {
    for (Index nghost = 1; nghost <= 4; ++nghost) {
      g_deg = order;
      Grid<Float> grid(X0, X1, N, Y0, Y1, N, nghost);
      auto u = grid.alloc_face_vector();
      fill(u, std::numeric_limits<Float>::quiet_NaN());
      init_interior(grid, u, U_poly, V_poly);
      apply_velocity_bconds(grid, dirichlet_poly_bconds(order), u);

      grid.foreach_range<Exec::SERIAL>(
          -u.x.nghost(),
          u.x.nx() + u.x.nghost(),
          -u.x.nghost(),
          u.x.ny() + u.x.nghost(),
          [&](Index i, Index j) {
            c.expect(u.x(i, j),
                     U_poly(x_of_ux(grid, i), y_of_ux(grid, j)),
                     tol,
                     "Dirichlet u.x",
                     i,
                     j);
          });
      grid.foreach_range<Exec::SERIAL>(
          -u.y.nghost(),
          u.y.nx() + u.y.nghost(),
          -u.y.nghost(),
          u.y.ny() + u.y.nghost(),
          [&](Index i, Index j) {
            c.expect(u.y(i, j),
                     V_poly(x_of_uy(grid, i), y_of_uy(grid, j)),
                     tol,
                     "Dirichlet u.y",
                     i,
                     j);
          });

      if (c.failures > 0) {
        Igor::Error("Polynomial exactness failed for order={}, nghost={}.", order, nghost);
        return false;
      }
    }
  }
  Igor::Info("Polynomial exactness: OK (worst deviation {:.3e}).", c.worst);
  return true;
}

// = 2. Discrete boundary data =====================================================================
// The operators the solver evaluates at the wall must return the prescribed data. `order == 1` is
// the case in which they do so exactly; see `test_accuracy` for the higher-order behaviour.
//
// Where two walls meet, the two boundary conditions generally disagree; `apply_velocity_bconds`
// resolves that in favour of the bottom/top condition. The wall rows and columns that are shared
// with another wall are therefore excluded below.
auto test_dirichlet_operators() -> bool {
  Checker c;
  constexpr Index N   = 32;
  constexpr Float tol = 1e-14;

  for (Index nghost = 1; nghost <= 4; ++nghost) {
    Grid<Float> grid(X0, X1, N, Y0, Y1, N, nghost);
    auto u = grid.alloc_face_vector();
    fill(u, std::numeric_limits<Float>::quiet_NaN());
    init_interior(grid, u, U_trig, V_trig);
    apply_velocity_bconds(grid, dirichlet_trig_bconds(1), u);

    const Index nx = grid.nx();
    const Index ny = grid.ny();

    // The normal component sits on the wall and must carry the boundary value exactly.
    for (Index j = 0; j < ny; ++j) {
      c.expect(u.x(0, j), U_trig(X0, grid.ym(j)), tol, "left u.x on wall", 0, j);
      c.expect(u.x(nx, j), U_trig(X1, grid.ym(j)), tol, "right u.x on wall", nx, j);
    }
    for (Index i = -nghost; i < nx + nghost; ++i) {
      c.expect(u.y(i, 0), V_trig(grid.xm(i), Y0), tol, "bottom u.y on wall", i, 0);
      c.expect(u.y(i, ny), V_trig(grid.xm(i), Y1), tol, "top u.y on wall", i, ny);
    }

    // The tangential component is reconstructed so that the two-point average across the wall --
    // the value the viscous stress at the boundary vertices sees -- is the boundary value.
    for (Index j = 1; j < ny; ++j) {
      c.expect(0.5 * (u.y(-1, j) + u.y(0, j)),
               V_trig(X0, grid.y(j)),
               tol,
               "left u.y wall average",
               -1,
               j);
      c.expect(0.5 * (u.y(nx, j) + u.y(nx - 1, j)),
               V_trig(X1, grid.y(j)),
               tol,
               "right u.y wall average",
               nx,
               j);
    }
    for (Index i = 1; i < nx; ++i) {
      c.expect(0.5 * (u.x(i, -1) + u.x(i, 0)),
               U_trig(grid.x(i), Y0),
               tol,
               "bottom u.x wall average",
               i,
               -1);
      c.expect(0.5 * (u.x(i, ny) + u.x(i, ny - 1)),
               U_trig(grid.x(i), Y1),
               tol,
               "top u.x wall average",
               i,
               ny);
    }

    if (c.failures > 0) {
      Igor::Error("Dirichlet boundary operators failed for nghost={}.", nghost);
      return false;
    }
  }
  Igor::Info("Dirichlet boundary operators: OK (worst deviation {:.3e}).", c.worst);
  return true;
}

// The discrete normal derivative across each wall must return the prescribed outward derivative.
auto test_neumann_operators() -> bool {
  Checker c;
  constexpr Index N   = 32;
  constexpr Float tol = 1e-13;

  for (const auto [dudn, dvdn] : std::array<std::array<Float, 2>, 2>{
           {{0.0, 0.0}, {1.7, -0.9}}
  }) {
    for (Index nghost = 1; nghost <= 4; ++nghost) {
      Grid<Float> grid(X0, X1, N, Y0, Y1, N, nghost);
      auto u = grid.alloc_face_vector();
      fill(u, std::numeric_limits<Float>::quiet_NaN());
      init_interior(grid, u, U_trig, V_trig);

      const Neumann bc{.clipped = false, .dUdn = dudn, .dVdn = dvdn, .order = 1};
      apply_velocity_bconds<Float, Layout::C>(grid, {bc, bc, bc, bc}, u);

      const Index nx = grid.nx();
      const Index ny = grid.ny();
      const Float dx = grid.dx();
      const Float dy = grid.dy();

      // Left/right: u.x lies on the wall (centred difference), u.y straddles it.
      for (Index j = 0; j < ny; ++j) {
        c.expect((u.x(1, j) - u.x(-1, j)) / (2.0 * dx), -dudn, tol, "left du/dn (u.x)", -1, j);
        c.expect((u.x(nx + 1, j) - u.x(nx - 1, j)) / (2.0 * dx),
                 dudn,
                 tol,
                 "right du/dn (u.x)",
                 nx + 1,
                 j);
      }
      for (Index j = 0; j <= ny; ++j) {
        c.expect((u.y(0, j) - u.y(-1, j)) / dx, -dvdn, tol, "left dv/dn (u.y)", -1, j);
        c.expect((u.y(nx, j) - u.y(nx - 1, j)) / dx, dvdn, tol, "right dv/dn (u.y)", nx, j);
      }

      // Bottom/top: u.y lies on the wall, u.x straddles it.
      for (Index i = -nghost; i < nx + nghost; ++i) {
        c.expect((u.y(i, 1) - u.y(i, -1)) / (2.0 * dy), -dvdn, tol, "bottom dv/dn (u.y)", i, -1);
        c.expect((u.y(i, ny + 1) - u.y(i, ny - 1)) / (2.0 * dy),
                 dvdn,
                 tol,
                 "top dv/dn (u.y)",
                 i,
                 ny + 1);
      }
      for (Index i = -nghost; i < nx + 1 + nghost; ++i) {
        c.expect((u.x(i, 0) - u.x(i, -1)) / dy, -dudn, tol, "bottom du/dn (u.x)", i, -1);
        c.expect((u.x(i, ny) - u.x(i, ny - 1)) / dy, dudn, tol, "top du/dn (u.x)", i, ny);
      }

      if (c.failures > 0) {
        Igor::Error("Neumann boundary operators failed for nghost={}, dUdn={}, dVdn={}.",
                    nghost,
                    dudn,
                    dvdn);
        return false;
      }
    }
  }
  Igor::Info("Neumann boundary operators: OK (worst deviation {:.3e}).", c.worst);
  return true;
}

// = 3. Structural invariants ======================================================================
// A homogeneous Neumann condition must mirror the interior samples about the wall. For the
// component sitting on the wall that means u.x(-1-m) == u.x(1+m), *not* u.x(-1-m) == u.x(0).
auto test_neumann_is_a_mirror() -> bool {
  Checker c;
  constexpr Index N = 24;

  for (Index nghost = 1; nghost <= 4; ++nghost) {
    Grid<Float> grid(X0, X1, N, Y0, Y1, N, nghost);
    auto u = grid.alloc_face_vector();
    fill(u, std::numeric_limits<Float>::quiet_NaN());
    init_interior(grid, u, U_trig, V_trig);

    const Neumann bc{.clipped = false, .dUdn = 0.0, .dVdn = 0.0, .order = 1};
    apply_velocity_bconds<Float, Layout::C>(grid, {bc, bc, bc, bc}, u);

    const Index nx = grid.nx();
    const Index ny = grid.ny();
    for (Index m = 0; m < nghost; ++m) {
      for (Index j = 0; j < ny; ++j) {
        c.expect(u.x(-1 - m, j), u.x(1 + m, j), 0.0, "left u.x mirror", -1 - m, j);
        c.expect(u.x(nx + 1 + m, j), u.x(nx - 1 - m, j), 0.0, "right u.x mirror", nx + 1 + m, j);
      }
      for (Index j = 0; j <= ny; ++j) {
        c.expect(u.y(-1 - m, j), u.y(m, j), 0.0, "left u.y mirror", -1 - m, j);
        c.expect(u.y(nx + m, j), u.y(nx - 1 - m, j), 0.0, "right u.y mirror", nx + m, j);
      }
      for (Index i = -nghost; i < nx + nghost; ++i) {
        c.expect(u.y(i, -1 - m), u.y(i, 1 + m), 0.0, "bottom u.y mirror", i, -1 - m);
        c.expect(u.y(i, ny + 1 + m), u.y(i, ny - 1 - m), 0.0, "top u.y mirror", i, ny + 1 + m);
      }
      for (Index i = -nghost; i < nx + 1 + nghost; ++i) {
        c.expect(u.x(i, -1 - m), u.x(i, m), 0.0, "bottom u.x mirror", i, -1 - m);
        c.expect(u.x(i, ny + m), u.x(i, ny - 1 - m), 0.0, "top u.x mirror", i, ny + m);
      }
    }
    if (c.failures > 0) {
      Igor::Error("Neumann mirror property failed for nghost={}.", nghost);
      return false;
    }
  }
  Igor::Info("Neumann mirror property: OK.");
  return true;
}

// Every ghost -- corners included -- must be a bit-exact copy of the interior sample it wraps onto.
auto test_periodic_wrap() -> bool {
  Checker c;
  constexpr Index N = 20;

  const auto wrap = [](Index i, Index n, Stagger stagger) -> Index {
    const Index period = stagger == Stagger::ALIGNED ? n - 1 : n;
    Index r            = i % period;
    if (r < 0) { r += period; }
    return r;
  };

  for (Index nghost = 1; nghost <= 4; ++nghost) {
    Grid<Float> grid(X0, X1, N, Y0, Y1, N, nghost);
    auto u = grid.alloc_face_vector();
    fill(u, std::numeric_limits<Float>::quiet_NaN());
    init_interior(grid, u, U_trig, V_trig);
    apply_velocity_bconds<Float, Layout::C>(grid, {Periodic{}, Periodic{}, Periodic{}, Periodic{}}, u);

    grid.foreach_range<Exec::SERIAL>(
        -nghost,
        u.x.nx() + nghost,
        -nghost,
        u.x.ny() + nghost,
        [&](Index i, Index j) {
          if (i >= 0 && i < u.x.nx() && j >= 0 && j < u.x.ny()) { return; }
          const Index wi = i >= 0 && i < u.x.nx() ? i : wrap(i, u.x.nx(), Stagger::ALIGNED);
          const Index wj = j >= 0 && j < u.x.ny() ? j : wrap(j, u.x.ny(), Stagger::STAGGERED);
          c.expect(u.x(i, j), u.x(wi, wj), 0.0, "periodic u.x", i, j);
        });
    grid.foreach_range<Exec::SERIAL>(
        -nghost,
        u.y.nx() + nghost,
        -nghost,
        u.y.ny() + nghost,
        [&](Index i, Index j) {
          if (i >= 0 && i < u.y.nx() && j >= 0 && j < u.y.ny()) { return; }
          const Index wi = i >= 0 && i < u.y.nx() ? i : wrap(i, u.y.nx(), Stagger::STAGGERED);
          const Index wj = j >= 0 && j < u.y.ny() ? j : wrap(j, u.y.ny(), Stagger::ALIGNED);
          c.expect(u.y(i, j), u.y(wi, wj), 0.0, "periodic u.y", i, j);
        });

    if (c.failures > 0) {
      Igor::Error("Periodic wrap failed for nghost={}.", nghost);
      return false;
    }
  }
  Igor::Info("Periodic wrap: OK.");
  return true;
}

// The values in the first ghost layer must not depend on how many further layers are filled.
auto test_ghost_count_independence() -> bool {
  Checker c;
  constexpr Index N = 24;

  for (Index order = 1; order <= 3; ++order) {
    Grid<Float> grid1(X0, X1, N, Y0, Y1, N, 1);
    Grid<Float> grid4(X0, X1, N, Y0, Y1, N, 4);
    auto u1 = grid1.alloc_face_vector();
    auto u4 = grid4.alloc_face_vector();
    fill(u1, std::numeric_limits<Float>::quiet_NaN());
    fill(u4, std::numeric_limits<Float>::quiet_NaN());
    init_interior(grid1, u1, U_trig, V_trig);
    init_interior(grid4, u4, U_trig, V_trig);
    apply_velocity_bconds(grid1, dirichlet_trig_bconds(order), u1);
    apply_velocity_bconds(grid4, dirichlet_trig_bconds(order), u4);

    const Index nx = grid1.nx();
    const Index ny = grid1.ny();
    for (Index j = -1; j <= ny; ++j) {
      c.expect(u1.x(-1, j), u4.x(-1, j), 0.0, "u.x left ghost vs nghost=4", -1, j);
      c.expect(u1.x(nx + 1, j), u4.x(nx + 1, j), 0.0, "u.x right ghost vs nghost=4", nx + 1, j);
      c.expect(u1.y(-1, j), u4.y(-1, j), 0.0, "u.y left ghost vs nghost=4", -1, j);
      c.expect(u1.y(nx, j), u4.y(nx, j), 0.0, "u.y right ghost vs nghost=4", nx, j);
    }
    for (Index i = -1; i <= nx; ++i) {
      c.expect(u1.x(i, -1), u4.x(i, -1), 0.0, "u.x bottom ghost vs nghost=4", i, -1);
      c.expect(u1.x(i, ny), u4.x(i, ny), 0.0, "u.x top ghost vs nghost=4", i, ny);
      c.expect(u1.y(i, -1), u4.y(i, -1), 0.0, "u.y bottom ghost vs nghost=4", i, -1);
      c.expect(u1.y(i, ny + 1), u4.y(i, ny + 1), 0.0, "u.y top ghost vs nghost=4", i, ny + 1);
    }

    if (c.failures > 0) {
      Igor::Error("Ghost-count independence failed for order={}.", order);
      return false;
    }
  }
  Igor::Info("Ghost-count independence: OK.");
  return true;
}

// = Cell-centered scalar boundary conditions ======================================================
auto test_scalar_bconds() -> bool {
  Checker c;
  constexpr Index N   = 24;
  constexpr Float tol = 1e-14;

  for (Index nghost = 1; nghost <= 4; ++nghost) {
    Grid<Float> grid(X0, X1, N, Y0, Y1, N, nghost);
    const Index nx = grid.nx();
    const Index ny = grid.ny();

    // Homogeneous Neumann is a mirror.
    auto f = grid.alloc_scalar();
    fill(f, std::numeric_limits<Float>::quiet_NaN());
    grid.foreach_i(FOREACH_FUNC { f(i, j) = U_trig(grid.xm(i), grid.ym(j)); });
    apply_neumann_bconds(grid, f);
    for (Index m = 0; m < nghost; ++m) {
      for (Index j = 0; j < ny; ++j) {
        c.expect(f(-1 - m, j), f(m, j), 0.0, "scalar Neumann left mirror", -1 - m, j);
        c.expect(f(nx + m, j), f(nx - 1 - m, j), 0.0, "scalar Neumann right mirror", nx + m, j);
      }
      for (Index i = -nghost; i < nx + nghost; ++i) {
        c.expect(f(i, -1 - m), f(i, m), 0.0, "scalar Neumann bottom mirror", i, -1 - m);
        c.expect(f(i, ny + m), f(i, ny - 1 - m), 0.0, "scalar Neumann top mirror", i, ny + m);
      }
    }

    // Dirichlet, CONSTANT: the ghost carries the boundary value, so an upwinded inflow flux
    // transports exactly that value across the boundary face.
    constexpr Float value = 0.375;
    auto g                = grid.alloc_scalar();
    fill(g, std::numeric_limits<Float>::quiet_NaN());
    grid.foreach_i(FOREACH_FUNC { g(i, j) = U_trig(grid.xm(i), grid.ym(j)); });
    apply_dirichlet_bconds(grid, g, value, DirichletFill::CONSTANT);
    for (Index m = 0; m < nghost; ++m) {
      for (Index j = 0; j < ny; ++j) {
        c.expect(g(-1 - m, j), value, tol, "scalar Dirichlet CONSTANT left", -1 - m, j);
        c.expect(g(nx + m, j), value, tol, "scalar Dirichlet CONSTANT right", nx + m, j);
      }
    }

    // Dirichlet, EXTRAPOLATE: the interpolant across the wall hits the boundary value.
    auto s = grid.alloc_scalar();
    fill(s, std::numeric_limits<Float>::quiet_NaN());
    grid.foreach_i(FOREACH_FUNC { s(i, j) = U_trig(grid.xm(i), grid.ym(j)); });
    apply_dirichlet_bconds(grid, s, value, DirichletFill::EXTRAPOLATE);
    for (Index j = 0; j < ny; ++j) {
      c.expect(0.5 * (s(-1, j) + s(0, j)), value, tol, "scalar Dirichlet left wall", -1, j);
      c.expect(
          0.5 * (s(nx, j) + s(nx - 1, j)), value, tol, "scalar Dirichlet right wall", nx, j);
    }
    for (Index i = 0; i < nx; ++i) {
      c.expect(0.5 * (s(i, -1) + s(i, 0)), value, tol, "scalar Dirichlet bottom wall", i, -1);
      c.expect(
          0.5 * (s(i, ny) + s(i, ny - 1)), value, tol, "scalar Dirichlet top wall", i, ny);
    }

    if (c.failures > 0) {
      Igor::Error("Scalar boundary conditions failed for nghost={}.", nghost);
      return false;
    }
  }
  Igor::Info("Scalar boundary conditions: OK.");
  return true;
}

// = 4. Order of accuracy ==========================================================================
// For a smooth, non-polynomial field the ghost values must converge at rate order + 1. Only the
// finest pair of resolutions is asserted on; the coarse ones are printed but are still
// pre-asymptotic.
//
// Note how large the absolute errors are for the higher orders: the worst ghost is always a corner,
// which is reached by two successive extrapolations and therefore squares their weight growth.
auto test_accuracy() -> bool {
  constexpr size_t NRES = 4;
  constexpr std::array<Index, NRES> resolutions{16, 32, 64, 128};
  constexpr Index nghost = 3;
  bool ok                = true;

  for (Index order = 1; order <= 3; ++order) {
    std::array<Float, NRES> err_dirichlet{};
    std::array<Float, NRES> err_neumann{};

    for (size_t n = 0; n < NRES; ++n) {
      const Index N = resolutions[n];

      Grid<Float> grid(X0, X1, N, Y0, Y1, N, nghost);
      auto u = grid.alloc_face_vector();
      fill(u, std::numeric_limits<Float>::quiet_NaN());
      init_interior(grid, u, U_trig, V_trig);
      apply_velocity_bconds(grid, dirichlet_trig_bconds(order), u);
      err_dirichlet[n] = max_error(grid, u, U_trig, V_trig);

      auto v = grid.alloc_face_vector();
      fill(v, std::numeric_limits<Float>::quiet_NaN());
      init_interior(grid, v, U_neum, V_neum);
      const Neumann bc{.clipped = false,
                       .dUdn    = 0.0,
                       .dVdn    = 0.0,
                       .fill    = NeumannFill::EXTRAPOLATE,
                       .order   = order};
      apply_velocity_bconds<Float, Layout::C>(grid, {bc, bc, bc, bc}, v);
      err_neumann[n] = max_error(grid, v, U_neum, V_neum);
    }

    const auto check_rates = [&](const std::array<Float, NRES>& err, const char* what) {
      constexpr Float floor_ = 1e-14;
      const Float expected   = static_cast<Float>(order + 1);
      for (size_t n = 1; n < NRES; ++n) {
        const Float rate = std::log2(std::max(err[n - 1], floor_) / std::max(err[n], floor_));
        Igor::Info("{} order={}: N={:4} err={:.3e} rate={:.2f}",
                   what,
                   order,
                   resolutions[n],
                   err[n],
                   rate);
      }
      const Float rate =
          std::log2(std::max(err[NRES - 2], floor_) / std::max(err[NRES - 1], floor_));
      // Errors already at round-off level cannot show a rate; accept those.
      if (err[NRES - 1] > 1e-12 && rate < expected - 0.25) {
        Igor::Error("{} order={}: expected rate >= {:.2f} between N={} and N={}, but got {:.2f}",
                    what,
                    order,
                    expected - 0.25,
                    resolutions[NRES - 2],
                    resolutions[NRES - 1],
                    rate);
        ok = false;
      }
    };
    check_rates(err_dirichlet, "Dirichlet");
    check_rates(err_neumann, "Neumann EXTRAPOLATE");
  }

  // A field that is even about all four walls is reproduced *exactly* by the reflection, corners
  // included -- reflection is not limited to a fixed polynomial degree.
  for (Index N : resolutions) {
    Grid<Float> grid(X0, X1, N, Y0, Y1, N, nghost);
    auto v = grid.alloc_face_vector();
    fill(v, std::numeric_limits<Float>::quiet_NaN());
    init_interior(grid, v, U_neum, V_neum);
    const Neumann bc{.clipped = false, .dUdn = 0.0, .dVdn = 0.0, .fill = NeumannFill::MIRROR};
    apply_velocity_bconds<Float, Layout::C>(grid, {bc, bc, bc, bc}, v);
    const Float err = max_error(grid, v, U_neum, V_neum);
    Igor::Info("Neumann MIRROR: N={:4} err={:.3e}", N, err);
    if (err > 1e-13) {
      Igor::Error("Neumann MIRROR is not exact for an even field: N={}, err={:.3e}", N, err);
      ok = false;
    }
  }
  return ok;
}

// =================================================================================================
auto main() -> int {
  bool any_failed = false;

  const auto run = [&any_failed](bool passed, const char* name) {
    if (!passed) {
      Igor::Error("{} failed.", name);
      any_failed = true;
    }
  };

  run(test_polynomial_exactness(), "Polynomial exactness");
  run(test_dirichlet_operators(), "Dirichlet boundary operators");
  run(test_neumann_operators(), "Neumann boundary operators");
  run(test_neumann_is_a_mirror(), "Neumann mirror property");
  run(test_periodic_wrap(), "Periodic wrap");
  run(test_ghost_count_independence(), "Ghost-count independence");
  run(test_scalar_bconds(), "Scalar boundary conditions");
  run(test_accuracy(), "Order of accuracy");

  return any_failed ? 1 : 0;
}
