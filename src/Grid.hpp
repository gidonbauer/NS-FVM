#ifndef NS_FVM_GRID_HPP_
#define NS_FVM_GRID_HPP_

#include <atomic>
#include <cstdint>
#include <cstdlib>
#include <memory>
#include <type_traits>

#ifdef NS_FVM_PARALLEL
#include <Kokkos_Core.hpp>
#endif  // NS_FVM_PARALLEL

#include <Igor/Logging.hpp>

#if defined(__clang__) || defined(__GNUC__)
#define NS_FVM_FOREACH_DEF __attribute__((flatten)) __attribute__((always_inline))
#else
#define NS_FVM_FOREACH_DEF
#warning "NS_FVM_FOREACH_DEF is not defined for this compiler; foreach kernels may not vectorize."
#endif

#define FOREACH_FUNC [=](Index i, Index j)

#ifndef NS_FVM_INDEX_TYPE
using Index = int32_t;
#else
static_assert(std::is_integral_v<NS_FVM_INDEX_TYPE> && std::is_signed_v<NS_FVM_INDEX_TYPE>,
              "NS_FVM_INDEX_TYPE must be a signed integer type.");
using Index = NS_FVM_INDEX_TYPE;
#endif  // FS_INDEX_TYPE

template <typename T>
void update_maximum_atomic(std::atomic<T>& maximum_value, T const& value) noexcept {
  T prev_value = maximum_value;
  while (prev_value < value && !maximum_value.compare_exchange_weak(prev_value, value)) {}
}

// TODO: Consider blocked layout. Each block should have its own ghost layer around it. This
//       requires a halo exchange after each update but opens up the path to efficient
//       parallelization.
enum class Layout { C, F };
enum class Dimension { X, Y };
enum class Coordinates { CARTESIAN, POLAR };
enum class Exec { PARALLEL, SERIAL };

template <typename Float, Layout LAYOUT>
requires(std::is_trivially_constructible_v<Float> && std::is_trivially_destructible_v<Float>)
class Scalar;

template <typename Float, Layout LAYOUT>
requires(std::is_trivially_constructible_v<Float> && std::is_trivially_destructible_v<Float>)
class VertexScalar;

template <typename Float, Layout LAYOUT>
requires(std::is_trivially_constructible_v<Float> && std::is_trivially_destructible_v<Float>)
class Vector;

template <typename Float, Layout LAYOUT>
requires(std::is_trivially_constructible_v<Float> && std::is_trivially_destructible_v<Float>)
class FaceVector;

#ifdef NS_FVM_PARALLEL
class KokkosRefCount {
  static size_t m_count;

 public:
  constexpr KokkosRefCount() noexcept {
    if (m_count == 0) { Kokkos::initialize(); }
    m_count += 1;
  }
  constexpr KokkosRefCount(const KokkosRefCount& /*other*/) noexcept {
    IGOR_ASSERT(m_count >= 1, "There must be at least one Kokkos user.");
    m_count += 1;
  }
  constexpr KokkosRefCount(KokkosRefCount&& /*other*/) noexcept {
    IGOR_ASSERT(m_count >= 1, "There must be at least one Kokkos user.");
    m_count += 1;
  }
  constexpr auto operator=(const KokkosRefCount& other) noexcept -> KokkosRefCount& = default;
  constexpr auto operator=(KokkosRefCount&& /*other*/) noexcept -> KokkosRefCount&  = default;
  constexpr ~KokkosRefCount() noexcept {
    IGOR_ASSERT(m_count >= 1, "There must be at least one Kokkos user.");
    if (m_count == 1) { Kokkos::finalize(); }
    m_count -= 1;
  }

  static constexpr auto count() noexcept -> size_t { return m_count; }
};
size_t KokkosRefCount::m_count = 0;
#endif  // NS_FVM_PARALLEL

template <typename Float, Layout LAYOUT = Layout::C>
class Grid {
  using Scalar       = Scalar<Float, LAYOUT>;
  using VertexScalar = VertexScalar<Float, LAYOUT>;
  using Vector       = Vector<Float, LAYOUT>;
  using FaceVector   = FaceVector<Float, LAYOUT>;

  Float m_x_min;
  Float m_x_max;
  Float m_dx;
  Index m_nx;

  Float m_y_min;
  Float m_y_max;
  Float m_dy;
  Index m_ny;

  Index m_nghost;

  Coordinates m_coords;

  std::shared_ptr<std::vector<Float*>> m_to_free = std::make_shared<std::vector<Float*>>();
#ifdef NS_FVM_PARALLEL
  KokkosRefCount m_kokkos_count;
#endif  // NS_FVM_PARALLEL

  constexpr auto alloc(Index nx, Index ny, Index nghost) const noexcept -> Scalar {
    IGOR_ASSERT(nx > 0 && ny > 0 && nghost >= 0,
                "Invalid dimensions: nx={}, ny={}, nghost={}",
                nx,
                ny,
                nghost);
    const auto size = static_cast<size_t>(nx + 2 * nghost) * static_cast<size_t>(ny + 2 * nghost);
    Float* data     = static_cast<Float*>(std::calloc(size, sizeof(Float)));  // NOLINT
    IGOR_ASSERT(data != nullptr, "Could not allocate scalar.");
    m_to_free->push_back(data);
    return {data, nx, ny, nghost};
  }

 public:
  constexpr Grid(Float x_min,
                 Float x_max,
                 Index nx,
                 Float y_min,
                 Float y_max,
                 Index ny,
                 Index nghost       = 1,
                 Coordinates coords = Coordinates::CARTESIAN) noexcept
      : m_x_min(x_min),
        m_x_max(x_max),
        m_dx((x_max - x_min) / nx),
        m_nx(nx),
        m_y_min(y_min),
        m_y_max(y_max),
        m_dy((y_max - y_min) / ny),
        m_ny(ny),
        m_nghost(nghost),
        m_coords(coords) {}

  constexpr Grid(const Grid& other) noexcept                    = default;
  constexpr Grid(Grid&& other) noexcept                         = default;
  constexpr auto operator=(const Grid& other) noexcept -> Grid& = default;
  constexpr auto operator=(Grid&& other) noexcept -> Grid&      = default;
  constexpr ~Grid() noexcept {
    if (m_to_free.use_count() == 1) {
      for (Float* ptr : *m_to_free) {
        std::free(ptr);  // NOLINT
      }
    }
  }

  [[nodiscard]] constexpr auto x_min() const noexcept -> Float { return m_x_min; }
  [[nodiscard]] constexpr auto x_max() const noexcept -> Float { return m_x_max; }
  [[nodiscard]] constexpr auto dx() const noexcept -> Float { return m_dx; }
  [[nodiscard]] constexpr auto nx() const noexcept -> Index { return m_nx; }
  [[nodiscard]] constexpr auto y_min() const noexcept -> Float { return m_y_min; }
  [[nodiscard]] constexpr auto y_max() const noexcept -> Float { return m_y_max; }
  [[nodiscard]] constexpr auto dy() const noexcept -> Float { return m_dy; }
  [[nodiscard]] constexpr auto ny() const noexcept -> Index { return m_ny; }
  [[nodiscard]] constexpr auto dv(Index /*i*/, Index j) const noexcept -> Float {
    switch (m_coords) {
      case Coordinates::CARTESIAN: return m_dx * m_dy;
      case Coordinates::POLAR:     return ym(j) * m_dy * m_dx;
    }
    Igor::Panic("Unreachable");
  }
  [[nodiscard]] constexpr auto nghost() const noexcept -> Index { return m_nghost; }
  [[nodiscard]] constexpr auto coords() const noexcept -> Coordinates { return m_coords; }

  // TODO: Maybe pass x and y as function argument to FOREACH_FUNC
  [[nodiscard]] constexpr auto x(Index i) const noexcept -> Float { return m_x_min + i * m_dx; }
  [[nodiscard]] constexpr auto y(Index j) const noexcept -> Float { return m_y_min + j * m_dy; }
  [[nodiscard]] constexpr auto xm(Index i) const noexcept -> Float {
    return m_x_min + (i + 0.5) * m_dx;
  }
  [[nodiscard]] constexpr auto ym(Index j) const noexcept -> Float {
    return m_y_min + (j + 0.5) * m_dy;
  }

  [[nodiscard]] constexpr auto alloc_scalar() const noexcept -> Scalar {
    return alloc(m_nx, m_ny, m_nghost);
  }
  [[nodiscard]] constexpr auto alloc_vertex_scalar() const noexcept -> VertexScalar {
    return alloc(m_nx + 1, m_ny + 1, m_nghost);
  }
  [[nodiscard]] constexpr auto alloc_vector() const noexcept -> Vector {
    auto x = alloc_scalar();
    auto y = alloc_scalar();
    return {x, y};
  }
  [[nodiscard]] constexpr auto alloc_face_vector() const noexcept -> FaceVector {
    auto x = alloc(m_nx + 1, m_ny, m_nghost);
    auto y = alloc(m_nx, m_ny + 1, m_nghost);
    return {x, y};
  }

  // Iterate the logical rectangle [ilo, ihi) x [jlo, jhi), innermost over the contiguous dimension.
  template <Exec EXEC = Exec::PARALLEL, typename FUNC>
  NS_FVM_FOREACH_DEF constexpr void
  foreach_range(Index ilo, Index ihi, Index jlo, Index jhi, const FUNC& func) const noexcept {
#ifdef NS_FVM_PARALLEL
    if constexpr (EXEC == Exec::PARALLEL) {
      Kokkos::parallel_for("foreach_range", Kokkos::MDRangePolicy({ilo, jlo}, {ihi, jhi}), func);
    } else
#endif  // NS_FVM_PARALLEL
      if constexpr (LAYOUT == Layout::F) {
        // Column-major: i is contiguous.
        for (Index j = jlo; j < jhi; ++j) {
          for (Index i = ilo; i < ihi; ++i) {
            func(i, j);
          }
        }
      } else {
        // Row-major: j is contiguous.
        for (Index i = ilo; i < ihi; ++i) {
          for (Index j = jlo; j < jhi; ++j) {
            func(i, j);
          }
        }
      }
  }

  template <Dimension DIM, Exec EXEC = Exec::PARALLEL, typename FUNC>
  NS_FVM_FOREACH_DEF constexpr void foreach_face_i(const FUNC& func) const noexcept {
    const Index ihi = (DIM == Dimension::X) ? nx() + 1 : nx();
    const Index jhi = (DIM == Dimension::X) ? ny() : ny() + 1;
    foreach_range<EXEC>(0, ihi, 0, jhi, func);
  }

  template <Dimension DIM, Exec EXEC = Exec::PARALLEL, typename FUNC>
  NS_FVM_FOREACH_DEF constexpr void foreach_face_a(const FUNC& func) const noexcept {
    const Index ihi = (DIM == Dimension::X) ? nx() + nghost() + 1 : nx() + nghost();
    const Index jhi = (DIM == Dimension::X) ? ny() + nghost() : ny() + nghost() + 1;
    foreach_range<EXEC>(-nghost(), ihi, -nghost(), jhi, func);
  }

  template <Exec EXEC = Exec::PARALLEL, typename FUNC>
  NS_FVM_FOREACH_DEF constexpr void foreach_i(const FUNC& func) const noexcept {
    foreach_range<EXEC>(0, nx(), 0, ny(), func);
  }

  template <Exec EXEC = Exec::PARALLEL, typename FUNC>
  NS_FVM_FOREACH_DEF constexpr void foreach_a(const FUNC& func) const noexcept {
    foreach_range<EXEC>(-nghost(), nx() + nghost(), -nghost(), ny() + nghost(), func);
  }

  template <Exec EXEC = Exec::PARALLEL, typename FUNC>
  NS_FVM_FOREACH_DEF constexpr void foreach_vertex_i(const FUNC& func) const noexcept {
    foreach_range<EXEC>(0, nx() + 1, 0, ny() + 1, func);
  }

  template <Exec EXEC = Exec::PARALLEL, typename FUNC>
  NS_FVM_FOREACH_DEF constexpr void foreach_vertex_a(const FUNC& func) const noexcept {
    foreach_range<EXEC>(-nghost(), nx() + 1 + nghost(), -nghost(), ny() + 1 + nghost(), func);
  }
};

template <typename Float, Layout LAYOUT = Layout::C>
requires(std::is_trivially_constructible_v<Float> && std::is_trivially_destructible_v<Float>)
class Scalar {
  Float* m_data;
  Index m_nx;
  Index m_ny;
  Index m_nghost;

  [[nodiscard]] constexpr auto get_idx(Index i, Index j) const noexcept -> Index {
    if constexpr (LAYOUT == Layout::C) {
      return (j + m_nghost) + (i + m_nghost) * (m_ny + 2 * m_nghost);
    } else {
      return (i + m_nghost) + (j + m_nghost) * (m_nx + 2 * m_nghost);
    }
  }

  constexpr Scalar(Float* data, Index nx, Index ny, Index nghost) noexcept
      : m_data(data),
        m_nx(nx),
        m_ny(ny),
        m_nghost(nghost) {}

 public:
  constexpr Scalar(const Scalar& other) noexcept                    = default;
  constexpr Scalar(Scalar&& other) noexcept                         = default;
  constexpr auto operator=(const Scalar& other) noexcept -> Scalar& = default;
  constexpr auto operator=(Scalar&& other) noexcept -> Scalar&      = default;
  constexpr ~Scalar() noexcept                                      = default;

  constexpr auto operator()(Index i, Index j) const noexcept -> Float& {
    IGOR_ASSERT(i >= -m_nghost && i < m_nx + m_nghost && j >= -m_nghost && j < m_ny + m_nghost,
                "Index ({}, {}) is out of bounds for Scalar of size {}:{}x{}:{}",
                i,
                j,
                -m_nghost,
                m_nx + m_nghost,
                -m_nghost,
                m_ny + m_nghost);
    return *(data() + get_idx(i, j));
  }

  [[nodiscard]] constexpr auto data() const noexcept -> Float* { return m_data; }

  [[nodiscard]] constexpr auto size() const noexcept -> Index {
    return (m_nx + 2 * m_nghost) * (m_ny + 2 * m_nghost);
  }
  [[nodiscard]] constexpr auto nx() const noexcept -> Index { return m_nx; }
  [[nodiscard]] constexpr auto ny() const noexcept -> Index { return m_ny; }
  [[nodiscard]] constexpr auto nghost() const noexcept -> Index { return m_nghost; }

  friend class Grid<Float, LAYOUT>;
};

template <typename Float, Layout LAYOUT = Layout::C>
requires(std::is_trivially_constructible_v<Float> && std::is_trivially_destructible_v<Float>)
class VertexScalar {
  using Scalar = Scalar<Float, LAYOUT>;
  Scalar m_s;

  constexpr VertexScalar(Scalar s) noexcept
      : m_s(s) {}

 public:
  constexpr VertexScalar(const VertexScalar& other) noexcept                    = default;
  constexpr VertexScalar(VertexScalar&& other) noexcept                         = default;
  constexpr auto operator=(const VertexScalar& other) noexcept -> VertexScalar& = default;
  constexpr auto operator=(VertexScalar&& other) noexcept -> VertexScalar&      = default;
  constexpr ~VertexScalar() noexcept                                            = default;

  constexpr auto operator()(Index i, Index j) const noexcept -> Float& { return m_s(i, j); }
  [[nodiscard]] constexpr auto data() const noexcept -> Float* { return m_s.data(); }
  [[nodiscard]] constexpr auto size() const noexcept -> Index { return m_s.size(); }
  [[nodiscard]] constexpr auto nx() const noexcept -> Index { return m_s.nx(); }
  [[nodiscard]] constexpr auto ny() const noexcept -> Index { return m_s.ny(); }
  [[nodiscard]] constexpr auto nghost() const noexcept -> Index { return m_s.nghost(); }

  friend class Grid<Float, LAYOUT>;
};

template <typename Float, Layout LAYOUT = Layout::C>
requires(std::is_trivially_constructible_v<Float> && std::is_trivially_destructible_v<Float>)
class Vector {
 public:
  using Scalar = Scalar<Float, LAYOUT>;
  Scalar x;
  Scalar y;

 private:
  constexpr Vector(Scalar x, Scalar y) noexcept
      : x(x),
        y(y) {}

 public:
  constexpr Vector(const Vector& other) noexcept                    = default;
  constexpr Vector(Vector&& other) noexcept                         = default;
  constexpr auto operator=(const Vector& other) noexcept -> Vector& = default;
  constexpr auto operator=(Vector&& other) noexcept -> Vector&      = default;
  constexpr ~Vector() noexcept                                      = default;

  friend class Grid<Float, LAYOUT>;
};

template <typename Float, Layout LAYOUT = Layout::C>
requires(std::is_trivially_constructible_v<Float> && std::is_trivially_destructible_v<Float>)
class FaceVector {
 public:
  using Scalar = Scalar<Float, LAYOUT>;
  Scalar x;
  Scalar y;

 private:
  constexpr FaceVector(Scalar x, Scalar y) noexcept
      : x(x),
        y(y) {}

 public:
  constexpr FaceVector(const FaceVector& other) noexcept                    = default;
  constexpr FaceVector(FaceVector&& other) noexcept                         = default;
  constexpr auto operator=(const FaceVector& other) noexcept -> FaceVector& = default;
  constexpr auto operator=(FaceVector&& other) noexcept -> FaceVector&      = default;
  constexpr ~FaceVector() noexcept                                          = default;

  constexpr auto left(Index i, Index j) const noexcept -> Float& { return x(i, j); }
  constexpr auto right(Index i, Index j) const noexcept -> Float& { return x(i + 1, j); }
  constexpr auto bottom(Index i, Index j) const noexcept -> Float& { return y(i, j); }
  constexpr auto top(Index i, Index j) const noexcept -> Float& { return y(i, j + 1); }

  friend class Grid<Float, LAYOUT>;
};

// = Copy ==========================================================================================
template <typename Float, Layout LAYOUT>
constexpr void copy(const Scalar<Float, LAYOUT> src, Scalar<Float, LAYOUT> dst) {
  std::copy_n(src.data(), src.size(), dst.data());
}

template <typename Float, Layout LAYOUT>
constexpr void copy(const Vector<Float, LAYOUT> src, Vector<Float, LAYOUT> dst) {
  copy(src.x, dst.x);
  copy(src.y, dst.y);
}

template <typename Float, Layout LAYOUT>
constexpr void copy(const FaceVector<Float, LAYOUT> src, FaceVector<Float, LAYOUT> dst) {
  copy(src.x, dst.x);
  copy(src.y, dst.y);
}

// = Fill ==========================================================================================
template <typename Float, Layout LAYOUT>
constexpr void fill(Scalar<Float, LAYOUT> s, Float value) {
  std::fill_n(s.data(), s.size(), value);
}

template <typename Float, Layout LAYOUT>
constexpr void fill(Vector<Float, LAYOUT> v, Float value) {
  fill(v.x, value);
  fill(v.y, value);
}

template <typename Float, Layout LAYOUT>
constexpr void fill(FaceVector<Float, LAYOUT> v, Float value) {
  fill(v.x, value);
  fill(v.y, value);
}

#endif  // NS_FVM_GRID_HPP_
