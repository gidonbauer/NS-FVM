#ifndef NS_FVM_GRID_HPP_
#define NS_FVM_GRID_HPP_

#include <cstdint>
#include <cstdlib>
#include <type_traits>

#include <Igor/Logging.hpp>

#ifndef NS_FVM_INDEX_TYPE
using Index = int32_t;
#else
static_assert(std::is_integral_v<NS_FVM_INDEX_TYPE> && std::is_signed_v<NS_FVM_INDEX_TYPE>,
              "NS_FVM_INDEX_TYPE must be a signed integer type.");
using Index = NS_FVM_INDEX_TYPE;
#endif  // FS_INDEX_TYPE

enum class Layout { C, F };

template <typename Float, Layout LAYOUT>
requires(std::is_trivially_constructible_v<Float> && std::is_trivially_destructible_v<Float>)
class Scalar;

template <typename Float, Layout LAYOUT>
requires(std::is_trivially_constructible_v<Float> && std::is_trivially_destructible_v<Float>)
class Vector;

template <typename Float, Layout LAYOUT>
requires(std::is_trivially_constructible_v<Float> && std::is_trivially_destructible_v<Float>)
class FaceVector;

template <typename Float, Layout LAYOUT = Layout::C>
class Grid {
  using Scalar     = Scalar<Float, LAYOUT>;
  using Vector     = Vector<Float, LAYOUT>;
  using FaceVector = FaceVector<Float, LAYOUT>;

  Float m_x_min;
  Float m_x_max;
  Float m_dx;
  Index m_nx;

  Float m_y_min;
  Float m_y_max;
  Float m_dy;
  Index m_ny;

  Float m_dv;

  Index m_nghost;

  constexpr auto alloc(Index nx, Index ny, Index nghost) noexcept -> Scalar {
    IGOR_ASSERT(nx > 0 && ny > 0 && nghost >= 0,
                "Invalid dimensions: nx={}, ny={}, nghost={}",
                nx,
                ny,
                nghost);
    const auto size = static_cast<size_t>(nx + 2 * nghost) * static_cast<size_t>(ny + 2 * nghost);
    Float* data     = static_cast<Float*>(std::calloc(size, sizeof(Float)));  // NOLINT
    IGOR_ASSERT(data != nullptr, "Could not allocate scalar.");
    return {data, nx, ny, nghost};
  }

 public:
  constexpr Grid(Float x_min,
                 Float x_max,
                 Index nx,
                 Float y_min,
                 Float y_max,
                 Index ny,
                 Index nghost = 1) noexcept
      : m_x_min(x_min),
        m_x_max(x_max),
        m_dx((x_max - x_min) / nx),
        m_nx(nx),
        m_y_min(y_min),
        m_y_max(y_max),
        m_dy((y_max - y_min) / ny),
        m_ny(ny),
        m_dv(m_dx * m_dy),
        m_nghost(nghost) {}

  [[nodiscard]] constexpr auto x_min() const noexcept -> Float { return m_x_min; }
  [[nodiscard]] constexpr auto x_max() const noexcept -> Float { return m_x_max; }
  [[nodiscard]] constexpr auto dx() const noexcept -> Float { return m_dx; }
  [[nodiscard]] constexpr auto nx() const noexcept -> Index { return m_nx; }
  [[nodiscard]] constexpr auto y_min() const noexcept -> Float { return m_y_min; }
  [[nodiscard]] constexpr auto y_max() const noexcept -> Float { return m_y_max; }
  [[nodiscard]] constexpr auto dy() const noexcept -> Float { return m_dy; }
  [[nodiscard]] constexpr auto ny() const noexcept -> Index { return m_ny; }
  [[nodiscard]] constexpr auto dv() const noexcept -> Float { return m_dv; }
  [[nodiscard]] constexpr auto nghost() const noexcept -> Index { return m_nghost; }

  [[nodiscard]] constexpr auto alloc_scalar() noexcept -> Scalar {
    return alloc(m_nx, m_ny, m_nghost);
  }
  constexpr void free_scalar(Scalar& s) noexcept { std::free(s.data()); /* NOLINT */ }

  [[nodiscard]] constexpr auto alloc_vector() noexcept -> Vector {
    auto x = alloc_scalar();
    auto y = alloc_scalar();
    return {x, y};
  }
  constexpr void free_vector(Vector& v) noexcept {
    free_scalar(v.x);
    free_scalar(v.y);
  }

  [[nodiscard]] constexpr auto alloc_face_vector() noexcept -> FaceVector {
    auto x = alloc(m_nx + 1, m_ny, m_nghost);
    auto y = alloc(m_nx, m_ny + 1, m_nghost);
    return {x, y};
  }
  constexpr void free_face_vector(FaceVector& v) noexcept {
    free_scalar(v.x);
    free_scalar(v.y);
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

  constexpr auto operator()(Index i, Index j) noexcept -> Float& {
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

  constexpr auto operator()(Index i, Index j) const noexcept -> const Float& {
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

  [[nodiscard]] constexpr auto data() noexcept -> Float* { return m_data; }
  [[nodiscard]] constexpr auto data() const noexcept -> const Float* { return m_data; }

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

  constexpr auto left(Index i, Index j) noexcept -> Float& { return x(i, j); }
  constexpr auto left(Index i, Index j) const noexcept -> const Float& { return x(i, j); }
  constexpr auto right(Index i, Index j) noexcept -> Float& { return x(i + 1, j); }
  constexpr auto right(Index i, Index j) const noexcept -> const Float& { return x(i + 1, j); }

  constexpr auto bottom(Index i, Index j) noexcept -> Float& { return y(i, j); }
  constexpr auto bottom(Index i, Index j) const noexcept -> const Float& { return y(i, j); }
  constexpr auto top(Index i, Index j) noexcept -> Float& { return y(i, j + 1); }
  constexpr auto top(Index i, Index j) const noexcept -> const Float& { return y(i, j + 1); }

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
constexpr void copy(const FaceVector<Float, LAYOUT> src, Vector<Float, LAYOUT> dst) {
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
