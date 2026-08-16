#ifndef NS_FVM_VTK_WRITER_HPP_
#define NS_FVM_VTK_WRITER_HPP_

#include <array>
#include <bit>
#include <cstring>
#include <fstream>
#include <type_traits>
#include <vector>

#include "Grid.hpp"

template <typename Float, Layout LAYOUT>
class VTKWriter {
  std::string m_output_dir;

  Float m_x_min;
  Float m_x_max;
  Float m_dx;
  Index m_nx;

  Float m_y_min;
  Float m_y_max;
  Float m_dy;
  Index m_ny;

  std::vector<std::string> m_scalar_names;
  std::vector<Scalar<Float, LAYOUT>> m_scalar_values{};

  std::vector<std::string> m_vector_names;
  std::vector<Vector<Float, LAYOUT>> m_vector_values{};

  template <typename T>
  requires(std::is_fundamental_v<T> && (sizeof(T) == 4 || sizeof(T) == 8))
  [[nodiscard]] constexpr auto interpret_as_big_endian_bytes(T value)
      -> std::array<const char, sizeof(T)> {
    if constexpr (std::endian::native == std::endian::big) {
      return std::bit_cast<std::array<const char, sizeof(T)>>(value);
    }
    using U = std::conditional_t<sizeof(T) == 4, std::uint32_t, std::uint64_t>;
    static_assert(sizeof(T) == sizeof(U));
    return std::bit_cast<std::array<const char, sizeof(T)>>(std::byteswap(std::bit_cast<U>(value)));
  }

  // -----------------------------------------------------------------------------------------------
  void write_header(std::ofstream& out, Float t) {
    // = Write VTK header ====
    out << "# vtk DataFile Version 2.0\n";
    out << "State of FluidSolver at time t=" << t << '\n';
    out << "BINARY\n";

    // = Write grid ==========
    out << "DATASET STRUCTURED_GRID\n";
    out << "DIMENSIONS " << m_nx + 1 << ' ' << m_ny + 1 << " 1\n";
    out << "POINTS " << (m_nx + 1) * (m_ny + 1) << " double\n";
    for (Index j = 0; j < m_ny + 1; ++j) {
      for (Index i = 0; i < m_nx + 1; ++i) {
        const double xi = m_x_min + i * m_dx;
        const double yj = m_y_min + j * m_dy;
        const double zk = 0.0;
        out.write(interpret_as_big_endian_bytes(xi).data(), sizeof(xi));
        out.write(interpret_as_big_endian_bytes(yj).data(), sizeof(yj));
        out.write(interpret_as_big_endian_bytes(zk).data(), sizeof(zk));
      }
    }
    out << "\n\n";

    // = Write cell data =====
    out << "CELL_DATA " << m_nx * m_ny << '\n';
  }

  // -----------------------------------------------------------------------------------------------
  void write_scalar(std::ofstream& out, const Scalar<Float, LAYOUT> s, const std::string& name) {
    out << "SCALARS " << name << " double 1\n";
    out << "LOOKUP_TABLE default\n";
    for (Index j = 0; j < s.ny(); ++j) {
      for (Index i = 0; i < s.nx(); ++i) {
        out.write(interpret_as_big_endian_bytes(s(i, j)).data(), sizeof(s(i, j)));
      }
    }
    out << "\n\n";
  }

  // -----------------------------------------------------------------------------------------------
  void write_vector(std::ofstream& out, const Vector<Float, LAYOUT> v, const std::string& name) {
    out << "VECTORS " << name << " double\n";
    for (Index j = 0; j < v.x.ny(); ++j) {
      for (Index i = 0; i < v.x.nx(); ++i) {
        constexpr double z_comp_ij = 0.0;
        out.write(interpret_as_big_endian_bytes(v.x(i, j)).data(), sizeof(v.x(i, j)));
        out.write(interpret_as_big_endian_bytes(v.y(i, j)).data(), sizeof(v.y(i, j)));
        out.write(interpret_as_big_endian_bytes(z_comp_ij).data(), sizeof(z_comp_ij));
      }
    }
    out << "\n\n";
  }

 public:
  constexpr VTKWriter(std::string output_dir, const Grid<Float, LAYOUT>& grid)
      : m_output_dir(std::move(output_dir)),
        m_x_min(grid.x_min()),
        m_x_max(grid.x_max()),
        m_dx(grid.dx()),
        m_nx(grid.nx()),
        m_y_min(grid.y_min()),
        m_y_max(grid.y_max()),
        m_dy(grid.dy()),
        m_ny(grid.ny()) {}

  constexpr VTKWriter(const VTKWriter& other) noexcept                    = delete;
  constexpr VTKWriter(VTKWriter&& other) noexcept                         = delete;
  constexpr auto operator=(const VTKWriter& other) noexcept -> VTKWriter& = delete;
  constexpr auto operator=(VTKWriter&& other) noexcept -> VTKWriter&      = delete;
  constexpr ~VTKWriter() noexcept                                         = default;

  constexpr void add_field(std::string name, const Scalar<Float, LAYOUT> s) {
    m_scalar_names.emplace_back(std::move(name));
    m_scalar_values.push_back(s);
  }

  constexpr void add_field(std::string name, const Vector<Float, LAYOUT> v) {
    m_vector_names.emplace_back(std::move(name));
    m_vector_values.push_back(v);
  }

  constexpr auto write(Float t = -1.0) -> bool {
    static Index write_counter = 0;
    const auto filename =
        Igor::detail::format("{}/state_{:06d}.vtk", m_output_dir, write_counter++);
    std::ofstream out(filename);
    if (!out) {
      Igor::Error("Could not open file `{}`: {}", filename, std::strerror(errno));
      return false;
    }

    write_header(out, t);
    if (!out) {
      Igor::Error("Could not write header to `{}`: {}", filename, std::strerror(errno));
      return false;
    }

    for (size_t i = 0; i < m_scalar_names.size(); ++i) {
      write_scalar(out, m_scalar_values[i], m_scalar_names[i]);
    }

    for (size_t i = 0; i < m_vector_names.size(); ++i) {
      write_vector(out, m_vector_values[i], m_vector_names[i]);
    }

    return out.good();
  }
};

#endif  // NS_FVM_VTK_WRITER_HPP_
