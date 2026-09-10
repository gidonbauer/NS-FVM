#pragma once

#include <array>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <limits>
#include <string>
#include <utility>
#include <vector>

#include <H5Cpp.h>

#include <Igor/Logging.hpp>

#include "Grid.hpp"

namespace detail {

struct Indent {
  size_t n;
};

constexpr auto operator<<(std::ostream& out, const Indent& indent) -> std::ostream& {
  for (size_t i = 0; i < indent.n; ++i) {
    out << ' ';
  }
  return out;
}

}  // namespace detail

template <typename Float, Layout LAYOUT>
class HDFWriter {
  constexpr static const char* HDF_FILENAME  = "solution.h5";
  constexpr static const char* XDMF_FILENAME = "solution.xdmf2";
  static constexpr const char* END_TAGS      = "    </Grid>\n  </Domain>\n</Xdmf>\n";

  H5::H5File m_hdf_file;
  std::ofstream m_xdmf_file;

  Float m_x_min;
  Float m_dx;
  Index m_nx;

  Float m_y_min;
  Float m_dy;
  Index m_ny;

  int write_counter = 0;
  std::array<char, 7> group_name_buffer{};
  std::ostream::pos_type next_write_pos = 0;

  Coordinates m_coords;

  std::vector<std::string> m_scalar_names;
  std::vector<Scalar<Float, LAYOUT>> m_scalar_values{};

  std::vector<std::string> m_vector_names;
  std::vector<Vector<Float, LAYOUT>> m_vector_values{};

  // ===============================================================================================
  static constexpr auto pred_type() noexcept -> H5::PredType {
    if constexpr (std::is_same_v<Float, double>) {
      return H5::PredType::NATIVE_DOUBLE;
    } else if constexpr (std::is_same_v<Float, float>) {
      return H5::PredType::NATIVE_FLOAT;
    } else {
      Igor::Panic("Invalid type.");
      std::unreachable();
    }
  }

  static constexpr auto precision() noexcept -> size_t { return sizeof(Float); }

  // ===============================================================================================
  constexpr void write_grid_unmaterialized(const H5::Group& group) {
    static_assert(LAYOUT == Layout::C, "Only Layout::C is supported.");
    const auto nxp = static_cast<hsize_t>(m_nx) + 1;
    const auto nyp = static_cast<hsize_t>(m_ny) + 1;

    const std::array<hsize_t, 2> file_size{nxp, nyp};
    H5::DataSpace file_space(2, file_size.data());
    H5::DataSet dataset_x = group.createDataSet("x", pred_type(), file_space);
    H5::DataSet dataset_y = group.createDataSet("y", pred_type(), file_space);

    H5::DataSpace mem_space(1, &nyp);
    const std::array<hsize_t, 2> mem_size{1, nyp};

    std::vector<Float> row_x(nyp);
    std::vector<Float> row_y(nyp);
    for (hsize_t i = 0; i < nxp; ++i) {
      const std::array<hsize_t, 2> offset{i, 0};
      file_space.selectHyperslab(H5S_SELECT_SET, mem_size.data(), offset.data());

      for (hsize_t j = 0; j < nyp; ++j) {
        switch (m_coords) {
          case Coordinates::CARTESIAN:
            row_x[j] = m_x_min + static_cast<Float>(i) * m_dx;
            row_y[j] = m_y_min + static_cast<Float>(j) * m_dy;
            break;
          case Coordinates::POLAR:
            const auto theta = m_x_min + static_cast<Float>(i) * m_dx;
            const auto r     = m_y_min + static_cast<Float>(j) * m_dy;
            row_x[j]         = r * std::cos(theta);
            row_y[j]         = r * std::sin(theta);
        }
      }

      dataset_x.write(row_x.data(), pred_type(), mem_space, file_space);
      dataset_y.write(row_y.data(), pred_type(), mem_space, file_space);
    }
  }

  // ===============================================================================================
  constexpr void
  write_scalar(const H5::Group& group, const std::string& name, Scalar<Float, LAYOUT> s) {
    static_assert(LAYOUT == Layout::C, "Only Layout::C is supported.");
    const auto nx = static_cast<size_t>(s.nx());
    const auto ny = static_cast<size_t>(s.ny());
    const auto ng = static_cast<size_t>(s.nghost());

    const std::array<hsize_t, 2> mem_size{nx + 2 * ng, ny + 2 * ng};
    const std::array<hsize_t, 2> file_size{nx, ny};
    const std::array<hsize_t, 2> offset{ng, ng};

    H5::DataSpace file_space(2, file_size.data());
    H5::DataSpace mem_space(2, mem_size.data());
    mem_space.selectHyperslab(H5S_SELECT_SET, file_size.data(), offset.data());

    H5::DataSet dataset = group.createDataSet(name, pred_type(), file_space);
    dataset.write(s.data(), pred_type(), mem_space, file_space);
  }

  // ===============================================================================================
  constexpr void write_vector_polar(const H5::Group& group,
                                    const std::string& name_x,
                                    const std::string& name_y,
                                    Vector<Float, LAYOUT> v) {
    const auto nx = static_cast<hsize_t>(v.x.nx());
    const auto ny = static_cast<hsize_t>(v.x.ny());
    const std::array<hsize_t, 2> file_size{nx, ny};
    H5::DataSpace file_space(2, file_size.data());
    H5::DataSet dataset_x = group.createDataSet(name_x, pred_type(), file_space);
    H5::DataSet dataset_y = group.createDataSet(name_y, pred_type(), file_space);

    H5::DataSpace mem_space(1, &ny);
    const std::array<hsize_t, 2> mem_size{1, ny};

    std::vector<Float> row_x(ny);
    std::vector<Float> row_y(ny);
    for (hsize_t i = 0; i < nx; ++i) {
      const std::array<hsize_t, 2> offset{i, 0};
      file_space.selectHyperslab(H5S_SELECT_SET, mem_size.data(), offset.data());

      for (hsize_t j = 0; j < ny; ++j) {
        const double vtheta = v.x(static_cast<Index>(i), static_cast<Index>(j));
        const double vr     = v.y(static_cast<Index>(i), static_cast<Index>(j));
        const double theta  = m_x_min + static_cast<Float>(i) * m_dx;
        row_x[j]            = vr * std::cos(theta) - vtheta * std::sin(theta);
        row_y[j]            = vr * std::sin(theta) + vtheta * std::cos(theta);
      }

      dataset_x.write(row_x.data(), pred_type(), mem_space, file_space);
      dataset_y.write(row_y.data(), pred_type(), mem_space, file_space);
    }
  }

  // ===============================================================================================
  [[nodiscard]] constexpr auto create_group(std::string_view& group_name) -> H5::Group {
    using detail::Indent;

    const auto bytes_written = std::snprintf(  // NOLINT
        group_name_buffer.data(),
        group_name_buffer.size(),
        "%06d",
        write_counter++);
    if (bytes_written < 0 || static_cast<size_t>(bytes_written) + 1 > group_name_buffer.size()) {
      Igor::Panic("Could not generate group name for write_counter = {}", write_counter);
    }
    IGOR_ASSERT(std::strlen(group_name_buffer.data()) == group_name_buffer.size() - 1,
                "Incorrect length of string in group_name_buffer: `{}`",
                group_name_buffer.data());
    group_name = std::string_view(group_name_buffer.data(), group_name_buffer.size() - 1);
    H5::Group group(m_hdf_file.createGroup(group_name.data()));  // NOLINT

    m_xdmf_file << Indent(6) << R"(<Grid Name=")" << group_name << R"(" GridType="Uniform">)"
                << '\n';
    return group;
  }

  constexpr void end_group() {
    using detail::Indent;
    m_xdmf_file << Indent(6) << R"(</Grid>)" << '\n';
  }

  // ===============================================================================================
  constexpr void write_time(const H5::Group& group, Float t) {
    using detail::Indent;

    const std::string dataset_name = "time";
    hsize_t size                   = 1;
    H5::DataSpace data_space(1, &size);
    H5::DataSet dataset = group.createDataSet(dataset_name, pred_type(), data_space);
    dataset.write(&t, pred_type());

    m_xdmf_file << Indent(8) << R"(<Time Value=")"
                << std::setprecision(std::numeric_limits<Float>::max_digits10) << t << R"("/>)"
                << '\n';
  }

  // ===============================================================================================
  constexpr void write_grid(const H5::Group& group,
                            std::string_view group_name,
                            const VertexScalar<Float, LAYOUT>* x = nullptr,
                            const VertexScalar<Float, LAYOUT>* y = nullptr) {
    using detail::Indent;
    using Igor::detail::format;

    const std::string grid_group_name = "Grid";
    H5::Group grid_group              = group.createGroup(grid_group_name);

    if (x == nullptr && y == nullptr) {
      write_grid_unmaterialized(grid_group);
    } else if (x != nullptr && y != nullptr) {
      write_scalar(grid_group, "x", x->scalar());
      write_scalar(grid_group, "y", y->scalar());
    } else {
      Igor::Panic("x and y must either both be nullptr or both not.");
    }

    m_xdmf_file << Indent(8)
                << format(R"(<Topology TopologyType="2DSMesh" Dimensions="{} {}"/>)",
                          m_nx + 1,
                          m_ny + 1)
                << '\n';
    m_xdmf_file << Indent(8) << R"(<Geometry GeometryType="X_Y">)" << '\n';
    m_xdmf_file
        << Indent(10)
        << format(
               R"(<DataItem Dimensions="{} {}" NumberType="Float" Precision="{}" Format="HDF">{}:/{}/{}/x</DataItem>)",
               m_nx + 1,
               m_ny + 1,
               precision(),
               HDF_FILENAME,
               group_name,
               grid_group_name)
        << '\n';
    m_xdmf_file
        << Indent(10)
        << format(
               R"(<DataItem Dimensions="{} {}" NumberType="Float" Precision="{}" Format="HDF">{}:/{}/{}/y</DataItem>)",
               m_nx + 1,
               m_ny + 1,
               precision(),
               HDF_FILENAME,
               group_name,
               grid_group_name)
        << '\n';
    m_xdmf_file << Indent(8) << R"(</Geometry>)" << '\n';
  }

  // ===============================================================================================
  constexpr void write_scalars(const H5::Group& group, std::string_view group_name) {
    using detail::Indent;
    using Igor::detail::format;

    const std::string scalar_group_name = "Scalar";
    H5::Group scalar_group              = group.createGroup(scalar_group_name);
    for (size_t i = 0; i < m_scalar_names.size(); ++i) {
      const auto& values = m_scalar_values[i];
      const auto& name   = m_scalar_names[i];
      write_scalar(scalar_group, name, values);

      m_xdmf_file << Indent(8)
                  << format(R"(<Attribute Name="{}" AttributeType="Scalar" Center="Cell">)", name)
                  << '\n';
      m_xdmf_file
          << Indent(10)
          << format(
                 R"(<DataItem Dimensions="{} {}" NumberType="Float" Precision="{}" Format="HDF">{}:/{}/Scalar/{}</DataItem>)",
                 values.nx(),
                 values.ny(),
                 precision(),
                 HDF_FILENAME,
                 group_name,
                 name)
          << '\n';
      m_xdmf_file << Indent(8) << R"(</Attribute>)" << '\n';
    }
  }

  // ===============================================================================================
  constexpr void write_vectors(const H5::Group& group, std::string_view group_name) {
    using detail::Indent;
    using Igor::detail::format;

    const std::string vector_group_name = "Vector";
    H5::Group vector_group              = group.createGroup(vector_group_name);
    for (size_t i = 0; i < m_vector_names.size(); ++i) {
      const auto& values = m_vector_values[i];
      const auto& name   = m_vector_names[i];
      const auto name_x  = name + "_x";
      const auto name_y  = name + "_y";

      switch (m_coords) {
        case Coordinates::CARTESIAN:
          write_scalar(vector_group, name_x, values.x);
          write_scalar(vector_group, name_y, values.y);
          break;
        case Coordinates::POLAR: write_vector_polar(vector_group, name_x, name_y, values); break;
      }

      m_xdmf_file << Indent(8)
                  << format(R"(<Attribute Name="{}" AttributeType="Vector" Center="Cell">)", name)
                  << '\n';
      m_xdmf_file
          << Indent(10)
          << format(
                 R"a(<DataItem ItemType="Function" Function="JOIN($0, $1, 0*$0)" Dimensions="{} {} 3" NumberType="Float" Precision="{}">)a",
                 values.x.nx(),
                 values.x.ny(),
                 precision())
          << '\n';
      m_xdmf_file
          << Indent(12)
          << format(
                 R"(<DataItem Dimensions="{} {}" NumberType="Float" Precision="{}" Format="HDF">{}:/{}/Vector/{}</DataItem>)",
                 values.x.nx(),
                 values.x.ny(),
                 precision(),
                 HDF_FILENAME,
                 group_name,
                 name_x)
          << '\n';
      m_xdmf_file
          << Indent(12)
          << format(
                 R"(<DataItem Dimensions="{} {}" NumberType="Float" Precision="{}" Format="HDF">{}:/{}/Vector/{}</DataItem>)",
                 values.x.nx(),
                 values.x.ny(),
                 precision(),
                 HDF_FILENAME,
                 group_name,
                 name_y)
          << '\n';
      m_xdmf_file << Indent(10) << R"(</DataItem>)" << '\n';
      m_xdmf_file << Indent(8) << R"(</Attribute>)" << '\n';
    }
  }

  // ===============================================================================================
  constexpr auto write(const VertexScalar<Float, LAYOUT>* x,
                       const VertexScalar<Float, LAYOUT>* y,
                       Float t = -1.0) -> bool {
    try {
      m_xdmf_file.seekp(next_write_pos);

      std::string_view group_name;
      H5::Group group = create_group(group_name);

      write_time(group, t);
      write_grid(group, group_name, x, y);
      write_scalars(group, group_name);
      write_vectors(group, group_name);

      end_group();

      next_write_pos = m_xdmf_file.tellp();
      m_xdmf_file << END_TAGS;

      m_hdf_file.flush(H5F_SCOPE_GLOBAL);
      m_xdmf_file.flush();

      return m_xdmf_file.good();
    } catch (const H5::Exception& e) {
      Igor::Error("Could not write data:");
      e.printErrorStack(stderr);

      // Restore last correct timestep
      m_xdmf_file.seekp(0, std::ofstream::end);
      const auto to_override = m_xdmf_file.tellp() - next_write_pos;
      if (to_override <= 0) {
        Igor::Panic("Expected positive amount of bytes written after last correct step but got {}",
                    to_override);
      }
      m_xdmf_file.seekp(next_write_pos);

      const size_t to_pad = static_cast<size_t>(to_override) - std::strlen(END_TAGS);
      m_xdmf_file << END_TAGS;
      for (size_t i = 0; i < to_pad; ++i) {
        m_xdmf_file << ' ';
      }

      m_hdf_file.flush(H5F_SCOPE_GLOBAL);
      m_xdmf_file.flush();

      return false;
    }
  }

 public:
  constexpr HDFWriter(const std::string& output_dir, const Grid<Float, LAYOUT>& grid)
      : m_hdf_file(output_dir + "/" + HDF_FILENAME, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT),
        m_xdmf_file(output_dir + "/" + XDMF_FILENAME),
        m_x_min(grid.x_min()),
        m_dx(grid.dx()),
        m_nx(grid.nx()),
        m_y_min(grid.y_min()),
        m_dy(grid.dy()),
        m_ny(grid.ny()),
        m_coords(grid.coords()) {
    using detail::Indent;

    H5::Exception::dontPrint();

    m_xdmf_file << R"(<?xml version="1.0" ?>)" << '\n';
    m_xdmf_file << R"(<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>)" << '\n';
    m_xdmf_file << R"(<Xdmf Version="2.0">)" << '\n';
    m_xdmf_file << Indent(2) << R"(<Domain>)" << '\n';
    m_xdmf_file << Indent(4)
                << R"(<Grid Name="state" GridType="Collection" CollectionType="Temporal">)" << '\n';

    next_write_pos = m_xdmf_file.tellp();
    m_xdmf_file << END_TAGS;
  }

  constexpr HDFWriter(const HDFWriter& other) noexcept                    = delete;
  constexpr HDFWriter(HDFWriter&& other) noexcept                         = delete;
  constexpr auto operator=(const HDFWriter& other) noexcept -> HDFWriter& = delete;
  constexpr auto operator=(HDFWriter&& other) noexcept -> HDFWriter&      = delete;
  constexpr ~HDFWriter() noexcept                                         = default;

  // ===============================================================================================
  constexpr void add_field(std::string name, const Scalar<Float, LAYOUT> s) {
    m_scalar_names.emplace_back(std::move(name));
    m_scalar_values.push_back(s);
  }

  constexpr void add_field(std::string name, const Vector<Float, LAYOUT> v) {
    m_vector_names.emplace_back(std::move(name));
    m_vector_values.push_back(v);
  }

  // ===============================================================================================
  constexpr auto write(Float t = -1.0) -> bool { return write(nullptr, nullptr, t); }
  constexpr auto write(const VertexScalar<Float, LAYOUT> x,
                       const VertexScalar<Float, LAYOUT> y,
                       Float t = -1.0) -> bool {
    IGOR_ASSERT(x.nx() == m_nx + 1 && x.ny() == m_ny + 1,
                "Invalid dimensions of x: expected ({}, {}) but got ({}, {})",
                m_nx + 1,
                m_ny + 1,
                x.nx(),
                x.ny());
    IGOR_ASSERT(y.nx() == m_nx + 1 && y.ny() == m_ny + 1,
                "Invalid dimensions of y: expected ({}, {}) but got ({}, {})",
                m_nx + 1,
                m_ny + 1,
                y.nx(),
                y.ny());
    return write(&x, &y, t);
  }
};
