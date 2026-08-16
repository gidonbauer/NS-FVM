#ifndef NS_FVM_IO_HPP_
#define NS_FVM_IO_HPP_

#include <cstring>
#include <filesystem>

#include <Igor/Logging.hpp>

// -------------------------------------------------------------------------------------------------
template <typename Float>
[[nodiscard]] constexpr auto should_save(Float t, Float dt, Float dt_write, Float t_end) -> bool {
  constexpr Float DT_SAFE      = 1e-6;
  static Float last_save_t     = -1.0;

  const bool dt_write_complete = std::fmod(t + DT_SAFE * dt, dt_write) < dt * (1.0 - DT_SAFE);
  const bool is_last           = std::abs(t - t_end) < DT_SAFE;
  const bool res               = dt_write_complete || is_last;
  if (res && is_last && std::abs(t - last_save_t) < DT_SAFE) { return false; }
  if (res) { last_save_t = t; }
  return res;
}

// -------------------------------------------------------------------------------------------------
#ifndef NS_FVM_BASE_DIR
#define NS_FVM_BASE_DIR "."
#endif  // FS_BASE_DIR
[[nodiscard]] auto
get_output_directory(std::string_view subdir   = "output",
                     std::string_view base_dir = NS_FVM_BASE_DIR,
                     std::source_location loc  = std::source_location::current()) noexcept
    -> std::string {
  // ===============
  auto strip_path = [](std::string_view full_path) -> std::string_view {
#if defined(WIN32) || defined(_WIN32)
#error "Not implemented yet: Requires different path separator ('\\') and potentially uses wchar"
#else
    constexpr char separator = '/';
#endif
    size_t counter = 0;
    for (char c : std::ranges::reverse_view(full_path)) {
      if (c == separator) { break; }
      ++counter;
    }

    return full_path.substr(full_path.size() - counter, counter);
  };

  // ===============
  auto strip_extension = [](std::string_view full_name) -> std::string_view {
    constexpr char separator = '.';
    size_t counter           = 0;
    for (char c : std::ranges::reverse_view(full_name)) {
      if (c == separator) { break; }
      ++counter;
    }

    return full_name.substr(0, full_name.size() - counter - 1);
  };

  // ===============
  try {
    return Igor::detail::format(
        "{}/{}/{}/", base_dir, subdir, strip_extension(strip_path(loc.file_name())));
  } catch (const std::exception& e) {
    Igor::Panic("Could not create output directory name.");
    std::unreachable();
  }
}

// -------------------------------------------------------------------------------------------------
[[nodiscard]] auto init_output_directory(std::string_view directory_name) noexcept -> bool {
  std::error_code ec;

  std::filesystem::remove_all(directory_name, ec);
  if (ec) {
    Igor::Warn("Could remove directory `{}`: {}", directory_name, ec.message());
    return false;
  }

  std::filesystem::create_directories(directory_name, ec);
  if (ec) {
    Igor::Warn("Could not create directory `{}`: {}", directory_name, ec.message());
    return false;
  }

  return true;
}

#endif  // NS_FVM_IO_HPP_
