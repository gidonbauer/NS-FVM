#pragma once

#include <array>

#include "Grid.hpp"

namespace Expected {

template <std::size_t N>
constexpr auto interp_n2(const std::array<Index, N>& xs, const std::array<double, N>& es, Index n)
    -> double {
  static_assert(N >= 2, "need at least two samples");

  // Exact hits: return the tabulated value bit-for-bit.
  for (std::size_t i = 0; i < N; ++i) {
    if (xs[i] == n) { return es[i]; }
  }

  const auto coeff = [&](std::size_t i) {
    const auto xi = static_cast<double>(xs[i]);
    return es[i] * xi * xi;
  };
  const auto nd = static_cast<double>(n);

  // Extrapolation: hold C at the nearest end -> pure quadratic scaling.
  if (n < xs.front()) { return coeff(0) / (nd * nd); }
  if (n > xs.back()) { return coeff(N - 1) / (nd * nd); }

  // Bracket: xs[i] < n < xs[i + 1]
  std::size_t i = 0;
  while (xs[i + 1] < n) {
    ++i;
  }

  const auto x0 = static_cast<double>(xs[i]);
  const auto x1 = static_cast<double>(xs[i + 1]);
  const auto t  = (nd - x0) / (x1 - x0);

  return ((1.0 - t) * coeff(i) + t * coeff(i + 1)) / (nd * nd);
}

}  // namespace Expected
