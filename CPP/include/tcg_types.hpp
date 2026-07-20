#pragma once

#include <array>
#include <cstddef>
#include <string>
#include <vector>

namespace tcg {

/// One edgel / contour vertex: x, y, direction, confidence, d2f (matches .cem columns).
struct Edge {
  double x{};
  double y{};
  double dir{};
  double conf{};
  double d2f{};
};

using Contour = std::vector<Edge>;
/// Per-fragment list of edge indices into the global edge table (0-based, C++ convention).
using ContourEdgeIndices = std::vector<int>;

struct EdgFile {
  int version{1};  // 1 or 2 (v3.0 header maps to version 2 layout in MATLAB loader)
  int width{};
  int height{};
  std::vector<Edge> edges;
  /// Optional maps (MATLAB row = y+1, col = x+1); stored row-major [y][x], size height x width.
  std::vector<double> edgemap;
  std::vector<double> thetamap;
};

struct CemFile {
  int width{};
  int height{};
  std::vector<Edge> edges;
  std::vector<Contour> contours;
  std::vector<ContourEdgeIndices> contour_edge_idx;
  /// Optional 8 properties per contour from [Contour Properties]; may be empty.
  std::vector<std::array<double, 8>> contour_props;
};

struct BreakerParams {
  int nbr_num_edges{20};
  /// Default corner threshold (radians); can be overridden per call.
  double corner_angle_th{3.14159265358979323846 / 6.0};
};

}  // namespace tcg
