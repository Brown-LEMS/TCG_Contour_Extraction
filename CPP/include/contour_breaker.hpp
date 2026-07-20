#pragma once

#include "tcg_types.hpp"

#include <string>

namespace tcg {

/// Port of contour_breaker_at_conner.m (+ filter_co_circular_along_cfrag, contour length, smooth).
/// @param ori_diff_th  if negative, uses params.corner_angle_th
struct CornerBreakResult {
  std::vector<Contour> contours;
  std::vector<ContourEdgeIndices> contour_edge_idx;
  std::vector<Edge> corner_pts;
};

CornerBreakResult contour_breaker_at_corner(const std::vector<Contour>& cfrags,
                                            const std::vector<ContourEdgeIndices>& cfrags_idx,
                                            const BreakerParams& params, double ori_diff_th = -1.0);

}  // namespace tcg
