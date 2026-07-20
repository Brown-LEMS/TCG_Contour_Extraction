#pragma once

#include "tcg_types.hpp"

namespace tcg {

struct GapFillResult {
  std::vector<Contour> contours;
  std::vector<ContourEdgeIndices> contour_edge_idx;
  std::vector<Edge> edges;
};

/// Port of util/contour_fill_gaps_DP.m (Step 4).
GapFillResult contour_fill_gaps_DP(const std::vector<Contour>& cfrags,
                                   const std::vector<ContourEdgeIndices>& cfrags_idx, int h, int w,
                                   const std::vector<Edge>& edges, const GradientMaps& grad,
                                   const GapFillParams& params);

}  // namespace tcg
