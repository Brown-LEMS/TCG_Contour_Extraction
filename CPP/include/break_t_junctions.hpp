#pragma once

#include "tcg_types.hpp"

namespace tcg {

/// Port of break_contours_at_T_junctions.m
/// Splits contours at internal edges with d2f (degree) > 1.
struct TBreakResult {
  std::vector<Contour> contours;
  std::vector<ContourEdgeIndices> contour_edge_idx;
  std::vector<Edge> edges;
};

TBreakResult break_contours_at_T_junctions(const std::vector<Contour>& cfrags,
                                           const std::vector<ContourEdgeIndices>& cfrags_idx,
                                           const std::vector<Edge>& edges);

}  // namespace tcg
