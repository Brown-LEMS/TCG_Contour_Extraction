#pragma once

#include "tcg_types.hpp"

#include <string>
#include <vector>

namespace tcg {

/// Write a simple junction file for visualization / downstream use.
/// Coordinates match .cem (0-based); add +1 when overlaying in MATLAB like draw_contours.
///
/// T_junctions: from classify_junction_type_wrt_graph_BP (degree-3 nodes where two
///   branches were merged as the through-contour).
/// Y_junctions: degree >= 3 nodes remaining in the contour graph (true multi-way
///   junctions that were not classified/merged as T).
bool write_jct(const std::string& path, const std::vector<Edge>& T_junctions,
               const std::vector<Edge>& Y_junctions, int h, int w, std::string& error);

}  // namespace tcg
