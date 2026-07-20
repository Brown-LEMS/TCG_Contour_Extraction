#pragma once

#include "tcg_types.hpp"

namespace tcg {

struct PruneParams {
  double noise_len_th{5.0};
  double noise_prob_th{0.05};
};

struct PruneResult {
  std::vector<Contour> contours;
  std::vector<ContourEdgeIndices> contour_edge_idx;
};

/// Port of prune_noise_curves.m
/// edgemap_soft: row-major h*w (from .edg load or imgradient).
PruneResult prune_noise_curves(const std::vector<Contour>& cfrags,
                               const std::vector<ContourEdgeIndices>& cfrags_idx, int h, int w,
                               const std::vector<double>& edgemap_soft, const PruneParams& params);

}  // namespace tcg
