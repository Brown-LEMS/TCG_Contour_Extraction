#pragma once

#include "fac_graph.hpp"
#include "tcg_types.hpp"

namespace tcg {

struct MergeGeomParams {
  double geom_merge_angle_th{kTcgParams.geom_merge_angle_th};
  int nbr_num_edges{kTcgParams.nbr_num_edges};
};

struct MergeGeomResult {
  FacGraph G;
  std::vector<Contour> contours;
  std::vector<ContourEdgeIndices> contour_edge_idx;
};

/// Port of merge_cfrags_graphical_model_geom.m
MergeGeomResult merge_cfrags_graphical_model_geom(const std::vector<Contour>& cfrags,
                                                  const std::vector<ContourEdgeIndices>& cfrags_idx,
                                                  const std::vector<Edge>& edges,
                                                  const MergeGeomParams& params,
                                                  double angle_diff_th = -1.0);

}  // namespace tcg
