#pragma once

#include "fac_graph.hpp"
#include "tcg_types.hpp"

namespace tcg {

struct ClassifyBPParams {
  double BP_merge_angle_th{kTcgParams.BP_merge_angle_th};
  int BP_nbr_num_edges{kTcgParams.BP_nbr_num_edges};
  double BP_clen_th{kTcgParams.BP_clen_th};
};

struct ClassifyBPResult {
  FacGraph G;
  std::vector<Contour> contours;
  std::vector<ContourEdgeIndices> contour_edge_idx;
  std::vector<Edge> T_junctions;
};

/// Port of classify_junction_type_wrt_graph_BP.m
ClassifyBPResult classify_junction_type_wrt_graph_BP(const std::vector<Contour>& cfrags,
                                                     const std::vector<ContourEdgeIndices>& cfrags_idx,
                                                     const std::vector<Edge>& edges,
                                                     const ClassifyBPParams& params);

}  // namespace tcg
