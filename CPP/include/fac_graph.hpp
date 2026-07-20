#pragma once

#include "tcg_types.hpp"

#include <vector>

namespace tcg {

struct FacVar {
  int id{};                  // 0-based
  int actual_edge_id{};      // 0-based into edges
  std::vector<int> nbrs_fac; // 0-based factor ids
  int dim{};                 // == nbrs_fac.size()
  int merged{0};
  double p{0};
  std::vector<std::vector<double>> angle_diff_mat;
  std::vector<Contour> cfrags_at_v;
};

struct FacFac {
  int id{};
  int nbrs_var[2]{};  // endpoint var ids (0-based)
  int removed{0};
};

struct FacGraph {
  std::vector<FacVar> vars;
  std::vector<FacFac> facs;
};

/// Port of construct_fac_graph_from_curve_fragments.m (0-based indices).
FacGraph construct_fac_graph(const std::vector<ContourEdgeIndices>& cf_idx);

double contour_length(const Contour& c);

/// Port of interpolate_cfrag.m — returns resampled (x,y) path.
void interpolate_cfrag(const Contour& cfrag, std::vector<double>& x_coords, std::vector<double>& y_coords);

/// Port of co_circular_cost.m. Expects orientation c1 ------> <--------- c2 (toward junction).
double co_circular_cost(const Contour& c1, const Contour& c2);

Contour reverse_contour(const Contour& c);
ContourEdgeIndices reverse_indices(const ContourEdgeIndices& ids);

/// Orient full fragments as c1 --> node --> c2 and concatenate (skip duplicate node in c2).
Contour merge_two_curve_fragments(const Contour& c1, const Contour& c2, const Edge& node);

void remove_empty_contours(std::vector<Contour>& cfrags, std::vector<ContourEdgeIndices>& cfrags_idx);

}  // namespace tcg
