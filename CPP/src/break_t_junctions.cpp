#include "break_t_junctions.hpp"

namespace tcg {

TBreakResult break_contours_at_T_junctions(const std::vector<Contour>& cfrags_in,
                                           const std::vector<ContourEdgeIndices>& cfrags_idx_in,
                                           const std::vector<Edge>& edges_in) {
  TBreakResult result;
  result.edges = edges_in;
  result.contours = cfrags_in;
  result.contour_edge_idx = cfrags_idx_in;

  // Process only the original set; appends go to the end (matching MATLAB).
  const size_t n_orig = result.contours.size();
  for (size_t c = 0; c < n_orig; ++c) {
    Contour& cur_c = result.contours[c];
    ContourEdgeIndices& cur_c_idx = result.contour_edge_idx[c];
    if (cur_c_idx.size() < 2) continue;

    std::vector<double> edge_degree(cur_c_idx.size(), 0.0);
    for (size_t i = 0; i < cur_c_idx.size(); ++i) {
      const int eid = cur_c_idx[i];
      if (eid >= 0 && static_cast<size_t>(eid) < result.edges.size()) {
        edge_degree[i] = result.edges[static_cast<size_t>(eid)].d2f;
      }
    }
    edge_degree.front() = 0.0;
    edge_degree.back() = 0.0;

    std::vector<int> breaks;  // 0-based positions along contour
    for (size_t i = 0; i < edge_degree.size(); ++i) {
      if (edge_degree[i] > 1.0) breaks.push_back(static_cast<int>(i));
    }
    if (breaks.empty()) continue;

    // MATLAB 1-based: cur_c(breaks(end):end), cur_c(1:breaks(1)), then mid segments
    const int last_b = breaks.back();
    const int first_b = breaks.front();

    Contour tail(cur_c.begin() + last_b, cur_c.end());
    ContourEdgeIndices tail_idx(cur_c_idx.begin() + last_b, cur_c_idx.end());
    Contour head(cur_c.begin(), cur_c.begin() + first_b + 1);
    ContourEdgeIndices head_idx(cur_c_idx.begin(), cur_c_idx.begin() + first_b + 1);

    cur_c = std::move(tail);
    cur_c_idx = std::move(tail_idx);
    result.contours.push_back(std::move(head));
    result.contour_edge_idx.push_back(std::move(head_idx));

    int prev_e = first_b;
    for (size_t j = 0; j + 1 < breaks.size(); ++j) {
      const int next_e = breaks[j + 1];
      Contour mid(cfrags_in[c].begin() + prev_e, cfrags_in[c].begin() + next_e + 1);
      ContourEdgeIndices mid_idx(cfrags_idx_in[c].begin() + prev_e,
                                 cfrags_idx_in[c].begin() + next_e + 1);
      result.contours.push_back(std::move(mid));
      result.contour_edge_idx.push_back(std::move(mid_idx));
      prev_e = next_e;
    }
  }
  return result;
}

}  // namespace tcg
