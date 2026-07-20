#include "merge_geom.hpp"
#include "fac_graph.hpp"

#include <algorithm>
#include <cmath>

namespace tcg {
namespace {

void merge_at_degree_2_node(FacGraph& G, std::vector<Contour>& merged_cem,
                            std::vector<ContourEdgeIndices>& merged_cf_idx, int v,
                            const Edge& actual_edge, const std::vector<int>& nbrs_fac,
                            const ContourEdgeIndices& c1_ids, const ContourEdgeIndices& c2_ids) {
  const int f1 = nbrs_fac[0];
  const int f2 = nbrs_fac[1];
  Contour c_merged =
      merge_two_curve_fragments(merged_cem[static_cast<size_t>(f1)], merged_cem[static_cast<size_t>(f2)],
                                actual_edge);
  ContourEdgeIndices c_ids_merged = c1_ids;
  if (c2_ids.size() > 1) c_ids_merged.insert(c_ids_merged.end(), c2_ids.begin() + 1, c2_ids.end());

  merged_cem[static_cast<size_t>(f1)] = std::move(c_merged);
  merged_cf_idx[static_cast<size_t>(f1)] = std::move(c_ids_merged);
  merged_cem[static_cast<size_t>(f2)].clear();
  merged_cf_idx[static_cast<size_t>(f2)].clear();

  G.vars[static_cast<size_t>(v)].nbrs_fac.clear();
  G.vars[static_cast<size_t>(v)].dim = 0;
  G.vars[static_cast<size_t>(v)].merged = 1;

  const int n1_id = (G.facs[static_cast<size_t>(f1)].nbrs_var[0] == v)
                        ? G.facs[static_cast<size_t>(f1)].nbrs_var[1]
                        : G.facs[static_cast<size_t>(f1)].nbrs_var[0];
  const int n3_id = (G.facs[static_cast<size_t>(f2)].nbrs_var[0] == v)
                        ? G.facs[static_cast<size_t>(f2)].nbrs_var[1]
                        : G.facs[static_cast<size_t>(f2)].nbrs_var[0];

  G.facs[static_cast<size_t>(f1)].nbrs_var[0] = n1_id;
  G.facs[static_cast<size_t>(f1)].nbrs_var[1] = n3_id;
  G.facs[static_cast<size_t>(f2)].removed = 1;

  auto& nf3 = G.vars[static_cast<size_t>(n3_id)].nbrs_fac;
  for (int& x : nf3) {
    if (x == f2) x = f1;
  }
  G.vars[static_cast<size_t>(n3_id)].dim = static_cast<int>(nf3.size());
}

}  // namespace

MergeGeomResult merge_cfrags_graphical_model_geom(const std::vector<Contour>& cfrags_in,
                                                  const std::vector<ContourEdgeIndices>& cfrags_idx_in,
                                                  const std::vector<Edge>& edges,
                                                  const MergeGeomParams& params, double angle_diff_th) {
  if (angle_diff_th < 0) angle_diff_th = params.geom_merge_angle_th;
  const int nbr_num_edges = params.nbr_num_edges;

  MergeGeomResult result;
  result.contours = cfrags_in;
  result.contour_edge_idx = cfrags_idx_in;
  result.G = construct_fac_graph(result.contour_edge_idx);

  auto& G = result.G;
  auto& merged_cem = result.contours;
  auto& merged_cf_idx = result.contour_edge_idx;

  for (size_t v = 0; v < G.vars.size(); ++v) {
    if (G.vars[v].dim != 2) continue;
    auto nbrs_fac = G.vars[v].nbrs_fac;
    if (nbrs_fac.size() != 2) continue;
    if (nbrs_fac[0] == nbrs_fac[1]) continue;

    ContourEdgeIndices c1_ids = merged_cf_idx[static_cast<size_t>(nbrs_fac[0])];
    ContourEdgeIndices c2_ids = merged_cf_idx[static_cast<size_t>(nbrs_fac[1])];
    Contour c1 = merged_cem[static_cast<size_t>(nbrs_fac[0])];
    Contour c2 = merged_cem[static_cast<size_t>(nbrs_fac[1])];
    if (c1.empty() || c2.empty() || c1_ids.empty() || c2_ids.empty()) continue;

    const int cur_nbr =
        std::min({nbr_num_edges, static_cast<int>(c1.size()), static_cast<int>(c2.size())});
    if (cur_nbr < 1) continue;

    const int eid = G.vars[v].actual_edge_id;

    // Orient truncated portions as c1 -> node <- c2; flip full index lists when needed.
    Contour c1_cost, c2_cost;
    ContourEdgeIndices c1_ids_orient = c1_ids;
    ContourEdgeIndices c2_ids_orient = c2_ids;

    if (c1_ids.front() == eid) {
      c1_cost.assign(c1.begin(), c1.begin() + cur_nbr);
      c1_cost = reverse_contour(c1_cost);
      c1_ids_orient = reverse_indices(c1_ids);
    } else if (c1_ids.back() == eid) {
      c1_cost.assign(c1.end() - cur_nbr, c1.end());
    } else {
      continue;
    }

    if (c2_ids.front() == eid) {
      c2_cost.assign(c2.begin(), c2.begin() + cur_nbr);
      c2_cost = reverse_contour(c2_cost);
      c2_ids_orient = reverse_indices(c2_ids);
    } else if (c2_ids.back() == eid) {
      c2_cost.assign(c2.end() - cur_nbr, c2.end());
    } else {
      continue;
    }

    const double angle_diff = co_circular_cost(c1_cost, c2_cost);
    bool is_merge = false;
    if (cur_nbr <= 5 && angle_diff < 2.0 * angle_diff_th) is_merge = true;
    else if (angle_diff < angle_diff_th) is_merge = true;
    G.vars[v].p = std::exp(-angle_diff);

    if (is_merge) {
      if (eid < 0 || static_cast<size_t>(eid) >= edges.size()) continue;
      merge_at_degree_2_node(G, merged_cem, merged_cf_idx, static_cast<int>(v), edges[static_cast<size_t>(eid)],
                             nbrs_fac, c1_ids_orient, reverse_indices(c2_ids_orient));
    }
  }

  remove_empty_contours(merged_cem, merged_cf_idx);
  result.G = construct_fac_graph(merged_cf_idx);
  return result;
}

}  // namespace tcg
