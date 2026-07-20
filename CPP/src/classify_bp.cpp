#include "classify_bp.hpp"
#include "fac_graph.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <set>
#include <unordered_set>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace tcg {
namespace {

Contour orient_toward_junction(const Contour& c, const ContourEdgeIndices& ids, int eid, int nbr) {
  Contour out;
  if (ids.empty() || c.empty()) return out;
  const int cur_nbr = std::min(nbr, static_cast<int>(ids.size()));
  if (ids.front() == eid) {
    out.assign(c.begin(), c.begin() + cur_nbr);
    out = reverse_contour(out);
  } else if (ids.back() == eid) {
    out.assign(c.end() - cur_nbr, c.end());
  }
  return out;
}

bool has_unique_nbrs(const std::vector<int>& nbrs) {
  std::unordered_set<int> s(nbrs.begin(), nbrs.end());
  return s.size() == nbrs.size();
}

}  // namespace

ClassifyBPResult classify_junction_type_wrt_graph_BP(const std::vector<Contour>& cfrags_in,
                                                     const std::vector<ContourEdgeIndices>& cfrags_idx_in,
                                                     const std::vector<Edge>& edges,
                                                     const ClassifyBPParams& params) {
  const double angle_diff_th = params.BP_merge_angle_th;
  const int nbr_num_edges = params.BP_nbr_num_edges;
  const double clen_th = params.BP_clen_th;
  const double w0 = 1.0;

  ClassifyBPResult result;
  result.contours = cfrags_in;
  result.contour_edge_idx = cfrags_idx_in;
  result.G = construct_fac_graph(result.contour_edge_idx);

  auto& G = result.G;
  auto& merged_cfrags = result.contours;
  auto& merged_cfrags_idx = result.contour_edge_idx;

  // Phase A: pairwise costs at junctions degree >= 3
  for (size_t vid = 0; vid < G.vars.size(); ++vid) {
    auto& var = G.vars[vid];
    const int v_degree = static_cast<int>(var.nbrs_fac.size());
    if (v_degree < 3) continue;
    if (!has_unique_nbrs(var.nbrs_fac)) continue;

    var.cfrags_at_v.assign(static_cast<size_t>(v_degree), Contour{});
    for (int i = 0; i < v_degree; ++i) {
      const int cid = var.nbrs_fac[static_cast<size_t>(i)];
      var.cfrags_at_v[static_cast<size_t>(i)] =
          orient_toward_junction(merged_cfrags[static_cast<size_t>(cid)],
                                 merged_cfrags_idx[static_cast<size_t>(cid)], var.actual_edge_id,
                                 nbr_num_edges);
    }

    var.angle_diff_mat.assign(static_cast<size_t>(v_degree),
                              std::vector<double>(static_cast<size_t>(v_degree), M_PI));
    for (int i = 0; i < v_degree; ++i) {
      for (int j = i + 1; j < v_degree; ++j) {
        var.angle_diff_mat[static_cast<size_t>(i)][static_cast<size_t>(j)] =
            co_circular_cost(var.cfrags_at_v[static_cast<size_t>(i)],
                             var.cfrags_at_v[static_cast<size_t>(j)]);
      }
    }
  }

  // Phase B: one BP-style message pass over short inter-junction fragments
  for (size_t cid = 0; cid < G.facs.size(); ++cid) {
    const int vid1 = G.facs[cid].nbrs_var[0];
    const int vid2 = G.facs[cid].nbrs_var[1];
    auto& var1 = G.vars[static_cast<size_t>(vid1)];
    auto& var2 = G.vars[static_cast<size_t>(vid2)];
    const int v_degree1 = static_cast<int>(var1.nbrs_fac.size());
    const int v_degree2 = static_cast<int>(var2.nbrs_fac.size());
    if (v_degree1 < 3 || v_degree2 < 3) continue;
    if (!has_unique_nbrs(var1.nbrs_fac) || !has_unique_nbrs(var2.nbrs_fac)) continue;

    std::unordered_set<int> s1(var1.nbrs_fac.begin(), var1.nbrs_fac.end());
    int inter = 0;
    for (int x : var2.nbrs_fac)
      if (s1.count(x)) ++inter;
    if (inter > 1) continue;

    if (contour_length(merged_cfrags[cid]) > clen_th) continue;
    if (var1.angle_diff_mat.empty() || var2.angle_diff_mat.empty()) continue;
    if (var1.cfrags_at_v.size() != static_cast<size_t>(v_degree1) ||
        var2.cfrags_at_v.size() != static_cast<size_t>(v_degree2))
      continue;

    const Contour& c0 = merged_cfrags[cid];

    auto find_ind = [](const std::vector<int>& nbrs, int cid_q) {
      for (size_t i = 0; i < nbrs.size(); ++i)
        if (nbrs[i] == cid_q) return static_cast<int>(i);
      return -1;
    };

    // Update J1
    const int ind_c1 = find_ind(var1.nbrs_fac, static_cast<int>(cid));
    if (ind_c1 < 0) continue;
    for (int i = 0; i < v_degree1; ++i) {
      if (var1.nbrs_fac[static_cast<size_t>(i)] == static_cast<int>(cid)) continue;
      double min_angle_diff = M_PI;
      const Contour& c1 = var1.cfrags_at_v[static_cast<size_t>(i)];
      for (int j = 0; j < v_degree2; ++j) {
        if (var2.nbrs_fac[static_cast<size_t>(j)] == static_cast<int>(cid)) continue;
        const Contour& c2 = var2.cfrags_at_v[static_cast<size_t>(j)];
        Contour c_merge = c2;
        if (!c_merge.empty()) c_merge.pop_back();
        Contour c0_rev = reverse_contour(c0);
        c_merge.insert(c_merge.end(), c0_rev.begin(), c0_rev.end());
        min_angle_diff = std::min(min_angle_diff, co_circular_cost(c1, c_merge) * w0);
      }
      if (i < ind_c1) {
        var1.angle_diff_mat[static_cast<size_t>(i)][static_cast<size_t>(ind_c1)] =
            std::min(var1.angle_diff_mat[static_cast<size_t>(i)][static_cast<size_t>(ind_c1)],
                     min_angle_diff);
      } else {
        var1.angle_diff_mat[static_cast<size_t>(ind_c1)][static_cast<size_t>(i)] =
            std::min(var1.angle_diff_mat[static_cast<size_t>(ind_c1)][static_cast<size_t>(i)],
                     min_angle_diff);
      }
    }

    // Update J2
    const int ind_c2 = find_ind(var2.nbrs_fac, static_cast<int>(cid));
    if (ind_c2 < 0) continue;
    for (int j = 0; j < v_degree2; ++j) {
      if (var2.nbrs_fac[static_cast<size_t>(j)] == static_cast<int>(cid)) continue;
      double min_angle_diff = M_PI;
      const Contour& c2 = var2.cfrags_at_v[static_cast<size_t>(j)];
      for (int i = 0; i < v_degree1; ++i) {
        if (var1.nbrs_fac[static_cast<size_t>(i)] == static_cast<int>(cid)) continue;
        const Contour& c1 = var1.cfrags_at_v[static_cast<size_t>(i)];
        Contour c_merge = c1;
        if (!c_merge.empty()) c_merge.pop_back();
        c_merge.insert(c_merge.end(), c0.begin(), c0.end());
        min_angle_diff = std::min(min_angle_diff, co_circular_cost(c2, c_merge) * w0);
      }
      if (j < ind_c2) {
        var2.angle_diff_mat[static_cast<size_t>(j)][static_cast<size_t>(ind_c2)] =
            std::min(var2.angle_diff_mat[static_cast<size_t>(j)][static_cast<size_t>(ind_c2)],
                     min_angle_diff);
      } else {
        var2.angle_diff_mat[static_cast<size_t>(ind_c2)][static_cast<size_t>(j)] =
            std::min(var2.angle_diff_mat[static_cast<size_t>(ind_c2)][static_cast<size_t>(j)],
                     min_angle_diff);
      }
    }
  }

  // Phase C: T-junction merges at degree-3 nodes
  std::set<int> cid_to_remove;
  for (size_t vid = 0; vid < G.vars.size(); ++vid) {
    auto& var = G.vars[vid];
    auto nbrs_fac = var.nbrs_fac;
    const int v_degree = static_cast<int>(nbrs_fac.size());
    if (v_degree != 3) continue;
    if (!has_unique_nbrs(nbrs_fac)) continue;
    if (var.angle_diff_mat.empty()) continue;

    const int eid = var.actual_edge_id;
    if (eid < 0 || static_cast<size_t>(eid) >= edges.size()) continue;

    std::vector<ContourEdgeIndices> cfrags_idx_at_v(static_cast<size_t>(v_degree));
    for (int i = 0; i < v_degree; ++i) {
      ContourEdgeIndices ids = merged_cfrags_idx[static_cast<size_t>(nbrs_fac[static_cast<size_t>(i)])];
      if (!ids.empty() && ids.front() == eid) ids = reverse_indices(ids);
      cfrags_idx_at_v[static_cast<size_t>(i)] = std::move(ids);
    }

    double min_angle_diff = std::numeric_limits<double>::infinity();
    int min_i = -1, min_j = -1;
    int n_min = 0;
    for (int i = 0; i < v_degree; ++i) {
      for (int j = i + 1; j < v_degree; ++j) {
        const double a = var.angle_diff_mat[static_cast<size_t>(i)][static_cast<size_t>(j)];
        if (a < min_angle_diff - 1e-15) {
          min_angle_diff = a;
          min_i = i;
          min_j = j;
          n_min = 1;
        } else if (std::abs(a - min_angle_diff) <= 1e-15) {
          ++n_min;
        }
      }
    }
    if (n_min > 1 || min_i < 0) continue;
    if (!(min_angle_diff < angle_diff_th)) continue;

    result.T_junctions.push_back(edges[static_cast<size_t>(eid)]);

    const int cid_1 = nbrs_fac[static_cast<size_t>(min_i)];
    const int cid_2 = nbrs_fac[static_cast<size_t>(min_j)];
    ContourEdgeIndices c1_ids = cfrags_idx_at_v[static_cast<size_t>(min_i)];
    ContourEdgeIndices c2_ids = reverse_indices(cfrags_idx_at_v[static_cast<size_t>(min_j)]);
    ContourEdgeIndices cfrag_idx = c1_ids;
    if (c2_ids.size() > 1) cfrag_idx.insert(cfrag_idx.end(), c2_ids.begin() + 1, c2_ids.end());

    Contour cfrag;
    cfrag.reserve(cfrag_idx.size());
    for (int id : cfrag_idx) {
      if (id >= 0 && static_cast<size_t>(id) < edges.size()) cfrag.push_back(edges[static_cast<size_t>(id)]);
    }
    merged_cfrags[static_cast<size_t>(cid_1)] = std::move(cfrag);
    merged_cfrags_idx[static_cast<size_t>(cid_1)] = std::move(cfrag_idx);

    // Detach c1,c2 from v (remove both indices; remove larger first)
    if (min_i > min_j) {
      var.nbrs_fac.erase(var.nbrs_fac.begin() + min_i);
      var.nbrs_fac.erase(var.nbrs_fac.begin() + min_j);
    } else {
      var.nbrs_fac.erase(var.nbrs_fac.begin() + min_j);
      var.nbrs_fac.erase(var.nbrs_fac.begin() + min_i);
    }
    var.dim = static_cast<int>(var.nbrs_fac.size());

    const int vid_2 = (G.facs[static_cast<size_t>(cid_2)].nbrs_var[0] == static_cast<int>(vid))
                          ? G.facs[static_cast<size_t>(cid_2)].nbrs_var[1]
                          : G.facs[static_cast<size_t>(cid_2)].nbrs_var[0];
    for (int& nv : G.facs[static_cast<size_t>(cid_1)].nbrs_var) {
      if (nv == static_cast<int>(vid)) nv = vid_2;
    }
    for (int& nf : G.vars[static_cast<size_t>(vid_2)].nbrs_fac) {
      if (nf == cid_2) nf = cid_1;
    }
    G.vars[static_cast<size_t>(vid_2)].dim =
        static_cast<int>(G.vars[static_cast<size_t>(vid_2)].nbrs_fac.size());

    cid_to_remove.insert(cid_2);
  }

  std::vector<int> rem(cid_to_remove.begin(), cid_to_remove.end());
  std::sort(rem.begin(), rem.end(), std::greater<int>());
  for (int cid : rem) {
    if (cid < 0 || static_cast<size_t>(cid) >= merged_cfrags.size()) continue;
    merged_cfrags.erase(merged_cfrags.begin() + cid);
    merged_cfrags_idx.erase(merged_cfrags_idx.begin() + cid);
  }
  remove_empty_contours(merged_cfrags, merged_cfrags_idx);
  result.G = construct_fac_graph(merged_cfrags_idx);
  return result;
}

}  // namespace tcg
