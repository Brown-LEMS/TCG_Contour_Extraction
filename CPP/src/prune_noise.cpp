#include "prune_noise.hpp"
#include "fac_graph.hpp"

#include <algorithm>
#include <cmath>
#include <set>
#include <unordered_set>

namespace tcg {
namespace {

inline size_t idx2(int y, int x, int w) {
  return static_cast<size_t>(y) * static_cast<size_t>(w) + static_cast<size_t>(x);
}

double mean_soft(const Contour& c, int h, int w, const std::vector<double>& soft) {
  if (c.empty() || soft.empty()) return 0.0;
  double s = 0.0;
  int n = 0;
  for (const auto& e : c) {
    int yi = static_cast<int>(std::lround(e.y)) + 1;
    int xi = static_cast<int>(std::lround(e.x)) + 1;
    yi = std::max(1, std::min(h, yi));
    xi = std::max(1, std::min(w, xi));
    s += soft[idx2(yi - 1, xi - 1, w)];
    ++n;
  }
  return (n > 0) ? s / static_cast<double>(n) : 0.0;
}

int other_var(const FacFac& f, int vid) {
  if (f.nbrs_var[0] == vid) return f.nbrs_var[1];
  if (f.nbrs_var[1] == vid) return f.nbrs_var[0];
  return -1;
}

}  // namespace

PruneResult prune_noise_curves(const std::vector<Contour>& cfrags_in,
                               const std::vector<ContourEdgeIndices>& cfrags_idx_in, int h, int w,
                               const std::vector<double>& edgemap_soft, const PruneParams& params) {
  PruneResult result;
  result.contours = cfrags_in;
  result.contour_edge_idx = cfrags_idx_in;

  const double len_th = params.noise_len_th;
  const double prob_th = params.noise_prob_th;

  FacGraph G = construct_fac_graph(result.contour_edge_idx);
  std::set<int> cid_to_remove;

  for (size_t vid = 0; vid < G.vars.size(); ++vid) {
    const auto& var = G.vars[vid];
    const int deg = static_cast<int>(var.nbrs_fac.size());

    if (deg == 1) {
      const int cid = var.nbrs_fac[0];
      if (contour_length(result.contours[static_cast<size_t>(cid)]) < len_th) {
        cid_to_remove.insert(cid);
      }
    } else if (deg == 2) {
      const int cid1 = var.nbrs_fac[0];
      const int cid2 = var.nbrs_fac[1];
      if (cid1 == cid2) continue;
      if (static_cast<int>(result.contour_edge_idx[static_cast<size_t>(cid1)].size()) > len_th ||
          static_cast<int>(result.contour_edge_idx[static_cast<size_t>(cid2)].size()) > len_th) {
        continue;
      }
      const int cvid1 = other_var(G.facs[static_cast<size_t>(cid1)], static_cast<int>(vid));
      const int cvid2 = other_var(G.facs[static_cast<size_t>(cid2)], static_cast<int>(vid));
      if (cvid1 < 0 || cvid2 < 0) continue;
      const auto& nf1 = G.vars[static_cast<size_t>(cvid1)].nbrs_fac;
      const auto& nf2 = G.vars[static_cast<size_t>(cvid2)].nbrs_fac;
      std::unordered_set<int> s1(nf1.begin(), nf1.end());
      bool inter = false;
      for (int x : nf2) {
        if (s1.count(x)) {
          inter = true;
          break;
        }
      }
      if (inter) {
        cid_to_remove.insert(cid1);
        cid_to_remove.insert(cid2);
      }
    } else if (deg == 3) {
      for (int cid : var.nbrs_fac) {
        if (static_cast<int>(result.contour_edge_idx[static_cast<size_t>(cid)].size()) > 2 * len_th)
          continue;
        const int cvid = other_var(G.facs[static_cast<size_t>(cid)], static_cast<int>(vid));
        if (cvid < 0) continue;
        for (int ccid : G.vars[static_cast<size_t>(cvid)].nbrs_fac) {
          if (ccid == cid) continue;
          if (static_cast<int>(result.contour_edge_idx[static_cast<size_t>(ccid)].size()) > 2 * len_th)
            continue;
          const int ccvid = other_var(G.facs[static_cast<size_t>(ccid)], cvid);
          if (ccvid == static_cast<int>(vid)) {
            const double p1 = mean_soft(result.contours[static_cast<size_t>(cid)], h, w, edgemap_soft);
            const double p2 = mean_soft(result.contours[static_cast<size_t>(ccid)], h, w, edgemap_soft);
            cid_to_remove.insert(p1 < p2 ? cid : ccid);
          }
        }
      }
    }
  }

  for (size_t cid = 0; cid < G.facs.size(); ++cid) {
    const auto& cur_c_idx = result.contour_edge_idx[cid];
    if (cur_c_idx.empty()) continue;
    if (cur_c_idx.front() == cur_c_idx.back()) {
      if (contour_length(result.contours[cid]) < 3.0 * len_th) cid_to_remove.insert(static_cast<int>(cid));
      continue;
    }
    const int v0 = G.facs[cid].nbrs_var[0];
    const int v1 = G.facs[cid].nbrs_var[1];
    if (G.vars[static_cast<size_t>(v0)].nbrs_fac.size() == 1 &&
        G.vars[static_cast<size_t>(v1)].nbrs_fac.size() == 1) {
      const double c_len = contour_length(result.contours[cid]);
      const double avg_prob = mean_soft(result.contours[cid], h, w, edgemap_soft);
      if ((avg_prob < prob_th && c_len < 3.0 * len_th) || c_len < len_th) {
        cid_to_remove.insert(static_cast<int>(cid));
      }
    }
  }

  // Remove in descending order so indices stay valid (MATLAB unique + delete).
  std::vector<int> rem(cid_to_remove.begin(), cid_to_remove.end());
  std::sort(rem.begin(), rem.end(), std::greater<int>());
  for (int cid : rem) {
    if (cid < 0 || static_cast<size_t>(cid) >= result.contours.size()) continue;
    result.contours.erase(result.contours.begin() + cid);
    result.contour_edge_idx.erase(result.contour_edge_idx.begin() + cid);
  }
  return result;
}

}  // namespace tcg
