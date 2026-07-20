#include "fac_graph.hpp"

#include <algorithm>
#include <cmath>
#include <unordered_map>
#include <utility>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace tcg {

FacGraph construct_fac_graph(const std::vector<ContourEdgeIndices>& cf_idx) {
  FacGraph G;
  std::unordered_map<int, int> edge_to_vid;
  for (size_t c = 0; c < cf_idx.size(); ++c) {
    if (cf_idx[c].empty()) continue;
    const int e1 = cf_idx[c].front();
    const int e2 = cf_idx[c].back();
    int id1, id2;
    auto it1 = edge_to_vid.find(e1);
    if (it1 == edge_to_vid.end()) {
      id1 = static_cast<int>(G.vars.size());
      FacVar v;
      v.id = id1;
      v.actual_edge_id = e1;
      G.vars.push_back(v);
      edge_to_vid[e1] = id1;
    } else {
      id1 = it1->second;
    }
    auto it2 = edge_to_vid.find(e2);
    if (it2 == edge_to_vid.end()) {
      id2 = static_cast<int>(G.vars.size());
      FacVar v;
      v.id = id2;
      v.actual_edge_id = e2;
      G.vars.push_back(v);
      edge_to_vid[e2] = id2;
    } else {
      id2 = it2->second;
    }
    FacFac f;
    f.id = static_cast<int>(G.facs.size());
    f.nbrs_var[0] = id1;
    f.nbrs_var[1] = id2;
    G.facs.push_back(f);
    G.vars[static_cast<size_t>(id1)].nbrs_fac.push_back(f.id);
    G.vars[static_cast<size_t>(id1)].dim = static_cast<int>(G.vars[static_cast<size_t>(id1)].nbrs_fac.size());
    G.vars[static_cast<size_t>(id2)].nbrs_fac.push_back(f.id);
    G.vars[static_cast<size_t>(id2)].dim = static_cast<int>(G.vars[static_cast<size_t>(id2)].nbrs_fac.size());
  }
  return G;
}

double contour_length(const Contour& c) {
  if (c.size() < 2) return 0.0;
  double sum = 0.0;
  for (size_t i = 0; i + 1 < c.size(); ++i) {
    const double dx = c[i + 1].x - c[i].x;
    const double dy = c[i + 1].y - c[i].y;
    sum += std::sqrt(dx * dx + dy * dy);
  }
  return sum;
}

void interpolate_cfrag(const Contour& cfrag, std::vector<double>& x_coords, std::vector<double>& y_coords) {
  x_coords.clear();
  y_coords.clear();
  const double c_len = contour_length(cfrag);
  std::vector<std::pair<double, double>> path;
  path.reserve(cfrag.size());
  for (const auto& e : cfrag) {
    if (!path.empty() && path.back().first == e.x && path.back().second == e.y) continue;
    path.emplace_back(e.x, e.y);
  }
  if (path.size() < 2) {
    if (!path.empty()) {
      x_coords.push_back(path[0].first);
      y_coords.push_back(path[0].second);
    }
    return;
  }
  std::vector<double> cum(path.size(), 0.0);
  for (size_t i = 1; i < path.size(); ++i) {
    const double dx = path[i].first - path[i - 1].first;
    const double dy = path[i].second - path[i - 1].second;
    cum[i] = cum[i - 1] + std::sqrt(dx * dx + dy * dy);
  }
  const int n_samp = std::max(2, static_cast<int>(std::lround(c_len) * 2));
  x_coords.resize(static_cast<size_t>(n_samp));
  y_coords.resize(static_cast<size_t>(n_samp));
  for (int s = 0; s < n_samp; ++s) {
    const double t = (n_samp == 1) ? 0.0 : (cum.back() * static_cast<double>(s) / (n_samp - 1));
    size_t i = 0;
    while (i + 1 < cum.size() && cum[i + 1] < t) ++i;
    double x = path[i].first, y = path[i].second;
    if (i + 1 < path.size()) {
      const double seg = cum[i + 1] - cum[i];
      const double a = (seg > 0) ? (t - cum[i]) / seg : 0.0;
      x = path[i].first + a * (path[i + 1].first - path[i].first);
      y = path[i].second + a * (path[i + 1].second - path[i].second);
    }
    x_coords[static_cast<size_t>(s)] = x;
    y_coords[static_cast<size_t>(s)] = y;
  }
}

double co_circular_cost(const Contour& c1_in, const Contour& c2_in) {
  // Flip c2: c1 ------> <---------c2  =>  c1 ------> -------> c2
  Contour c2 = reverse_contour(c2_in);
  const Contour& c1 = c1_in;

  auto local_angle = [](const Contour& a, const Contour& b, int len) {
    const int la = static_cast<int>(a.size());
    const int lb = static_cast<int>(b.size());
    const int L = std::min({la, lb, len});
    if (L < 1) return M_PI;
    const double dx1 = a.back().x - a[static_cast<size_t>(la - L)].x;
    const double dy1 = a.back().y - a[static_cast<size_t>(la - L)].y;
    const double dx2 = b[static_cast<size_t>(L - 1)].x - b[0].x;
    const double dy2 = b[static_cast<size_t>(L - 1)].y - b[0].y;
    const double n1 = std::sqrt(dx1 * dx1 + dy1 * dy1);
    const double n2 = std::sqrt(dx2 * dx2 + dy2 * dy2);
    if (n1 < 1e-12 || n2 < 1e-12) return M_PI;
    double cosv = (dx1 * dx2 + dy1 * dy2) / (n1 * n2);
    cosv = std::max(-1.0, std::min(1.0, cosv));
    return std::acos(cosv);
  };

  double local_diff_alpha = local_angle(c1, c2, 5);
  if (static_cast<int>(c1.size()) <= 5 || static_cast<int>(c2.size()) <= 5) {
    return local_diff_alpha;
  }

  std::vector<double> x1s, y1s, x2s, y2s;
  interpolate_cfrag(c1, x1s, y1s);
  interpolate_cfrag(c2, x2s, y2s);
  Contour ic1(x1s.size()), ic2(x2s.size());
  for (size_t i = 0; i < x1s.size(); ++i) {
    ic1[i].x = x1s[i];
    ic1[i].y = y1s[i];
  }
  for (size_t i = 0; i < x2s.size(); ++i) {
    ic2[i].x = x2s[i];
    ic2[i].y = y2s[i];
  }
  local_diff_alpha = local_angle(ic1, ic2, 15);

  const int len1 = static_cast<int>(ic1.size()) - 2;
  const int len2 = static_cast<int>(ic2.size()) - 2;
  if (len1 <= 0 || len2 <= 0) return local_diff_alpha;

  double sum = 0.0;
  int count = 0;
  for (int i = 0; i < len1; ++i) {
    const double x1 = ic1[static_cast<size_t>(i)].x;
    const double y1 = ic1[static_cast<size_t>(i)].y;
    const double dx1 = ic1[static_cast<size_t>(i + 2)].x - x1;
    const double dy1 = ic1[static_cast<size_t>(i + 2)].y - y1;
    for (int j = 0; j < len2; ++j) {
      const double x2 = ic2[static_cast<size_t>(j + 2)].x;
      const double y2 = ic2[static_cast<size_t>(j + 2)].y;
      const double dx2 = x2 - ic2[static_cast<size_t>(j)].x;
      const double dy2 = y2 - ic2[static_cast<size_t>(j)].y;
      const double dx = x2 - x1;
      const double dy = y2 - y1;
      double a = std::abs(2.0 * std::atan2(dy, dx) - std::atan2(dy1, dx1) - std::atan2(dy2, dx2));
      a = std::fmod(a, 2.0 * M_PI);
      if (a < 0) a += 2.0 * M_PI;
      if (a > M_PI) a = 2.0 * M_PI - a;
      sum += a;
      ++count;
    }
  }
  double diff_alpha = (count > 0) ? (sum / static_cast<double>(count)) : local_diff_alpha;
  return std::min(local_diff_alpha, diff_alpha);
}

Contour reverse_contour(const Contour& c) {
  Contour out = c;
  std::reverse(out.begin(), out.end());
  return out;
}

ContourEdgeIndices reverse_indices(const ContourEdgeIndices& ids) {
  ContourEdgeIndices out = ids;
  std::reverse(out.begin(), out.end());
  return out;
}

static bool xy_equal(const Edge& a, const Edge& node) {
  return a.x == node.x && a.y == node.y;
}

Contour merge_two_curve_fragments(const Contour& c1_in, const Contour& c2_in, const Edge& node) {
  Contour c1 = c1_in;
  Contour c2 = c2_in;
  if (!c1.empty() && xy_equal(c1.front(), node)) c1 = reverse_contour(c1);
  if (!c2.empty() && xy_equal(c2.back(), node)) c2 = reverse_contour(c2);
  Contour out = c1;
  if (c2.size() > 1) out.insert(out.end(), c2.begin() + 1, c2.end());
  return out;
}

void remove_empty_contours(std::vector<Contour>& cfrags, std::vector<ContourEdgeIndices>& cfrags_idx) {
  std::vector<Contour> c2;
  std::vector<ContourEdgeIndices> i2;
  c2.reserve(cfrags.size());
  i2.reserve(cfrags_idx.size());
  for (size_t i = 0; i < cfrags.size(); ++i) {
    if (cfrags_idx[i].empty()) continue;
    c2.push_back(std::move(cfrags[i]));
    i2.push_back(std::move(cfrags_idx[i]));
  }
  cfrags = std::move(c2);
  cfrags_idx = std::move(i2);
}

}  // namespace tcg
