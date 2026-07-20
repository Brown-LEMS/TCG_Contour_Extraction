#include "contour_breaker.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>

namespace tcg {

namespace {

double edgel_dist(const Edge& a, const Edge& b) {
  double dx = a.x - b.x;
  double dy = a.y - b.y;
  return std::sqrt(dx * dx + dy * dy);
}

double contour_length(const Contour& c) {
  if (c.size() < 2) return 0;
  double sum = 0;
  for (size_t i = 0; i + 1 < c.size(); ++i) sum += edgel_dist(c[i], c[i + 1]);
  return sum;
}

//> MATLAB smooth(y): default method 'moving' with span 5 (not 0.05*n).
//  Endpoints use a truncated window. Even span is reduced by 1 (MATLAB rule).
std::vector<double> smooth_moving_default(const std::vector<double>& y) {
  int n = static_cast<int>(y.size());
  if (n == 0) return y;
  int span = 5;
  if (span % 2 == 0) --span;
  if (span < 1) span = 1;
  if (span > n) {
    span = n;
    if (span % 2 == 0 && span > 1) --span;
  }
  int half = span / 2;
  std::vector<double> out(static_cast<size_t>(n));
  for (int i = 0; i < n; ++i) {
    int lo = std::max(0, i - half);
    int hi = std::min(n - 1, i + half);
    double s = 0;
    int cnt = 0;
    for (int j = lo; j <= hi; ++j) {
      s += y[static_cast<size_t>(j)];
      ++cnt;
    }
    out[static_cast<size_t>(i)] = s / static_cast<double>(cnt);
  }
  return out;
}

/// filter_co_circular_along_cfrag.m
void filter_co_circular_along_cfrag(const Contour& cfrag, double nbr_range_th_d,
                                    std::vector<double>& ori_diff_vec, std::vector<double>& dtheta) {
  int len = static_cast<int>(cfrag.size());
  ori_diff_vec.assign(static_cast<size_t>(len), 0.0);
  dtheta.assign(static_cast<size_t>(len), 0.0);
  if (len < 3) return;

  int nbr = static_cast<int>(std::round(nbr_range_th_d));
  if (nbr < 1) nbr = 1;

  std::vector<double> X(static_cast<size_t>(len)), Y(static_cast<size_t>(len));
  for (int i = 0; i < len; ++i) {
    X[static_cast<size_t>(i)] = cfrag[static_cast<size_t>(i)].x + 1.0;
    Y[static_cast<size_t>(i)] = cfrag[static_cast<size_t>(i)].y + 1.0;
  }

  std::vector<double> DX0(static_cast<size_t>(len - 2)), DY0(static_cast<size_t>(len - 2));
  for (int i = 0; i < len - 2; ++i) {
    DX0[static_cast<size_t>(i)] = X[static_cast<size_t>(i + 2)] - X[static_cast<size_t>(i)];
    DY0[static_cast<size_t>(i)] = Y[static_cast<size_t>(i + 2)] - Y[static_cast<size_t>(i)];
  }
  std::vector<double> DX(static_cast<size_t>(len)), DY(static_cast<size_t>(len));
  DX[0] = DX0[0];
  DY[0] = DY0[0];
  for (int i = 1; i < len - 1; ++i) {
    DX[static_cast<size_t>(i)] = DX0[static_cast<size_t>(i - 1)];
    DY[static_cast<size_t>(i)] = DY0[static_cast<size_t>(i - 1)];
  }
  DX[static_cast<size_t>(len - 1)] = DX0[static_cast<size_t>(len - 3)];
  DY[static_cast<size_t>(len - 1)] = DY0[static_cast<size_t>(len - 3)];

  std::vector<double> theta_vec(static_cast<size_t>(len));
  for (int i = 0; i < len; ++i) {
    theta_vec[static_cast<size_t>(i)] = std::atan2(DY[static_cast<size_t>(i)], DX[static_cast<size_t>(i)]);
  }

  dtheta[0] = 0;
  for (int i = 1; i < len; ++i) {
    double dt = std::abs(theta_vec[static_cast<size_t>(i)] - theta_vec[static_cast<size_t>(i - 1)]);
    if (dt > M_PI) dt = 2 * M_PI - dt;
    dtheta[static_cast<size_t>(i)] = dt;
  }

  for (int i = nbr; i <= len - nbr - 1; ++i) {
    double sum_a = 0;
    for (int k = 0; k < nbr; ++k) {
      double dx1 = DX[static_cast<size_t>(i - nbr + k)];
      double dy1 = DY[static_cast<size_t>(i - nbr + k)];
      double dx2 = DX[static_cast<size_t>(i + 1 + k)];
      double dy2 = DY[static_cast<size_t>(i + 1 + k)];
      double dx = X[static_cast<size_t>(i + 1 + k)] - X[static_cast<size_t>(i - nbr + k)];
      double dy = Y[static_cast<size_t>(i + 1 + k)] - Y[static_cast<size_t>(i - nbr + k)];

      double a = std::abs(2 * std::atan2(dy, dx) - std::atan2(dy1, dx1) - std::atan2(dy2, dx2));
      a = std::fmod(a, 2 * M_PI);
      if (a > M_PI) a = 2 * M_PI - a;
      sum_a += a;
    }
    ori_diff_vec[static_cast<size_t>(i)] = sum_a / static_cast<double>(nbr);
  }
}

/// MATLAB [~, id] = max(v) returns the first maximum index only.
int argmax_first(const std::vector<double>& v) {
  if (v.empty()) return 0;
  int best_i = 0;
  double best = v[0];
  for (int i = 1; i < static_cast<int>(v.size()); ++i) {
    if (v[static_cast<size_t>(i)] > best) {
      best = v[static_cast<size_t>(i)];
      best_i = i;
    }
  }
  return best_i;
}

bool rows_all_equal_5(const Edge& a, const Edge& b) {
  return a.x == b.x && a.y == b.y && a.dir == b.dir && a.conf == b.conf && a.d2f == b.d2f;
}

}  // namespace

CornerBreakResult contour_breaker_at_corner(const std::vector<Contour>& cfrags_in,
                                            const std::vector<ContourEdgeIndices>& cfrags_idx_in,
                                            const BreakerParams& params, double ori_diff_th) {
  CornerBreakResult result;
  if (ori_diff_th < 0) ori_diff_th = params.corner_angle_th;

  const int nbr_num_edges = params.nbr_num_edges;
  const int iter = 1;

  std::vector<Contour> new_cfrags = cfrags_in;
  std::vector<ContourEdgeIndices> new_cfrags_idx = cfrags_idx_in;

  size_t i = 0;
  while (i < new_cfrags.size()) {
    Contour& cur_c = new_cfrags[i];
    ContourEdgeIndices& cur_c_idx = new_cfrags_idx[i];

    double c_len = contour_length(cur_c);

    if (!cur_c.empty() && c_len < (2 * nbr_num_edges) && rows_all_equal_5(cur_c.front(), cur_c.back())) {
      ++i;
      continue;
    }

    bool split_done = false;
    if (c_len > nbr_num_edges && static_cast<int>(cur_c.size()) > nbr_num_edges) {
      double cur_nbr_num =
          std::min(std::ceil(nbr_num_edges / 2.0), std::ceil(static_cast<double>(cur_c.size()) / 6.0)) /
          static_cast<double>(iter);
      std::vector<double> ori_diff_vec, Dtheta;
      filter_co_circular_along_cfrag(cur_c, cur_nbr_num, ori_diff_vec, Dtheta);

      int nbr_i = static_cast<int>(std::round(cur_nbr_num));
      int two_n = 2 * nbr_i;
      const int len = static_cast<int>(ori_diff_vec.size());
      for (int k = 0; k < two_n && k < len; ++k) ori_diff_vec[static_cast<size_t>(k)] = 0;
      if (len > 0 && two_n >= 0) {
        int start0 = len - two_n - 1;
        if (start0 < 0) start0 = 0;
        for (int k = start0; k < len; ++k) ori_diff_vec[static_cast<size_t>(k)] = 0;
      }

      ori_diff_vec = smooth_moving_default(ori_diff_vec);

      for (size_t k = 0; k < Dtheta.size(); ++k) {
        if (Dtheta[k] < (ori_diff_th / 2.0)) ori_diff_vec[k] = Dtheta[k];
      }

      int id0 = argmax_first(ori_diff_vec);
      double max_ori_diff = ori_diff_vec[static_cast<size_t>(id0)];
      int id_matlab = id0 + 1;

      if (max_ori_diff > ori_diff_th && id_matlab >= 1 && id_matlab <= static_cast<int>(cur_c.size())) {
        size_t split = static_cast<size_t>(id_matlab - 1);

        Contour head(cur_c.begin(), cur_c.begin() + split + 1);
        Contour tail(cur_c.begin() + split, cur_c.end());

        ContourEdgeIndices head_idx, tail_idx;
        if (!cur_c_idx.empty() && cur_c_idx.size() == cur_c.size()) {
          head_idx.assign(cur_c_idx.begin(), cur_c_idx.begin() + static_cast<std::ptrdiff_t>(split + 1));
          tail_idx.assign(cur_c_idx.begin() + static_cast<std::ptrdiff_t>(split), cur_c_idx.end());
        }

        // Capture corner before push_back: reallocation would invalidate cur_c / cur_c_idx refs.
        Edge corner = tail.front();
        cur_c = std::move(tail);
        cur_c_idx = std::move(tail_idx);
        new_cfrags.push_back(std::move(head));
        new_cfrags_idx.push_back(std::move(head_idx));

        result.corner_pts.push_back(corner);
        split_done = true;
      }
    }
    if (split_done) continue;
    ++i;
  }

  result.contours = std::move(new_cfrags);
  result.contour_edge_idx = std::move(new_cfrags_idx);
  return result;
}

}  // namespace tcg
