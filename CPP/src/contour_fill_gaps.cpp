#include "contour_fill_gaps.hpp"
#include "dp_gap.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <opencv2/core.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv2/ximgproc.hpp>
#include <unordered_map>
#include <utility>

namespace tcg {
namespace {

inline size_t idx2(int y, int x, int w) {
  return static_cast<size_t>(y) * static_cast<size_t>(w) + static_cast<size_t>(x);
}

inline int clampi(int v, int lo, int hi) {
  return std::max(lo, std::min(hi, v));
}

double contour_length(const Contour& c) {
  if (c.size() < 2) return 0;
  double sum = 0;
  for (size_t i = 0; i + 1 < c.size(); ++i) {
    const double dx = c[i + 1].x - c[i].x;
    const double dy = c[i + 1].y - c[i].y;
    sum += std::sqrt(dx * dx + dy * dy);
  }
  return sum;
}

cv::Mat morph_clean(const cv::Mat& bin) {
  // Remove isolated foreground pixels (no 8-neighbors).
  cv::Mat out = bin.clone();
  for (int y = 0; y < bin.rows; ++y) {
    for (int x = 0; x < bin.cols; ++x) {
      if (!bin.at<uchar>(y, x)) continue;
      int nbr = 0;
      for (int dy = -1; dy <= 1; ++dy) {
        for (int dx = -1; dx <= 1; ++dx) {
          if (dx == 0 && dy == 0) continue;
          const int yy = y + dy, xx = x + dx;
          if (yy < 0 || yy >= bin.rows || xx < 0 || xx >= bin.cols) continue;
          if (bin.at<uchar>(yy, xx)) ++nbr;
        }
      }
      if (nbr == 0) out.at<uchar>(y, x) = 0;
    }
  }
  return out;
}

cv::Mat morph_skel(const cv::Mat& bin) {
  cv::Mat src = bin.clone();
  // OpenCV thinning expects 0/255.
  for (int y = 0; y < src.rows; ++y)
    for (int x = 0; x < src.cols; ++x)
      if (src.at<uchar>(y, x)) src.at<uchar>(y, x) = 255;
  cv::Mat dst;
  // Guo-Hall is the closer match to MATLAB bwmorph(...,'skel').
  cv::ximgproc::thinning(src, dst, cv::ximgproc::THINNING_GUOHALL);
  for (int y = 0; y < dst.rows; ++y)
    for (int x = 0; x < dst.cols; ++x)
      dst.at<uchar>(y, x) = dst.at<uchar>(y, x) ? 1 : 0;
  return dst;
}

bool is_junction_nbhd(const uchar* p /* 3x3 row-major */) {
  // Kovesi numbering mapped to row-major: 0 1 2 / 3 4 5 / 6 7 8
  // Perimeter order matching findendsjunctions: 1,2,3,6,9,8,7,4 -> indices 0,3,6,7,8,5,2,1
  const int order[8] = {0, 3, 6, 7, 8, 5, 2, 1};
  int crossings = 0;
  for (int i = 0; i < 8; ++i) {
    const int a = p[order[i]] ? 1 : 0;
    const int b = p[order[(i + 1) % 8]] ? 1 : 0;
    crossings += std::abs(a - b);
  }
  return p[4] && crossings >= 6;
}

std::vector<int> find_junction_map(const cv::Mat& edge_bin /* 0/1 */, int h, int w) {
  cv::Mat sk = morph_skel(edge_bin);
  std::vector<int> jct(static_cast<size_t>(h) * static_cast<size_t>(w), 0);
  for (int y = 1; y + 1 < h; ++y) {
    for (int x = 1; x + 1 < w; ++x) {
      uchar nb[9];
      int k = 0;
      for (int dy = -1; dy <= 1; ++dy)
        for (int dx = -1; dx <= 1; ++dx) nb[k++] = sk.at<uchar>(y + dy, x + dx);
      if (is_junction_nbhd(nb)) jct[idx2(y, x, w)] = 1;
    }
  }
  // Dilate by 3x3 ones like conv2(...,'same')
  std::vector<int> dil = jct;
  for (int y = 0; y < h; ++y) {
    for (int x = 0; x < w; ++x) {
      if (!jct[idx2(y, x, w)]) continue;
      for (int dy = -1; dy <= 1; ++dy) {
        for (int dx = -1; dx <= 1; ++dx) {
          const int yy = y + dy, xx = x + dx;
          if (yy < 0 || yy >= h || xx < 0 || xx >= w) continue;
          dil[idx2(yy, xx, w)] = 1;
        }
      }
    }
  }
  return dil;
}

std::vector<int> convert_cfrags_to_EdgeGroupMap(const std::vector<Contour>& cfrags, int h, int w) {
  std::vector<int> label(static_cast<size_t>(h) * static_cast<size_t>(w), 0);
  const double ratio = 2.0;
  for (size_t ci = 0; ci < cfrags.size(); ++ci) {
    const Contour& c = cfrags[ci];
    if (c.size() < 2) continue;
    // unique rows stable on (x,y)
    std::vector<std::pair<double, double>> path;
    path.reserve(c.size());
    for (const auto& e : c) {
      if (!path.empty() && path.back().first == e.x && path.back().second == e.y) continue;
      path.emplace_back(e.x, e.y);
    }
    if (path.size() < 2) continue;

    const double c_len = contour_length(c);
    std::vector<double> cum(path.size(), 0.0);
    for (size_t i = 1; i < path.size(); ++i) {
      const double dx = path[i].first - path[i - 1].first;
      const double dy = path[i].second - path[i - 1].second;
      cum[i] = cum[i - 1] + std::sqrt(dx * dx + dy * dy);
    }
    const int n_samp = std::max(2, static_cast<int>(std::lround(c_len * ratio)));
    cv::Mat binary = cv::Mat::zeros(h, w, CV_8UC1);
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
      const int xi = clampi(static_cast<int>(std::lround(x)) + 1, 1, w) - 1;
      const int yi = clampi(static_cast<int>(std::lround(y)) + 1, 1, h) - 1;
      binary.at<uchar>(yi, xi) = 1;
    }
    cv::Mat sk = morph_skel(binary);
    const int cid1 = static_cast<int>(ci) + 1;  // MATLAB 1-based contour id
    for (int y = 0; y < h; ++y)
      for (int x = 0; x < w; ++x)
        if (sk.at<uchar>(y, x)) label[idx2(y, x, w)] = cid1;
  }
  return label;
}

struct FacVar {
  int id{};                 // 0-based
  int actual_edge_id{};     // 0-based into edges
  std::vector<int> nbrs_fac;  // 0-based contour/factor ids
};

struct FacFac {
  int id{};
  int nbrs_var[2]{};  // endpoint var ids (0-based)
};

struct FacGraph {
  std::vector<FacVar> vars;
  std::vector<FacFac> facs;
};

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
    G.vars[static_cast<size_t>(id2)].nbrs_fac.push_back(f.id);
  }
  return G;
}

void build_endpoint_maps(const FacGraph& G, const std::vector<Edge>& edges, int h, int w,
                         std::vector<int>& end_point_eid_map, std::vector<int>& end_point_vid_map) {
  // Store MATLAB-style 1-based ids (0 = empty).
  end_point_eid_map.assign(static_cast<size_t>(h) * static_cast<size_t>(w), 0);
  end_point_vid_map.assign(static_cast<size_t>(h) * static_cast<size_t>(w), 0);
  for (const auto& v : G.vars) {
    const Edge& e = edges[static_cast<size_t>(v.actual_edge_id)];
    int ey = clampi(static_cast<int>(std::lround(e.y)) + 1, 1, h) - 1;
    int ex = clampi(static_cast<int>(std::lround(e.x)) + 1, 1, w) - 1;
    end_point_eid_map[idx2(ey, ex, w)] = v.actual_edge_id + 1;
    end_point_vid_map[idx2(ey, ex, w)] = v.id + 1;
  }
}

void make_nbr_grids(int R, std::vector<int>& xx, std::vector<int>& yy, std::vector<double>& rr,
                    std::vector<double>& th) {
  // Match MATLAB: nbr_xx = repmat(-R:R,[2R+1,1]); nbr_yy = nbr_xx'; then (:) column-major.
  const int n = 2 * R + 1;
  xx.clear();
  yy.clear();
  rr.clear();
  th.clear();
  xx.reserve(static_cast<size_t>(n * n));
  yy.reserve(static_cast<size_t>(n * n));
  rr.reserve(static_cast<size_t>(n * n));
  th.reserve(static_cast<size_t>(n * n));
  for (int c = 0; c < n; ++c) {      // column
    for (int r = 0; r < n; ++r) {    // row
      const int dx = c - R;
      const int dy = r - R;
      xx.push_back(dx);
      yy.push_back(dy);
      rr.push_back(std::sqrt(static_cast<double>(dx * dx + dy * dy)));
      th.push_back(std::atan2(static_cast<double>(dy), static_cast<double>(dx)));
    }
  }
}

bool angle_in_fan(double nbr_th, double angle_min, double angle_max) {
  if (angle_max > M_PI) {
    return (nbr_th >= angle_min && nbr_th <= M_PI) || (nbr_th > -M_PI && nbr_th <= angle_max - 2 * M_PI);
  }
  if (angle_min < -M_PI) {
    return (nbr_th >= angle_min + 2 * M_PI && nbr_th <= M_PI) || (nbr_th > -M_PI && nbr_th <= angle_max);
  }
  return nbr_th >= angle_min && nbr_th <= angle_max;
}

enum class GeomPass { First, Second };

void geometric_completion_pass(FacGraph& G, std::vector<Contour>& cfrags,
                               std::vector<ContourEdgeIndices>& cfrags_idx, std::vector<Edge>& edges,
                               std::vector<int>& EdgeGroupMap, std::vector<int>& jct_map,
                               std::vector<int>& end_point_eid_map, std::vector<int>& end_point_vid_map,
                               int h, int w, int search_nbr_range, double search_ori_range, GeomPass pass) {
  std::vector<int> xx, yy;
  std::vector<double> rr, th;
  make_nbr_grids(search_nbr_range, xx, yy, rr, th);
  const double line_scale = (pass == GeomPass::First) ? 0.5 : 0.75;
  const double dist_scale = (pass == GeomPass::First) ? 1.0 : 2.0;

  for (size_t vid = 0; vid < G.vars.size(); ++vid) {
    FacVar& var = G.vars[vid];
    if (var.nbrs_fac.size() != 1) {
      edges[static_cast<size_t>(var.actual_edge_id)].d2f = static_cast<double>(var.nbrs_fac.size());
      continue;
    }
    edges[static_cast<size_t>(var.actual_edge_id)].d2f =
        std::max(edges[static_cast<size_t>(var.actual_edge_id)].d2f, 1.0);
    if (edges[static_cast<size_t>(var.actual_edge_id)].d2f > 1.0) continue;

    int cur_x = clampi(static_cast<int>(std::lround(edges[static_cast<size_t>(var.actual_edge_id)].x)) + 1, 1, w);
    int cur_y = clampi(static_cast<int>(std::lround(edges[static_cast<size_t>(var.actual_edge_id)].y)) + 1, 1, h);

    if (jct_map[idx2(cur_y - 1, cur_x - 1, w)]) {
      if (pass == GeomPass::First)
        edges[static_cast<size_t>(var.actual_edge_id)].d2f = 3.0;
      else
        edges[static_cast<size_t>(var.actual_edge_id)].d2f = 1.0;
      continue;
    }

    const int cid = var.nbrs_fac[0];
    Contour& cur_c = cfrags[static_cast<size_t>(cid)];
    const double c_len0 = contour_length(cur_c);
    const int e_size = static_cast<int>(cur_c.size());
    if (pass == GeomPass::First && e_size < 10) continue;

    bool is_start = true;
    double end_dx = 0, end_dy = 0;
    if (G.facs[static_cast<size_t>(cid)].nbrs_var[0] == static_cast<int>(vid)) {
      const int j = std::min(e_size, 5) - 1;
      end_dx = cur_c[0].x - cur_c[static_cast<size_t>(j)].x;
      end_dy = cur_c[0].y - cur_c[static_cast<size_t>(j)].y;
    } else {
      is_start = false;
      const int j = std::max(e_size - 5, 0);
      end_dx = cur_c.back().x - cur_c[static_cast<size_t>(j)].x;
      end_dy = cur_c.back().y - cur_c[static_cast<size_t>(j)].y;
    }

    const double angle = std::atan2(end_dy, end_dx);
    std::vector<char> qualify(xx.size(), 0);
    auto mark_fan = [&](double amax, double amin, double rmax) {
      for (size_t k = 0; k < xx.size(); ++k) {
        if (rr[k] <= rmax && angle_in_fan(th[k], amin, amax)) qualify[k] = 1;
      }
    };
    mark_fan(angle + search_ori_range, angle - search_ori_range, search_nbr_range);
    mark_fan(angle + M_PI / 3.0, angle - M_PI / 3.0, 2.0);

    int add_end_edge_idx1 = -1;
    double dist_1 = 0;
    std::vector<int> reach_vids;
    for (size_t k = 0; k < xx.size(); ++k) {
      if (!qualify[k]) continue;
      const int xe = clampi(xx[k] + cur_x, 1, w);
      const int ye = clampi(yy[k] + cur_y, 1, h);
      const int rvid1 = end_point_vid_map[idx2(ye - 1, xe - 1, w)];
      if (rvid1 <= 0) continue;
      const int rvid0 = rvid1 - 1;
      if (rvid0 == static_cast<int>(vid)) continue;
      reach_vids.push_back(rvid0);
    }
    if (!reach_vids.empty()) {
      add_end_edge_idx1 = G.vars[static_cast<size_t>(reach_vids[0])].actual_edge_id;
      const Edge& re = edges[static_cast<size_t>(add_end_edge_idx1)];
      dist_1 = std::sqrt(std::pow(cur_x - re.x - 1.0, 2) + std::pow(cur_y - re.y - 1.0, 2));
      if (dist_1 > dist_scale * c_len0) add_end_edge_idx1 = -1;
    }

    int add_end_edge_idx2 = -1;
    double dist_2 = 0;
    const double line_end_x = cur_x + std::cos(angle) * search_nbr_range * line_scale;
    const double line_end_y = cur_y + std::sin(angle) * search_nbr_range * line_scale;
    const int nline = search_nbr_range;
    for (int i = 0; i < nline; ++i) {
      const double t = (nline == 1) ? 0.0 : static_cast<double>(i) / (nline - 1);
      const double lx = std::min(std::max(cur_x + t * (line_end_x - cur_x), 1.0), static_cast<double>(w));
      const double ly = std::min(std::max(cur_y + t * (line_end_y - cur_y), 1.0), static_cast<double>(h));
      const int ix = static_cast<int>(std::lround(lx));
      const int iy = static_cast<int>(std::lround(ly));
      int reach_cid1 = EdgeGroupMap[idx2(iy - 1, ix - 1, w)];
      if (reach_cid1 == cid + 1) continue;
      if (reach_cid1 <= 0) continue;
      const int reach_cid0 = reach_cid1 - 1;
      const Contour& attach = cfrags[static_cast<size_t>(reach_cid0)];
      const auto& attach_idx = cfrags_idx[static_cast<size_t>(reach_cid0)];
      double best = std::numeric_limits<double>::infinity();
      int r_idx = 0;
      for (size_t j = 0; j < attach.size(); ++j) {
        const double ddx = (attach[j].x + 1.0) - lx;
        const double ddy = (attach[j].y + 1.0) - ly;
        const double d = std::sqrt(ddx * ddx + ddy * ddy);
        if (d < best) {
          best = d;
          r_idx = static_cast<int>(j);
        }
      }
      add_end_edge_idx2 = attach_idx[static_cast<size_t>(r_idx)];
      const Edge& re = edges[static_cast<size_t>(add_end_edge_idx2)];
      dist_2 = std::sqrt(std::pow(cur_x - re.x - 1.0, 2) + std::pow(cur_y - re.y - 1.0, 2));
      const double c_len_attach = static_cast<double>(attach_idx.size());
      if (dist_2 > dist_scale * c_len_attach) add_end_edge_idx2 = -1;
      break;
    }

    int c_type = 0;
    int add_end_edge_idx = -1;
    if (add_end_edge_idx1 >= 0 && add_end_edge_idx2 < 0) {
      c_type = 1;
      add_end_edge_idx = add_end_edge_idx1;
    } else if (add_end_edge_idx1 < 0 && add_end_edge_idx2 >= 0) {
      c_type = 2;
      add_end_edge_idx = add_end_edge_idx2;
    } else if (add_end_edge_idx1 >= 0 && add_end_edge_idx2 >= 0) {
      if (dist_1 < dist_2) {
        c_type = 1;
        add_end_edge_idx = add_end_edge_idx1;
      } else {
        c_type = 2;
        add_end_edge_idx = add_end_edge_idx2;
      }
    }

    if (c_type == 1) {
      const Edge& ae = edges[static_cast<size_t>(add_end_edge_idx)];
      if (is_start) {
        cur_c.insert(cur_c.begin(), ae);
        cfrags_idx[static_cast<size_t>(cid)].insert(cfrags_idx[static_cast<size_t>(cid)].begin(),
                                                    add_end_edge_idx);
        G.facs[static_cast<size_t>(cid)].nbrs_var[0] = reach_vids[0];
      } else {
        cur_c.push_back(ae);
        cfrags_idx[static_cast<size_t>(cid)].push_back(add_end_edge_idx);
        G.facs[static_cast<size_t>(cid)].nbrs_var[1] = reach_vids[0];
      }
      G.vars[static_cast<size_t>(reach_vids[0])].nbrs_fac.push_back(cid);
      auto& nf = var.nbrs_fac;
      nf.erase(std::remove(nf.begin(), nf.end(), cid), nf.end());
      edges[static_cast<size_t>(add_end_edge_idx)].d2f = 2.0;
      if (pass == GeomPass::First) end_point_eid_map[idx2(cur_y - 1, cur_x - 1, w)] = 0;
    } else if (c_type == 2) {
      const Edge& ae = edges[static_cast<size_t>(add_end_edge_idx)];
      if (is_start) {
        cur_c.insert(cur_c.begin(), ae);
        cfrags_idx[static_cast<size_t>(cid)].insert(cfrags_idx[static_cast<size_t>(cid)].begin(),
                                                    add_end_edge_idx);
      } else {
        cur_c.push_back(ae);
        cfrags_idx[static_cast<size_t>(cid)].push_back(add_end_edge_idx);
      }
      edges[static_cast<size_t>(add_end_edge_idx)].d2f = 3.0;
      const int jx = clampi(static_cast<int>(std::lround(ae.x)) + 1, 1, w) - 1;
      const int jy = clampi(static_cast<int>(std::lround(ae.y)) + 1, 1, h) - 1;
      jct_map[idx2(jy, jx, w)] = 1;
      if (pass == GeomPass::First) end_point_eid_map[idx2(cur_y - 1, cur_x - 1, w)] = 0;
    }
  }
}

std::vector<double> crop_colmaj(const std::vector<double>& src, int H, int W, int y0, int x0, int h_ref,
                                int w_ref) {
  std::vector<double> out(static_cast<size_t>(h_ref) * static_cast<size_t>(w_ref), 0.0);
  for (int x = 0; x < w_ref; ++x) {
    for (int y = 0; y < h_ref; ++y) {
      out[static_cast<size_t>(x * h_ref + y)] = src[idx2(y0 + y, x0 + x, W)];
    }
  }
  (void)H;
  return out;
}

std::vector<double> crop_int_as_double(const std::vector<int>& src, int W, int y0, int x0, int h_ref,
                                       int w_ref) {
  std::vector<double> out(static_cast<size_t>(h_ref) * static_cast<size_t>(w_ref), 0.0);
  for (int x = 0; x < w_ref; ++x) {
    for (int y = 0; y < h_ref; ++y) {
      out[static_cast<size_t>(x * h_ref + y)] =
          static_cast<double>(src[idx2(y0 + y, x0 + x, W)] > 0 ? 1 : 0);
    }
  }
  return out;
}

}  // namespace

GapFillResult contour_fill_gaps_DP(const std::vector<Contour>& cfrags_in,
                                   const std::vector<ContourEdgeIndices>& cfrags_idx_in, int h, int w,
                                   const std::vector<Edge>& edges_in, const GradientMaps& grad,
                                   const GapFillParams& params) {
  GapFillResult result;
  result.contours = cfrags_in;
  result.contour_edge_idx = cfrags_idx_in;
  result.edges = edges_in;

  auto& cfrags = result.contours;
  auto& cfrags_idx = result.contour_edge_idx;
  auto& edges = result.edges;

  const int r_range = params.DP_gap_range;
  const double theta_range = params.DP_angle_th;
  const double contrast_th = params.DP_contrast_th;
  const int search_nbr_range = params.shape_gap_range;
  const double search_ori_range = params.shape_ori_range;
  const double cost_th = 2.0;

  FacGraph G = construct_fac_graph(cfrags_idx);
  std::vector<int> end_point_eid_map, end_point_vid_map;
  build_endpoint_maps(G, edges, h, w, end_point_eid_map, end_point_vid_map);

  std::vector<int> EdgeGroupMap = convert_cfrags_to_EdgeGroupMap(cfrags, h, w);

  cv::Mat edgeim = cv::Mat::zeros(h, w, CV_8UC1);
  for (int y = 0; y < h; ++y)
    for (int x = 0; x < w; ++x)
      if (EdgeGroupMap[idx2(y, x, w)]) edgeim.at<uchar>(y, x) = 1;
  edgeim = morph_clean(edgeim);
  edgeim = morph_skel(edgeim);
  std::vector<int> jct_map = find_junction_map(edgeim, h, w);

  geometric_completion_pass(G, cfrags, cfrags_idx, edges, EdgeGroupMap, jct_map, end_point_eid_map,
                            end_point_vid_map, h, w, search_nbr_range, search_ori_range, GeomPass::First);

  // DT_map = exp(-bwdist.^2 / 8)
  cv::Mat bin255 = cv::Mat::zeros(h, w, CV_8UC1);
  for (int y = 0; y < h; ++y)
    for (int x = 0; x < w; ++x)
      if (EdgeGroupMap[idx2(y, x, w)]) bin255.at<uchar>(y, x) = 255;
  cv::Mat dt;
  cv::distanceTransform(bin255, dt, cv::DIST_L2, cv::DIST_MASK_PRECISE);
  std::vector<double> DT_map(static_cast<size_t>(h) * static_cast<size_t>(w), 0.0);
  for (int y = 0; y < h; ++y) {
    for (int x = 0; x < w; ++x) {
      const double d = dt.at<float>(y, x);
      DT_map[idx2(y, x, w)] = std::exp(-(d * d) / 8.0);
    }
  }

  for (size_t vid = 0; vid < G.vars.size(); ++vid) {
    FacVar& var = G.vars[vid];
    if (var.nbrs_fac.size() != 1) {
      edges[static_cast<size_t>(var.actual_edge_id)].d2f = static_cast<double>(var.nbrs_fac.size());
      continue;
    }
    if (edges[static_cast<size_t>(var.actual_edge_id)].d2f > 1.0) continue;

    int cur_x = clampi(static_cast<int>(std::lround(edges[static_cast<size_t>(var.actual_edge_id)].x)) + 1, 1, w);
    int cur_y = clampi(static_cast<int>(std::lround(edges[static_cast<size_t>(var.actual_edge_id)].y)) + 1, 1, h);
    if (jct_map[idx2(cur_y - 1, cur_x - 1, w)]) {
      edges[static_cast<size_t>(var.actual_edge_id)].d2f = 3.0;
      continue;
    }
    if (end_point_eid_map[idx2(cur_y - 1, cur_x - 1, w)] == 0) continue;

    const int cid = var.nbrs_fac[0];
    const Contour& cur_c = cfrags[static_cast<size_t>(cid)];
    const int e_size = static_cast<int>(cur_c.size());
    bool is_start = true;
    double end_dx = 0, end_dy = 0;
    if (G.facs[static_cast<size_t>(cid)].nbrs_var[0] == static_cast<int>(vid)) {
      const int j = std::min(e_size, 4) - 1;
      end_dx = cur_c[0].x - cur_c[static_cast<size_t>(j)].x;
      end_dy = cur_c[0].y - cur_c[static_cast<size_t>(j)].y;
    } else {
      is_start = false;
      // MATLAB: cur_c(end,:) - cur_c(max(e_size-3,1),:)
      const int j = std::max(e_size - 4, 0);
      end_dx = cur_c.back().x - cur_c[static_cast<size_t>(j)].x;
      end_dy = cur_c.back().y - cur_c[static_cast<size_t>(j)].y;
    }

    const int x_min = std::max(cur_x - r_range, 1);
    const int x_max = std::min(cur_x + r_range, w);
    const int y_min = std::max(cur_y - r_range, 1);
    const int y_max = std::min(cur_y + r_range, h);
    const int h_ref = y_max - y_min + 1;
    const int w_ref = x_max - x_min + 1;
    const double start_x_ref = static_cast<double>(cur_x - x_min);  // 0-based in crop
    const double start_y_ref = static_cast<double>(cur_y - y_min);
    const double start_theta = std::atan2(end_dy, end_dx);

    std::vector<double> E_col = crop_int_as_double(EdgeGroupMap, w, y_min - 1, x_min - 1, h_ref, w_ref);
    std::vector<double> soft_col = crop_colmaj(grad.edgemap_soft, h, w, y_min - 1, x_min - 1, h_ref, w_ref);
    std::vector<double> O_col = crop_colmaj(grad.thetamap, h, w, y_min - 1, x_min - 1, h_ref, w_ref);
    std::vector<double> DT_col = crop_colmaj(DT_map, h, w, y_min - 1, x_min - 1, h_ref, w_ref);
    for (size_t i = 0; i < E_col.size(); ++i) {
      if (E_col[i] > 0) DT_col[i] = 1.0;
    }

    DpGapResult dp = dp_gap_cpt(DT_col, soft_col, O_col, start_x_ref, start_y_ref, start_theta, h_ref, w_ref,
                                theta_range, contrast_th);

    const int N = h_ref * w_ref;
    double best_s = cost_th + 1.0;
    int best_i = -1;
    for (int i = 0; i < N; ++i) {
      double cn = dp.cost[static_cast<size_t>(i)] / std::max(dp.len[static_cast<size_t>(i)], 1e-12);
      if (dp.len[static_cast<size_t>(i)] > r_range) cn = 1000.0;
      const int yy = i % h_ref;
      const int xx = i / h_ref;
      if (jct_map[idx2(y_min - 1 + yy, x_min - 1 + xx, w)] > 0) cn = 1000.0;
      if (cn < best_s) {
        best_s = cn;
        best_i = i;
      }
    }
    if (best_i < 0 || best_s > cost_th) continue;

    // Backtrack (0-based linear). MATLAB while(back_p_mat(cur_ind)) stops on 0.
    std::vector<int> opt_path;
    int cur_ind = best_i;
    opt_path.push_back(cur_ind);
    while (dp.back_p[static_cast<size_t>(cur_ind)] != 0.0) {
      cur_ind = static_cast<int>(dp.back_p[static_cast<size_t>(cur_ind)]);
      opt_path.push_back(cur_ind);
    }

    std::vector<int> opt_x, opt_y;
    opt_x.reserve(opt_path.size());
    opt_y.reserve(opt_path.size());
    for (int lin : opt_path) {
      const int yy = lin % h_ref;
      const int xx = lin / h_ref;
      opt_y.push_back(yy + y_min);  // 1-based image
      opt_x.push_back(xx + x_min);
    }

    std::vector<double> path_prob;
    if (opt_x.size() >= 2) {
      path_prob.resize(opt_x.size() - 1);
      for (size_t i = 0; i + 1 < opt_x.size(); ++i) {
        path_prob[i] = grad.edgemap_soft[idx2(opt_y[i] - 1, opt_x[i] - 1, w)];
      }
    }
    if (!path_prob.empty()) {
      double mean_p = 0;
      for (double p : path_prob) mean_p += p;
      mean_p /= static_cast<double>(path_prob.size());
      if (mean_p < 0.05 && path_prob.size() > 5) continue;
    }

    int replace_end_edge_idx = -1;
    const int path_end_x = opt_x[0];
    const int path_end_y = opt_y[0];
    if (EdgeGroupMap[idx2(path_end_y - 1, path_end_x - 1, w)] > 0) {
      const int reach_cid0 = EdgeGroupMap[idx2(path_end_y - 1, path_end_x - 1, w)] - 1;
      const Contour& attach = cfrags[static_cast<size_t>(reach_cid0)];
      const auto& attach_idx = cfrags_idx[static_cast<size_t>(reach_cid0)];
      double best = std::numeric_limits<double>::infinity();
      int r_idx = 0;
      for (size_t j = 0; j < attach.size(); ++j) {
        const double ddx = (attach[j].x + 1.0) - path_end_x;
        const double ddy = (attach[j].y + 1.0) - path_end_y;
        const double d2 = ddx * ddx + ddy * ddy;
        if (d2 < best) {
          best = d2;
          r_idx = static_cast<int>(j);
        }
      }
      opt_x[0] = static_cast<int>(std::lround(attach[static_cast<size_t>(r_idx)].x + 1.0));
      opt_y[0] = static_cast<int>(std::lround(attach[static_cast<size_t>(r_idx)].y + 1.0));
      replace_end_edge_idx = attach_idx[static_cast<size_t>(r_idx)];
      if (r_idx != 0 && r_idx != static_cast<int>(attach_idx.size()) - 1) {
        edges[static_cast<size_t>(replace_end_edge_idx)].d2f = 3.0;
        jct_map[idx2(path_end_y - 1, path_end_x - 1, w)] = 1;
      } else {
        edges[static_cast<size_t>(replace_end_edge_idx)].d2f = 2.0;
      }
    }

    end_point_eid_map[idx2(cur_y - 1, cur_x - 1, w)] = 0;
    const int x2update = static_cast<int>(std::lround(static_cast<double>(opt_x[0])));
    const int y2update = static_cast<int>(std::lround(static_cast<double>(opt_y[0])));
    int eid2update = 0;  // 1-based like MATLAB map; 0 = none

    if (is_start) {
      std::vector<Edge> add_edges;
      for (size_t i = 0; i + 1 < opt_x.size(); ++i) {
        Edge e;
        e.x = opt_x[i] - 1.0;
        e.y = opt_y[i] - 1.0;
        e.dir = std::atan2(static_cast<double>(opt_y[i + 1] - opt_y[i]),
                           static_cast<double>(opt_x[i + 1] - opt_x[i]));
        e.conf = (i < path_prob.size()) ? path_prob[i] : 0.0;
        e.d2f = 0.0;
        add_edges.push_back(e);
      }
      Contour cfrag_new = add_edges;
      cfrag_new.push_back(cfrags[static_cast<size_t>(cid)].front());
      ContourEdgeIndices cfrag_idx_new;
      if (replace_end_edge_idx >= 0) {
        if (!add_edges.empty()) add_edges.erase(add_edges.begin());
        cfrag_new = add_edges;
        cfrag_new.insert(cfrag_new.begin(), edges[static_cast<size_t>(replace_end_edge_idx)]);
        cfrag_new.push_back(cfrags[static_cast<size_t>(cid)].front());
      }
      const int insert_begin = static_cast<int>(edges.size());
      edges.insert(edges.end(), add_edges.begin(), add_edges.end());
      const int insert_end = static_cast<int>(edges.size());
      for (int e = insert_begin; e < insert_end; ++e) cfrag_idx_new.push_back(e);
      cfrag_idx_new.push_back(cfrags_idx[static_cast<size_t>(cid)].front());
      eid2update = insert_begin + 1;  // MATLAB 1-based
      if (replace_end_edge_idx >= 0) {
        cfrag_idx_new.insert(cfrag_idx_new.begin(), replace_end_edge_idx);
        eid2update = 0;
      }
      cfrags.push_back(std::move(cfrag_new));
      cfrags_idx.push_back(std::move(cfrag_idx_new));
    } else {
      std::vector<int> ox = opt_x, oy = opt_y;
      std::vector<double> pp = path_prob;
      std::reverse(ox.begin(), ox.end());
      std::reverse(oy.begin(), oy.end());
      std::reverse(pp.begin(), pp.end());
      std::vector<Edge> add_edges;
      for (size_t i = 1; i < ox.size(); ++i) {
        Edge e;
        e.x = ox[i] - 1.0;
        e.y = oy[i] - 1.0;
        e.dir = std::atan2(static_cast<double>(oy[i] - oy[i - 1]), static_cast<double>(ox[i] - ox[i - 1]));
        e.conf = (i - 1 < pp.size()) ? pp[i - 1] : 0.0;
        e.d2f = 0.0;
        add_edges.push_back(e);
      }
      Contour cfrag_new;
      cfrag_new.push_back(cfrags[static_cast<size_t>(cid)].back());
      cfrag_new.insert(cfrag_new.end(), add_edges.begin(), add_edges.end());
      if (replace_end_edge_idx >= 0) {
        if (!add_edges.empty()) add_edges.pop_back();
        cfrag_new.clear();
        cfrag_new.push_back(cfrags[static_cast<size_t>(cid)].back());
        cfrag_new.insert(cfrag_new.end(), add_edges.begin(), add_edges.end());
        cfrag_new.push_back(edges[static_cast<size_t>(replace_end_edge_idx)]);
      }
      const int insert_begin = static_cast<int>(edges.size());
      edges.insert(edges.end(), add_edges.begin(), add_edges.end());
      ContourEdgeIndices cfrag_idx_new;
      cfrag_idx_new.push_back(cfrags_idx[static_cast<size_t>(cid)].back());
      for (int e = insert_begin; e < static_cast<int>(edges.size()); ++e) cfrag_idx_new.push_back(e);
      eid2update = static_cast<int>(edges.size());  // MATLAB 1-based last
      if (replace_end_edge_idx >= 0) {
        cfrag_idx_new.push_back(replace_end_edge_idx);
        eid2update = 0;
      }
      cfrags.push_back(std::move(cfrag_new));
      cfrags_idx.push_back(std::move(cfrag_idx_new));
    }

    if (y2update >= 1 && y2update <= h && x2update >= 1 && x2update <= w)
      end_point_eid_map[idx2(y2update - 1, x2update - 1, w)] = eid2update;
    const int new_cid1 = static_cast<int>(cfrags.size());
    for (size_t i = 0; i < opt_x.size(); ++i) {
      const int yy = clampi(static_cast<int>(std::lround(static_cast<double>(opt_y[i]))), 1, h) - 1;
      const int xx = clampi(static_cast<int>(std::lround(static_cast<double>(opt_x[i]))), 1, w) - 1;
      EdgeGroupMap[idx2(yy, xx, w)] = new_cid1;
    }
  }

  // Second geometric pass after rebuilding graph
  G = construct_fac_graph(cfrags_idx);
  build_endpoint_maps(G, edges, h, w, end_point_eid_map, end_point_vid_map);
  geometric_completion_pass(G, cfrags, cfrags_idx, edges, EdgeGroupMap, jct_map, end_point_eid_map,
                            end_point_vid_map, h, w, search_nbr_range, search_ori_range, GeomPass::Second);

  // Remove empty
  std::vector<Contour> c2;
  std::vector<ContourEdgeIndices> i2;
  for (size_t i = 0; i < cfrags.size(); ++i) {
    if (cfrags_idx[i].empty()) continue;
    c2.push_back(std::move(cfrags[i]));
    i2.push_back(std::move(cfrags_idx[i]));
  }
  result.contours = std::move(c2);
  result.contour_edge_idx = std::move(i2);
  return result;
}

}  // namespace tcg
