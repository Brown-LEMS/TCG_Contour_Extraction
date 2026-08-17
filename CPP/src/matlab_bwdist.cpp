#include "matlab_bwdist.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace tcg {
namespace {

// Large finite "infinity" so intersection arithmetic stays well-defined.
constexpr float kInf = 1e10f;

// Felzenszwalb & Huttenlocher 1D squared-distance transform.
void edt_1d(const float* f, float* d, int n, std::vector<int>& v, std::vector<float>& z) {
  v.resize(static_cast<size_t>(n));
  z.resize(static_cast<size_t>(n) + 1);
  int k = 0;
  v[0] = 0;
  z[0] = -std::numeric_limits<float>::infinity();
  z[1] = std::numeric_limits<float>::infinity();
  for (int q = 1; q < n; ++q) {
    auto sep = [&](int q_, int r_) {
      return ((f[q_] + static_cast<float>(q_) * static_cast<float>(q_)) -
              (f[r_] + static_cast<float>(r_) * static_cast<float>(r_))) /
             (2.f * static_cast<float>(q_ - r_));
    };
    float s = sep(q, v[static_cast<size_t>(k)]);
    while (s <= z[static_cast<size_t>(k)]) {
      --k;
      s = sep(q, v[static_cast<size_t>(k)]);
    }
    ++k;
    v[static_cast<size_t>(k)] = q;
    z[static_cast<size_t>(k)] = s;
    z[static_cast<size_t>(k) + 1] = std::numeric_limits<float>::infinity();
  }
  k = 0;
  for (int q = 0; q < n; ++q) {
    while (z[static_cast<size_t>(k) + 1] < static_cast<float>(q)) ++k;
    const int r = v[static_cast<size_t>(k)];
    const float dx = static_cast<float>(q - r);
    d[q] = dx * dx + f[r];
  }
}

}  // namespace

cv::Mat bwdist_euclidean(const cv::Mat& feature) {
  CV_Assert(feature.type() == CV_8UC1);
  const int h = feature.rows;
  const int w = feature.cols;
  cv::Mat out(h, w, CV_32FC1);

  std::vector<float> col_f(static_cast<size_t>(h));
  cv::Mat sq(h, w, CV_32FC1);

  bool any_feature = false;
  for (int x = 0; x < w; ++x) {
    float prev = kInf;
    for (int y = 0; y < h; ++y) {
      if (feature.at<uchar>(y, x)) {
        prev = 0.f;
        any_feature = true;
      } else if (prev < kInf) {
        prev += 1.f;
      }
      col_f[static_cast<size_t>(y)] = prev;
    }
    prev = kInf;
    for (int y = h - 1; y >= 0; --y) {
      if (feature.at<uchar>(y, x)) {
        prev = 0.f;
      } else if (prev < kInf) {
        prev += 1.f;
      }
      const float d = std::min(col_f[static_cast<size_t>(y)], prev);
      col_f[static_cast<size_t>(y)] = (d >= kInf * 0.5f) ? kInf : d * d;
      sq.at<float>(y, x) = col_f[static_cast<size_t>(y)];
    }
  }

  if (!any_feature) {
    out.setTo(std::numeric_limits<float>::infinity());
    return out;
  }

  std::vector<float> row_f(static_cast<size_t>(w));
  std::vector<float> row_d(static_cast<size_t>(w));
  std::vector<int> v;
  std::vector<float> z;
  for (int y = 0; y < h; ++y) {
    for (int x = 0; x < w; ++x) row_f[static_cast<size_t>(x)] = sq.at<float>(y, x);
    edt_1d(row_f.data(), row_d.data(), w, v, z);
    for (int x = 0; x < w; ++x) {
      out.at<float>(y, x) = std::sqrt(row_d[static_cast<size_t>(x)]);
    }
  }
  return out;
}

}  // namespace tcg
