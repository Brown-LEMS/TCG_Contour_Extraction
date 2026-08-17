#include "image_gradient.hpp"

#include <algorithm>
#include <cmath>
#include <opencv2/imgcodecs.hpp>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace tcg {
namespace {

double wrap_to_pi(double a) {
  const double two_pi = 2.0 * M_PI;
  a = std::fmod(a + M_PI, two_pi);
  if (a < 0) a += two_pi;
  return a - M_PI;
}

// MATLAB rgb2gray on uint8: weighted sum then round back to uint8;
// imgradientxy then promotes with double(I).
cv::Mat rgb2gray_matlab(const cv::Mat& bgr_u8) {
  CV_Assert(bgr_u8.type() == CV_8UC3);
  cv::Mat gray(bgr_u8.rows, bgr_u8.cols, CV_64FC1);
  for (int y = 0; y < bgr_u8.rows; ++y) {
    const cv::Vec3b* row = bgr_u8.ptr<cv::Vec3b>(y);
    double* out = gray.ptr<double>(y);
    for (int x = 0; x < bgr_u8.cols; ++x) {
      const double b = row[x][0];
      const double g = row[x][1];
      const double r = row[x][2];
      double yv = 0.2989 * r + 0.5870 * g + 0.1140 * b;
      yv = std::round(yv);
      if (yv < 0.0) yv = 0.0;
      if (yv > 255.0) yv = 255.0;
      out[x] = yv;
    }
  }
  return gray;
}

inline double sample_replicate(const cv::Mat& gray /* CV_64F */, int y, int x) {
  y = std::max(0, std::min(gray.rows - 1, y));
  x = std::max(0, std::min(gray.cols - 1, x));
  return gray.at<double>(y, x);
}

// MATLAB imgradientxy(...,'sobel'):
//   h = -fspecial('sobel');
//   Gx = imfilter(I, h', 'replicate');
//   Gy = imfilter(I, h,  'replicate');
// where -fspecial('sobel') = [-1 -2 -1; 0 0 0; 1 2 1]
void sobel_matlab(const cv::Mat& gray, cv::Mat& gx, cv::Mat& gy) {
  const int h = gray.rows;
  const int w = gray.cols;
  gx.create(h, w, CV_64FC1);
  gy.create(h, w, CV_64FC1);

  // Gy kernel h = [-1 -2 -1; 0 0 0; 1 2 1]
  // Gx kernel h' = [-1 0 1; -2 0 2; -1 0 1]
  for (int y = 0; y < h; ++y) {
    double* gx_row = gx.ptr<double>(y);
    double* gy_row = gy.ptr<double>(y);
    for (int x = 0; x < w; ++x) {
      const double i00 = sample_replicate(gray, y - 1, x - 1);
      const double i01 = sample_replicate(gray, y - 1, x);
      const double i02 = sample_replicate(gray, y - 1, x + 1);
      const double i10 = sample_replicate(gray, y, x - 1);
      const double i12 = sample_replicate(gray, y, x + 1);
      const double i20 = sample_replicate(gray, y + 1, x - 1);
      const double i21 = sample_replicate(gray, y + 1, x);
      const double i22 = sample_replicate(gray, y + 1, x + 1);

      gx_row[x] = -i00 + i02 - 2.0 * i10 + 2.0 * i12 - i20 + i22;
      gy_row[x] = -i00 - 2.0 * i01 - i02 + i20 + 2.0 * i21 + i22;
    }
  }
}

}  // namespace

bool load_image_gradient(const std::string& image_path, GradientMaps& out, std::string& error) {
  cv::Mat bgr = cv::imread(image_path, cv::IMREAD_COLOR);
  if (bgr.empty()) {
    error = "Cannot read image: " + image_path;
    return false;
  }

  // Port of:
  //   [edgemap_soft, thetamap] = imgradient(rgb2gray(img));
  //   edgemap_soft = edgemap_soft/max(edgemap_soft(:));
  //   thetamap = wrapToPi(-thetamap/180*pi + pi/2);
  cv::Mat gray = rgb2gray_matlab(bgr);
  cv::Mat gx, gy;
  sobel_matlab(gray, gx, gy);

  const int h = gray.rows;
  const int w = gray.cols;
  out.height = h;
  out.width = w;
  out.edgemap_soft.assign(static_cast<size_t>(h) * static_cast<size_t>(w), 0.0);
  out.thetamap.assign(static_cast<size_t>(h) * static_cast<size_t>(w), 0.0);

  double max_mag = 0.0;
  for (int y = 0; y < h; ++y) {
    const double* gx_row = gx.ptr<double>(y);
    const double* gy_row = gy.ptr<double>(y);
    for (int x = 0; x < w; ++x) {
      const double gxv = gx_row[x];
      const double gyv = gy_row[x];
      const double mag = std::hypot(gxv, gyv);
      const size_t i = static_cast<size_t>(y) * static_cast<size_t>(w) + static_cast<size_t>(x);
      out.edgemap_soft[i] = mag;
      if (mag > max_mag) max_mag = mag;
      // MATLAB: Gdir = atan2(-Gy, Gx) * 180/pi  (degrees)
      // then:   wrapToPi(-Gdir/180*pi + pi/2)
      const double gdir_rad = std::atan2(-gyv, gxv);
      out.thetamap[i] = wrap_to_pi(-gdir_rad + M_PI / 2.0);
    }
  }
  if (max_mag > 0) {
    for (double& v : out.edgemap_soft) v /= max_mag;
  }
  return true;
}

}  // namespace tcg
