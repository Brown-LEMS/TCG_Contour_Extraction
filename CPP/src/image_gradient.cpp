#include "image_gradient.hpp"

#include <cmath>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/imgproc.hpp>

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

}  // namespace

bool load_image_gradient(const std::string& image_path, GradientMaps& out, std::string& error) {
  cv::Mat bgr = cv::imread(image_path, cv::IMREAD_COLOR);
  if (bgr.empty()) {
    error = "Cannot read image: " + image_path;
    return false;
  }
  cv::Mat gray;
  cv::cvtColor(bgr, gray, cv::COLOR_BGR2GRAY);
  gray.convertTo(gray, CV_64F);

  cv::Mat gx, gy;
  cv::Sobel(gray, gx, CV_64F, 1, 0, 3);
  cv::Sobel(gray, gy, CV_64F, 0, 1, 3);

  const int h = gray.rows;
  const int w = gray.cols;
  out.height = h;
  out.width = w;
  out.edgemap_soft.assign(static_cast<size_t>(h) * static_cast<size_t>(w), 0.0);
  out.thetamap.assign(static_cast<size_t>(h) * static_cast<size_t>(w), 0.0);

  double max_mag = 0.0;
  for (int y = 0; y < h; ++y) {
    for (int x = 0; x < w; ++x) {
      const double dx = gx.at<double>(y, x);
      const double dy = gy.at<double>(y, x);
      const double mag = std::sqrt(dx * dx + dy * dy);
      out.edgemap_soft[static_cast<size_t>(y) * static_cast<size_t>(w) + static_cast<size_t>(x)] = mag;
      if (mag > max_mag) max_mag = mag;
      // MATLAB imgradient angle (degrees) -> radians, then wrapToPi(-th + pi/2)
      const double ang = std::atan2(dy, dx);
      out.thetamap[static_cast<size_t>(y) * static_cast<size_t>(w) + static_cast<size_t>(x)] =
          wrap_to_pi(-ang + M_PI / 2.0);
    }
  }
  if (max_mag > 0) {
    for (double& v : out.edgemap_soft) v /= max_mag;
  }
  return true;
}

}  // namespace tcg
