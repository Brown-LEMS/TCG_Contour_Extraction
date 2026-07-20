#pragma once

#include <vector>

namespace tcg {

struct DpGapResult {
  std::vector<double> cost;     // column-major, size h*w
  std::vector<double> back_p;   // 0-based linear index x*h+y
  std::vector<double> len;
};

/// Port of util/DP_gap_cpt.cpp (MEX stripped).
/// Maps are column-major (MATLAB layout): index = x*h + y.
/// start_xy_theta = {x, y, theta} in 0-based crop coordinates.
DpGapResult dp_gap_cpt(const std::vector<double>& E_map, const std::vector<double>& E_map_soft,
                       const std::vector<double>& O_map, double start_x, double start_y, double start_theta,
                       int h, int w, double theta_th, double contrast_th);

}  // namespace tcg
