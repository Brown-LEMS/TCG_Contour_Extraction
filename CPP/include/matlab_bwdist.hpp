#pragma once

#include <opencv2/core.hpp>
#include <string>
#include <vector>

namespace tcg {

/// Exact Euclidean distance transform matching MATLAB bwdist(BW,'euclidean').
/// `feature` is a 0/1 (or 0/nonzero) CV_8U image; output is CV_32F, same size.
/// Each pixel gets distance to the nearest nonzero pixel (0 on features).
/// If there are no features, all entries are +Inf (as in MATLAB).
cv::Mat bwdist_euclidean(const cv::Mat& feature);

}  // namespace tcg
