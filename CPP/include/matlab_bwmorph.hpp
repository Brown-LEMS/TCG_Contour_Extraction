#pragma once

#include <opencv2/core.hpp>

namespace tcg {

/// MATLAB bwmorph(BW,'skel',Inf) using images.internal.lutskel LUTs + bwlookup.
/// Input/output are 0/1 binary images (CV_8U). Outside pixels are treated as 0.
cv::Mat morph_skel_matlab(const cv::Mat& bin);

}  // namespace tcg
