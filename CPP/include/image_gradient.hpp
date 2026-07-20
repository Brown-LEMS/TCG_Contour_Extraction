#pragma once

#include "tcg_types.hpp"

#include <string>

namespace tcg {

/// Port of MATLAB: [edgemap_soft, thetamap] = imgradient(rgb2gray(img));
/// then normalize magnitude and wrapToPi(-deg2rad(angle)+pi/2).
bool load_image_gradient(const std::string& image_path, GradientMaps& out, std::string& error);

}  // namespace tcg
