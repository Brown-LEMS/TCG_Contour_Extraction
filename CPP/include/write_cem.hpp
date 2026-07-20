#pragma once

#include "tcg_types.hpp"

#include <string>

namespace tcg {

/// Port of util/io/write_cem.m — writes .CEM v2.0 loadable by MATLAB load_contours.
/// Contour vertices are concatenated into the edgemap; contour lines store 0-based ids.
bool write_cem(const std::string& path, const std::vector<Contour>& contours, int h, int w,
               std::string& error);

}  // namespace tcg
