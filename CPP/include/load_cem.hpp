#pragma once

#include "tcg_types.hpp"

#include <string>

namespace tcg {

/// Load .CEM v2.0 (same structure as MATLAB load_contours).
bool load_cem(const std::string& path, CemFile& out, std::string& error);

}  // namespace tcg
