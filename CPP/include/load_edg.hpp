#pragma once

#include "tcg_types.hpp"

#include <string>

namespace tcg {

/// Load an EDGE_MAP .edg file (v1 / v2 / v3 header). Fills edgemap and thetamap like MATLAB load_edg opt==1.
bool load_edg(const std::string& path, EdgFile& out, std::string& error);

}  // namespace tcg
