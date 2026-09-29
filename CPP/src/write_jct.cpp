#include "write_jct.hpp"

#include <fstream>
#include <iomanip>

namespace tcg {

bool write_jct(const std::string& path, const std::vector<Edge>& T_junctions,
               const std::vector<Edge>& Y_junctions, int h, int w, std::string& error) {
  std::ofstream out(path);
  if (!out) {
    error = "Cannot open for write: " + path;
    return false;
  }

  out << std::setprecision(10);
  out << "# TCG Junctions v1.0\n";
  out << "# Coordinates are 0-based (same as .cem). Add +1 for MATLAB image overlays.\n";
  out << "# T: degree-3 nodes where two branches were merged (classify-BP).\n";
  out << "# Y: degree>=3 graph nodes remaining after classify-BP / on final graph.\n";
  out << "size=[" << w << " " << h << "]\n";

  out << "[T_junctions]\n";
  out << "count=" << T_junctions.size() << "\n";
  for (const auto& e : T_junctions) {
    out << e.x << " " << e.y << " " << e.dir << " " << e.conf << "\n";
  }

  out << "[Y_junctions]\n";
  out << "count=" << Y_junctions.size() << "\n";
  for (const auto& e : Y_junctions) {
    out << e.x << " " << e.y << " " << e.dir << " " << e.conf << "\n";
  }

  return true;
}

}  // namespace tcg
