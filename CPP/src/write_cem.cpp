#include "write_cem.hpp"

#include <fstream>
#include <iomanip>

namespace tcg {

bool write_cem(const std::string& path, const std::vector<Contour>& contours, int h, int w,
               std::string& error) {
  std::ofstream out(path);
  if (!out) {
    error = "Cannot open for write: " + path;
    return false;
  }

  out << std::setprecision(10);
  out << ".CEM v2.0\n";
  out << "size=[" << w << " " << h << "]\n";
  out << "[Edgemap]\n";

  size_t total_edges = 0;
  for (const auto& c : contours) total_edges += c.size();
  out << "count=" << total_edges << "\n";

  for (const auto& c : contours) {
    for (const auto& e : c) {
      out << "(" << e.x << ", " << e.y << ")\t" << e.dir << "\t" << e.conf << "\t" << e.d2f << "\n";
    }
  }

  out << "[Contours]\n";
  out << "count=" << contours.size() << "\n";
  int next_id = 0;
  for (const auto& c : contours) {
    out << "[";
    for (size_t k = 0; k < c.size(); ++k) {
      if (k) out << " ";
      out << next_id++;
    }
    out << " ]\n";
  }

  out << "[Contour Properties]\n";
  out << "# <len> <avg. str> <mean con> <Lstd> <Rstd> <avg. d2f> <avg. k> <max k>\n";
  for (const auto& c : contours) {
    out << c.size() << " 2 0 0 0 0 0.5 0.5\n";
  }

  return true;
}

}  // namespace tcg
