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

bool write_cemv(const std::string& path, const std::vector<Contour>& contours, std::string& error) {
  // Port of det_save_cemv.m (used by Main_cem2cemv.m).
  std::ofstream out(path);
  if (!out) {
    error = "Cannot open for write: " + path;
    return false;
  }

  out << "# CONTOUR_EDGE_MAP : Logical-Linear + Shock_Grouping\n";
  out << "# .cem files\n";
  out << "#\n";
  out << "# Format :\n";
  out << "# Each contour block will consist of the following\n";
  out << "# [BEGIN CONTOUR]\n";
  out << "# EDGE_COUNT=num_of_edges\n";
  out << "# [Pixel_Pos]  Pixel_Dir Pixel_Conf  [Sub_Pixel_Pos] Sub_Pixel_Dir Sub_Pixel_Conf\n";
  out << "# ...\n";
  out << "# ...\n";
  out << "# [END CONTOUR]\n";
  out << "\n";

  // MATLAB leaves these empty (det_save_cemv.m).
  out << "CONTOUR_COUNT=\n";
  out << "TOTAL_EDGE_COUNT=\n";

  out << std::fixed << std::setprecision(6);
  for (const auto& c : contours) {
    out << "[BEGIN CONTOUR]\n";
    out << "EDGE_COUNT=" << c.size() << "\n";
    for (const auto& e : c) {
      // Pixel fields unused (zeros); sub-pixel stores x, y, dir, conf.
      out << " [0, 0] 0.000000 0.000000 [" << e.x << ", " << e.y << "] " << e.dir << " " << e.conf
          << "\n";
    }
    out << "[END CONTOUR]\n\n";
  }

  return true;
}

}  // namespace tcg
