#include "contour_breaker.hpp"
#include "load_cem.hpp"
#include "load_edg.hpp"

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <iomanip>

static void print_usage(const char* argv0) {
  std::cerr << "Usage: " << argv0 << " <input.edg> <input.cem> [ori_diff_th] [corners_out.txt]\n"
            << "  ori_diff_th: corner angle threshold in radians (default: pi/18, matching main_TCG_CH.m)\n"
            << "  corners_out.txt: optional path to write detected corners (x y dir conf d2f)\n";
}

int main(int argc, char** argv) {

  if (argc < 3) {
    print_usage(argv[0]);
    return 1;
  }

  const std::string edg_path = argv[1];
  const std::string cem_path = argv[2];
  double ori_th = M_PI / 18.0;
  if (argc >= 4) ori_th = std::strtod(argv[3], nullptr);

  tcg::EdgFile edg;
  tcg::CemFile cem;
  std::string err;

  if (!tcg::load_edg(edg_path, edg, err)) {
    std::cerr << "load_edg: " << err << std::endl;
    return 1;
  }
  if (!tcg::load_cem(cem_path, cem, err)) {
    std::cerr << "load_cem: " << err << std::endl;
    return 1;
  }

  //> Print summary of edg and cem files
  std::cout << "[Summary] EDG and CEM files" << std::endl;
  std::cout << " - EDG: size = " << edg.width << "x" << edg.height
            << " / number of edges = " << edg.edges.size() << std::endl;
  std::cout << " - CEM: image size = " << cem.width << "x" << cem.height 
            << " / number of edges = " << cem.edges.size()
            << " / number of contours = " << cem.contours.size() << std::endl;

  if (edg.edges.size() != cem.edges.size()) {
    std::cout << "Note: EDG and CEM edge counts differ (expected if files are not paired).\n";
  }

  tcg::BreakerParams params;
  params.nbr_num_edges = 20;
  params.corner_angle_th = M_PI / 6.0;

  //>========================== MATLAB contour_breaker_at_corner function ==========================
  tcg::CornerBreakResult out = tcg::contour_breaker_at_corner(cem.contours, cem.contour_edge_idx, params, ori_th);

  std::cout << "[contour_breaker_at_corner] number of curve fragments = " << out.contours.size()
            << " / number of corners detected = " << out.corner_pts.size() << "\n";
  // std::cout << "ori_diff_th=" << ori_th << " rad\n";

  // Evaluation (Optional): write corners for MATLAB comparison (x y dir conf d2f).
  if (argc >= 5) {
    std::ofstream ofs(argv[4]);
    if (!ofs) {
      std::cerr << "Cannot write corner file: " << argv[4] << "\n";
      return 1;
    }
    ofs << std::setprecision(17);
    for (const auto& e : out.corner_pts) {
      ofs << e.x << '\t' << e.y << '\t' << e.dir << '\t' << e.conf << '\t' << e.d2f << '\n';
    }
    std::cout << "Wrote " << out.corner_pts.size() << " corners to " << argv[4] << "\n";
  }
  //>========================== MATLAB contour_breaker_at_corner function ==========================

  return 0;
}
