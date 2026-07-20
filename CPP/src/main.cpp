#include "contour_breaker.hpp"
#include "contour_fill_gaps.hpp"
#include "image_gradient.hpp"
#include "load_cem.hpp"
#include "load_edg.hpp"

#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

static void print_usage(const char* argv0) {
  std::cerr << "Usage: " << argv0
            << " <input.edg> <input.cem> <input.image> [ori_diff_th] [corners_out.txt]\n"
            << "  Runs Step 3 (corner break) then Step 4 (gap fill DP).\n"
            << "  ori_diff_th: radians (default: pi/18)\n";
}

int main(int argc, char** argv) {
  if (argc < 4) {
    print_usage(argv[0]);
    return 1;
  }

  const std::string edg_path = argv[1];
  const std::string cem_path = argv[2];
  const std::string img_path = argv[3];
  double ori_th = M_PI / 18.0;
  if (argc >= 5) ori_th = std::strtod(argv[4], nullptr);

  tcg::EdgFile edg;
  tcg::CemFile cem;
  std::string err;

  if (!tcg::load_edg(edg_path, edg, err)) {
    std::cerr << "load_edg: " << err << "\n";
    return 1;
  }
  if (!tcg::load_cem(cem_path, cem, err)) {
    std::cerr << "load_cem: " << err << "\n";
    return 1;
  }

  tcg::GradientMaps grad;
  if (!tcg::load_image_gradient(img_path, grad, err)) {
    std::cerr << "load_image_gradient: " << err << "\n";
    return 1;
  }

  std::cout << "[Summary] EDG " << edg.width << "x" << edg.height << " edges = " << edg.edges.size()
            << " | CEM contours = " << cem.contours.size() << " | IMG size = " << grad.width << "x" << grad.height
            << "\n";

  tcg::BreakerParams bparams;
  bparams.nbr_num_edges = 20;
  bparams.corner_angle_th = M_PI / 6.0;

  //>========================== MATLAB contour_breaker_at_corner function ==========================
  auto t0 = std::chrono::steady_clock::now();
  tcg::CornerBreakResult broken = tcg::contour_breaker_at_corner(cem.contours, cem.contour_edge_idx, bparams, ori_th);
  auto t1 = std::chrono::steady_clock::now();
  const double contour_breaker_time_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
  std::cout << "[Step 3] number of curve fragments = " << broken.contours.size()
            << " / number of corners = " << broken.corner_pts.size()
            << " / time = " << contour_breaker_time_ms << " ms" << std::endl;
  //>========================== MATLAB contour_breaker_at_corner function ==========================

  //> Optional: evaluate the result consistency between MATLAB and C++ on Step 3
  // if (argc >= 6) {
  //   std::ofstream ofs(argv[5]);
  //   if (!ofs) {
  //     std::cerr << "Cannot write corner file: " << argv[5] << "\n";
  //     return 1;
  //   }
  //   ofs << std::setprecision(17);
  //   for (const auto& e : broken.corner_pts) {
  //     ofs << e.x << '\t' << e.y << '\t' << e.dir << '\t' << e.conf << '\t' << e.d2f << '\n';
  //   }
  // }

  //>========================== MATLAB contour_fill_gaps_DP function ==========================
  tcg::GapFillParams gparams;
  gparams.DP_gap_range = 15;
  gparams.DP_angle_th = M_PI / 4.0;
  gparams.DP_contrast_th = 0.1;
  gparams.shape_gap_range = 8;
  gparams.shape_ori_range = M_PI / 9.0;

  const int h = grad.height;
  const int w = grad.width;

  t0 = std::chrono::steady_clock::now();
  tcg::GapFillResult filled = tcg::contour_fill_gaps_DP(broken.contours, broken.contour_edge_idx, h, w, cem.edges, grad, gparams);
  t1 = std::chrono::steady_clock::now();
  const double contour_fill_gaps_time_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();

  std::cout << "[Step 4] number of curve fragments = " << filled.contours.size() 
            << " / number of edges = " << filled.edges.size()
            << " / time = " << contour_fill_gaps_time_ms << " ms" << std::endl;
  //>========================== MATLAB contour_fill_gaps_DP function ==========================

  return 0;
}
