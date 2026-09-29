#include <cmath>

#pragma once

namespace tcg {

/// Every tunable parameter and decision threshold in the C++ TCG pipeline.
//> Defaults match the original MATLAB code

struct TcgParams {

  bool vis{false};

  // ========================= Corner break (contour_breaker_at_corner) =========================
  int nbr_num_edges{20};                      //> Neighborhood length used when measuring orientation change.
  double corner_angle_th{M_PI / 6.0};         //> Default corner threshold (radians). Used on the second corner-break pass.
  double corner_ori_diff_th{M_PI / 18.0};     //> First corner-break pass only. The original MATLAB code passes pi/18 here; later passes fall back to `corner_angle_th`.
  int corner_smooth_span{5};                  //> `smooth` default moving-average span.
  int corner_break_iters{1};                  //> Divides the co-circular neighborhood (`cur_nbr_num / iter`).
  double corner_nbr_divisor{2.0};             //> `ceil(nbr_num_edges / corner_nbr_divisor)` in the neighborhood size.
  double corner_len_divisor{6.0};             //> `ceil(contour_size / corner_len_divisor)` in the neighborhood size.
  int corner_closed_len_factor{2};            //> Closed fragments shorter than `factor * nbr_num_edges` are left intact.
  double corner_dtheta_factor{2.0};           //> Locally smooth samples (`Dtheta < ori_diff_th / factor`) are not corners.

  // ========================= Gap fill (contour_fill_gaps_DP) ==========================
  int DP_gap_range{15};                        //> Gap fill range (in number of edgels).
  double DP_angle_th{M_PI / 4.0};              //> Angle threshold (in radians)
  double DP_contrast_th{0.1};                  //> Contrast threshold (in [0, 1])
  int shape_gap_range{8};                      //> The range, in number of edgels, to fill gaps with shape-based DP
  double shape_ori_range{M_PI / 9.0};          //> Shape orientation range (radians).
  double gap_cost_th{2.0};                     //> Reject a DP gap whose normalized cost exceeds this.
  double gap_path_mean_prob_th{0.05};          //> Reject a completed gap whose mean edge probability is below this when the path is longer than `gap_path_prob_min_len`.
  int gap_path_prob_min_len{5};                //> Minimum length of a gap path (in number of edgels).
  double gap_invalid_cost{1000.0};             //> Cost assigned to DP cells past the gap range or on a junction.
  double dt_sq_denom{8.0};                     //> `DT_map = exp(-bwdist^2 / dt_sq_denom)`. MATLAB uses `2^2/2` (= 8).

  // ========================= Geometric completion (contour_completion_geometric) ==========================
  double geom_line_scale_pass1{0.5};           //> Extending-line length, as a fraction of `shape_gap_range`, per geometric pass.
  double geom_line_scale_pass2{0.75};          //> Extending-line length, as a fraction of `shape_gap_range`, per geometric pass.
  double geom_dist_scale_pass1{1.0};           //> Max link distance, as a multiple of the fragment length, per geometric pass.
  double geom_dist_scale_pass2{2.0};           //> Max link distance, as a multiple of the fragment length, per geometric pass.
  int geom_min_edges_pass1{10};                //> First geometric pass skips fragments shorter than this (in edgels).
  int geom_tangent_edges{5};                   //> Edgels used to estimate the endpoint tangent in geometric completion.
  int dp_tangent_edges{4};                     //> Edgels used to estimate the endpoint tangent before the DP search.
  double geom_wide_fan_half_angle{M_PI / 3.0}; //> Extra search fan: `angle +/- geom_wide_fan_half_angle`, radius `geom_wide_fan_radius`.
  double geom_wide_fan_radius{2.0};            //> Extra search fan: `angle +/- geom_wide_fan_half_angle`, radius `geom_wide_fan_radius`.
  double edgegroup_sample_ratio{2.0};          //> Samples per unit length when rasterizing fragments onto the edge-group map.

  // ========================= DP step cost ==========================
  double dp_weight_gradient{0.65};           //> Weight on gradient magnitude.
  double dp_weight_shape{0.35};              //> Weight on the shape term.
  double dp_grad_denom_eps{0.01};            //> Soft gradient response `fG = fG^2 / (fG^2 + dp_grad_denom_eps)`.
  double dp_ori_kappa{4.0};                  //> Orientation agreement `fO = exp(-acos(dq)^2 * kappa/pi * kappa/pi / 2)`.
  double dp_local_dir_cos_min{0.1};          //> Skip a neighbor step whose direction cosine with the current heading is below this.
  double dp_start_shape_weight_scale{2.0};   //> Shape-term weight multiplier on the first step away from the start pixel.

  // ========================= Prune noisy curves ==========================
  double noise_len_th{5.0};                  //> Length threshold for pruning noisy curves.
  double noise_prob_th{0.05};                //> Probability threshold for pruning noisy curves.
  double prune_branch_len_factor{2.0};       //> Degree-3 branches longer than `factor * noise_len_th` are kept.
  double prune_isolated_len_factor{3.0};     //> Isolated or loop fragments use `factor * noise_len_th` as the length gate.

  // ========================= Geometric curve merging ==========================
  double geom_merge_angle_th{M_PI / 6.0};    //> Angle threshold for geometric merge.

  // ========================= Junction classification ==========================
  double BP_merge_angle_th{M_PI / 9.0};
  int BP_nbr_num_edges{20};                  //> Neighborhood length used when measuring orientation change.
  double BP_clen_th{15.0};                  //> Length threshold for junction classification.
  double BP_cost_weight{1.0};                //> Preference on the co-circular cost (`w0`; 1 means no extra preference).

  // ========================= Co-circular cost and contour resampling ==========================
  int cocirc_local_len{5};                   //> Local tangent length (in edgels) before interpolation.
  int cocirc_interp_len{15};                 //> Local tangent length (in samples) after interpolation.
  int cfrag_resample_factor{2};              //> `interpolate_cfrag` sample count is `round(length) * cfrag_resample_factor`.
};

inline constexpr TcgParams kTcgParams{};

}  // namespace tcg
