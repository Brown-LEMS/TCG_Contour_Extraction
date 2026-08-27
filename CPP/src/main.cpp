#include "break_t_junctions.hpp"
#include "classify_bp.hpp"
#include "contour_breaker.hpp"
#include "contour_fill_gaps.hpp"
#include "image_gradient.hpp"
#include "load_cem.hpp"
#include "load_edg.hpp"
#include "merge_geom.hpp"
#include "prune_noise.hpp"
#include "write_cem.hpp"

#include <chrono>
#include <iostream>
#include <string>

enum class OutputFormat { Cem, Cemv, Both };

static void print_usage(const char* argv0) {
    std::cerr << "Usage: " << argv0 << " <input.edg> <input.cem> <input.image> [format] [output]\n"
              << "  Runs Step 3 (corner break) through final post-process.\n"
              << "  Corner ori_diff_th is fixed at pi/18 (matching main_TCG_CH.m).\n"
              << "  format: cem | cemv | both   (default: both)\n"
              << "  output: output path or stem\n"
              << "    - format=cem:  writes .cem  (default: "
                 "<image_stem>_tcg_cpp.cem)\n"
              << "    - format=cemv: writes .cemv (default: "
                 "<image_stem>_tcg_cpp.cemv)\n"
              << "    - format=both: writes .cem and sibling .cemv\n"
              << "                   (default stem: <image_stem>_tcg_cpp;\n"
              << "                    if output ends in .cem/.cemv, the other is "
                 "derived)\n"
              << "  Backward compatible: a 4th arg that is not cem|cemv|both is "
                 "treated as\n"
              << "  output path with format=both.\n";
}

static std::string image_stem(const std::string& img_path) {
    std::string stem = img_path;
    const auto slash = stem.find_last_of("/\\");
    if (slash != std::string::npos)
        stem = stem.substr(slash + 1);
    const auto dot = stem.find_last_of('.');
    if (dot != std::string::npos)
        stem = stem.substr(0, dot);
    return stem;
}

static bool ends_with_ci(const std::string& s, const std::string& ext) {
    if (s.size() < ext.size())
        return false;
    for (size_t i = 0; i < ext.size(); ++i) {
        char a = s[s.size() - ext.size() + i];
        char b = ext[i];
        if (a >= 'A' && a <= 'Z')
            a = static_cast<char>(a - 'A' + 'a');
        if (b >= 'A' && b <= 'Z')
            b = static_cast<char>(b - 'A' + 'a');
        if (a != b)
            return false;
    }
    return true;
}

static std::string replace_or_append_ext(const std::string& path, const std::string& new_ext) {
    // new_ext includes leading '.', e.g. ".cem"
    if (ends_with_ci(path, ".cemv"))
        return path.substr(0, path.size() - 5) + new_ext;
    if (ends_with_ci(path, ".cem"))
        return path.substr(0, path.size() - 4) + new_ext;
    return path + new_ext;
}

static bool parse_format(const std::string& s, OutputFormat& out) {
    if (s == "cem") {
        out = OutputFormat::Cem;
        return true;
    }
    if (s == "cemv") {
        out = OutputFormat::Cemv;
        return true;
    }
    if (s == "both") {
        out = OutputFormat::Both;
        return true;
    }
    return false;
}

static void resolve_output_paths(OutputFormat format, const std::string& img_path,
                                 const std::string& output_arg, bool have_output,
                                 std::string& out_cem, std::string& out_cemv, bool& do_cem,
                                 bool& do_cemv) {
    do_cem = (format == OutputFormat::Cem || format == OutputFormat::Both);
    do_cemv = (format == OutputFormat::Cemv || format == OutputFormat::Both);

    const std::string stem = image_stem(img_path) + "_tcg_cpp";
    if (!have_output) {
        out_cem = stem + ".cem";
        out_cemv = stem + ".cemv";
        return;
    }

    if (format == OutputFormat::Cem) {
        out_cem = ends_with_ci(output_arg, ".cem") ? output_arg
                                                   : replace_or_append_ext(output_arg, ".cem");
        out_cemv.clear();
    } else if (format == OutputFormat::Cemv) {
        out_cem.clear();
        out_cemv = ends_with_ci(output_arg, ".cemv") ? output_arg
                                                     : replace_or_append_ext(output_arg, ".cemv");
    } else {
        // both
        if (ends_with_ci(output_arg, ".cemv")) {
            out_cemv = output_arg;
            out_cem = replace_or_append_ext(output_arg, ".cem");
        } else if (ends_with_ci(output_arg, ".cem")) {
            out_cem = output_arg;
            out_cemv = replace_or_append_ext(output_arg, ".cemv");
        } else {
            out_cem = output_arg + ".cem";
            out_cemv = output_arg + ".cemv";
        }
    }
}

int main(int argc, char** argv) {
    if (argc < 4) {
        print_usage(argv[0]);
        return 1;
    }

    const std::string edg_path = argv[1];
    const std::string cem_path = argv[2];
    const std::string img_path = argv[3];

    OutputFormat format = OutputFormat::Both;
    std::string output_arg;
    bool have_output = false;
    if (argc >= 5) {
        if (parse_format(argv[4], format)) {
            if (argc >= 6) {
                output_arg = argv[5];
                have_output = true;
            }
        } else {
            // Backward compatible: 4th arg is output path, format=both.
            format = OutputFormat::Both;
            output_arg = argv[4];
            have_output = true;
        }
    }

    std::string out_cem, out_cemv;
    bool do_cem = false, do_cemv = false;
    resolve_output_paths(format, img_path, output_arg, have_output, out_cem, out_cemv, do_cem,
                         do_cemv);

    //> This corner orientation threshold is made constant, corresponding to
    // main_TCG_CH.m Step 3
    const double ori_th = M_PI / 18.0;

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

    const int h = grad.height;
    const int w = grad.width;

    std::cout << "[Summary] EDG " << edg.width << "x" << edg.height
              << " edges = " << edg.edges.size() << " | CEM contours = " << cem.contours.size()
              << " | IMG size = " << w << "x" << h << std::endl;

    //> Soft map from .edg (MATLAB edgemap_soft0)
    const std::vector<double>& edgemap_soft0 = edg.edgemap;

    tcg::BreakerParams bparams;
    bparams.nbr_num_edges = 20;
    bparams.corner_angle_th = M_PI / 6.0;

    //>=========== MATLAB contour_breaker_at_corner function =============
    auto t0 = std::chrono::steady_clock::now();
    tcg::CornerBreakResult broken =
        tcg::contour_breaker_at_corner(cem.contours, cem.contour_edge_idx, bparams, ori_th);
    auto t1 = std::chrono::steady_clock::now();
    const double contour_breaker_time_ms =
        std::chrono::duration<double, std::milli>(t1 - t0).count();
    std::cout << "[Step 3] number of curve fragments = " << broken.contours.size()
              << " / number of corners = " << broken.corner_pts.size()
              << " / time = " << contour_breaker_time_ms << " ms" << std::endl;
    //>============= MATLAB contour_breaker_at_corner function ==============

    //>============== MATLAB contour_fill_gaps_DP function ==========================
    tcg::GapFillParams gparams;
    gparams.DP_gap_range = 15;
    gparams.DP_angle_th = M_PI / 4.0;
    gparams.DP_contrast_th = 0.1;
    gparams.shape_gap_range = 8;
    gparams.shape_ori_range = M_PI / 9.0;

    t0 = std::chrono::steady_clock::now();
    tcg::GapFillResult filled = tcg::contour_fill_gaps_DP(broken.contours, broken.contour_edge_idx,
                                                          h, w, cem.edges, grad, gparams);
    t1 = std::chrono::steady_clock::now();
    const double contour_fill_gaps_time_ms =
        std::chrono::duration<double, std::milli>(t1 - t0).count();
    std::cout << "[Step 4] number of curve fragments = " << filled.contours.size()
              << " / number of edges = " << filled.edges.size()
              << " / time = " << contour_fill_gaps_time_ms << " ms" << std::endl;
    //>============= MATLAB contour_fill_gaps_DP function ==============

    //>============ MATLAB break_contours_at_T_junctions function =============
    t0 = std::chrono::steady_clock::now();
    tcg::TBreakResult tbroken =
        tcg::break_contours_at_T_junctions(filled.contours, filled.contour_edge_idx, filled.edges);
    t1 = std::chrono::steady_clock::now();
    const double break_T_time_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    std::cout << "[T-break] number of curve fragments = " << tbroken.contours.size()
              << " / time = " << break_T_time_ms << " ms" << std::endl;
    //>============ MATLAB break_contours_at_T_junctions function =============

    //>============MATLAB prune_noise_curves function ============
    tcg::PruneParams pparams;
    pparams.noise_len_th = 5.0;
    pparams.noise_prob_th = 0.05;

    t0 = std::chrono::steady_clock::now();
    tcg::PruneResult pruned = tcg::prune_noise_curves(tbroken.contours, tbroken.contour_edge_idx, h,
                                                      w, edgemap_soft0, pparams);
    t1 = std::chrono::steady_clock::now();
    const double prune_time_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    std::cout << "[Prune] number of curve fragments = " << pruned.contours.size()
              << " / time = " << prune_time_ms << " ms" << std::endl;
    //>=========== MATLAB prune_noise_curves function ============

    //>============= MATLAB merge_cfrags_graphical_model_geomfunction ============
    tcg::MergeGeomParams mparams;
    mparams.geom_merge_angle_th = M_PI / 6.0;
    mparams.nbr_num_edges = 20;

    t0 = std::chrono::steady_clock::now();
    tcg::MergeGeomResult merged = tcg::merge_cfrags_graphical_model_geom(
        pruned.contours, pruned.contour_edge_idx, tbroken.edges, mparams);
    t1 = std::chrono::steady_clock::now();
    const double merge_time_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    std::cout << "[Merge-geom] number of curve fragments = " << merged.contours.size()
              << " / time = " << merge_time_ms << " ms" << std::endl;
    //>============ MATLAB merge_cfrags_graphical_model_geom function ============

    //>============== MATLAB classify_junction_type_wrt_graph_BP function ==========
    tcg::ClassifyBPParams bpparams;
    bpparams.BP_merge_angle_th = M_PI / 9.0;
    bpparams.BP_nbr_num_edges = 20;
    bpparams.BP_clen_th = 15.0;

    t0 = std::chrono::steady_clock::now();
    tcg::ClassifyBPResult classified = tcg::classify_junction_type_wrt_graph_BP(
        merged.contours, merged.contour_edge_idx, tbroken.edges, bpparams);
    t1 = std::chrono::steady_clock::now();
    const double classify_bp_time_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    std::cout << "[Classify-BP] number of curve fragments = " << classified.contours.size()
              << " / number of T junctions = " << classified.T_junctions.size()
              << " / time = " << classify_bp_time_ms << " ms" << std::endl;
    //>============== MATLAB classify_junction_type_wrt_graph_BP
    // function ==============

    //>============== MATLAB contour_breaker_at_corner function (2nd pass) ==============
    t0 = std::chrono::steady_clock::now();
    tcg::CornerBreakResult broken2 =
        tcg::contour_breaker_at_corner(classified.contours, classified.contour_edge_idx, bparams);
    t1 = std::chrono::steady_clock::now();
    const double contour_breaker2_time_ms =
        std::chrono::duration<double, std::milli>(t1 - t0).count();
    std::cout << "[Corner2] number of curve fragments = " << broken2.contours.size()
              << " / number of corners = " << broken2.corner_pts.size()
              << " / time = " << contour_breaker2_time_ms << " ms" << std::endl;
    //>============== MATLAB contour_breaker_at_corner function (2nd pass) ==============

    //>============== MATLAB prune_noise_curves function (2nd pass) ==============
    t0 = std::chrono::steady_clock::now();
    tcg::PruneResult finalp = tcg::prune_noise_curves(broken2.contours, broken2.contour_edge_idx, h,
                                                      w, edgemap_soft0, pparams);
    t1 = std::chrono::steady_clock::now();
    const double prune2_time_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    std::cout << "[Final] number of curve fragments = " << finalp.contours.size()
              << " / time = " << prune2_time_ms << " ms" << std::endl;
    //>============== MATLAB prune_noise_curves function (2nd pass) ==============

    //> Write final contours (format: cem | cemv | both)
    if (do_cem) {
        if (!tcg::write_cem(out_cem, finalp.contours, h, w, err)) {
            std::cerr << "write_cem: " << err << "\n";
            return 1;
        }
        std::cout << "[Write] final contours -> " << out_cem << " (" << finalp.contours.size()
                  << " fragments)" << std::endl;
    }
    if (do_cemv) {
        // Port of util/io/det_save_cemv.m / Main_cem2cemv.m
        if (!tcg::write_cemv(out_cemv, finalp.contours, err)) {
            std::cerr << "write_cemv: " << err << "\n";
            return 1;
        }
        std::cout << "[Write] final contours (cemv) -> " << out_cemv << " ("
                  << finalp.contours.size() << " fragments)" << std::endl;
    }

    return 0;
}
