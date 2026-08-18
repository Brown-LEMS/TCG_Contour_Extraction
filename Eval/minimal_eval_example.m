addpath(genpath('util/'));

current_script_path = fileparts(mfilename('fullpath'));
gt_conout_path  = fullfile(current_script_path, "../example_data/n03425413_14351_to_tcg_matlab.cem");
cp_contour_path = fullfile(current_script_path, "../example_data/n03425413_14351_to_tcg_cpp.cem");

[TP_GT_L, TP_CP_L, GT_L, CP_L] = Compare_Curve_Fragment_Maps(gt_conout_path, cp_contour_path, 0, 1);

%> Calculate the recall, precision, and F-score
recall    = TP_GT_L / GT_L;
precision = TP_CP_L / CP_L;
F         = 2 * recall * precision / (recall + precision);
fprintf('Recall:    %f\n', recall);
fprintf('Precision: %f\n', precision);
fprintf('F-score:   %f\n', F);
