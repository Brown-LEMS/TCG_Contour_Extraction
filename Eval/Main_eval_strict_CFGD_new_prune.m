clear all; close all;
addpath (genpath('util/'));

GT_dir_path = 'Data/CFGD/';
img_path = 'Data/CFGD_img/';
cem_path_root = 'Data/TO_SEL_CFGD/final_curves/';
% cem_path_root = '/media/New_Volume/Research/Project_contour/SEL_cf_grouping/Tests/PB_TO_Kovesi/selection_test/';
% cem_path_root = '/media/New_Volume/Research/Project_contour/SEL_cf_grouping/Tests/PB_Kokkinos/selection_test/';
prefix = 'TO_SEL_CFGD_strict_eval_';

p_th_vec = [0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.96 0.99];

cem_folders = dir(cem_path_root);
cem_folders = cem_folders(end-length(p_th_vec)+1:end);
img_files = dir([img_path '*.jpg']);

for f = 1:length(cem_folders)
    cem_path = [cem_path_root cem_folders(f).name '/'];
    cem_folders(f).name
    
    sum_TP_GT = 0;
    sum_TP_CP = 0;
    sum_GT = 0;
    sum_CP = 0;
    for i = 1:length(img_files)
        disp([num2str(i) '/' num2str(length(img_files))]);
        img_name = img_files(i).name
        cp_file_path = [cem_path img_name(1:end-4) '.cem'];
        CP_CEM = load_contours(cp_file_path);
        GT_1_CEM = load_contours([GT_dir_path img_name(1:end-4) '_s1.cem']);
        GT_2_CEM = load_contours([GT_dir_path img_name(1:end-4) '_s2.cem']);
        GT_3_CEM = load_contours([GT_dir_path img_name(1:end-4) '_s3.cem']);
        
        [TP_GT_L, TP_CP_L, GT_L, CP_L] = Compare_Curve_Fragment_Maps_v2_new_score(GT_1_CEM, CP_CEM);
        sum_TP_GT = sum_TP_GT + TP_GT_L;
        sum_TP_CP = sum_TP_CP + TP_CP_L;
        sum_GT = sum_GT + GT_L;
        sum_CP = sum_CP + CP_L;
        
        [TP_GT_L, TP_CP_L, GT_L, CP_L] = Compare_Curve_Fragment_Maps_v2_new_score(GT_2_CEM, CP_CEM);
        sum_TP_GT = sum_TP_GT + TP_GT_L;
        sum_TP_CP = sum_TP_CP + TP_CP_L;
        sum_GT = sum_GT + GT_L;
        sum_CP = sum_CP + CP_L;

        [TP_GT_L, TP_CP_L, GT_L, CP_L] = Compare_Curve_Fragment_Maps_v2_new_score(GT_3_CEM, CP_CEM);
        sum_TP_GT = sum_TP_GT + TP_GT_L;
        sum_TP_CP = sum_TP_CP + TP_CP_L;
        sum_GT = sum_GT + GT_L;
        sum_CP = sum_CP + CP_L;

    end
    
    Recall(1,f) = sum_TP_GT/sum_GT;
    Precision(1,f) = sum_TP_CP/sum_CP;
    F(1,f) = 2*Recall(1,f)*Precision(1,f)/(Recall(1,f)+Precision(1,f));

end

% evaluate the unpruned
cem_path = cem_path_root;
sum_TP_GT = 0;
sum_TP_CP = 0;
sum_GT = 0;
sum_CP = 0;
for i = 1:length(img_files)
    disp([num2str(i) '/' num2str(length(img_files))]);
    img_name = img_files(i).name
    cp_file_path = [cem_path img_name(1:end-4) '.cem'];
    CP_CEM = load_contours(cp_file_path);
    GT_1_CEM = load_contours([GT_dir_path img_name(1:end-4) '_s1.cem']);
    GT_2_CEM = load_contours([GT_dir_path img_name(1:end-4) '_s2.cem']);
    GT_3_CEM = load_contours([GT_dir_path img_name(1:end-4) '_s3.cem']);

    [TP_GT_L, TP_CP_L, GT_L, CP_L] = Compare_Curve_Fragment_Maps_v2_new_score(GT_1_CEM, CP_CEM);
    sum_TP_GT = sum_TP_GT + TP_GT_L;
    sum_TP_CP = sum_TP_CP + TP_CP_L;
    sum_GT = sum_GT + GT_L;
    sum_CP = sum_CP + CP_L;

    [TP_GT_L, TP_CP_L, GT_L, CP_L] = Compare_Curve_Fragment_Maps_v2_new_score(GT_2_CEM, CP_CEM);
    sum_TP_GT = sum_TP_GT + TP_GT_L;
    sum_TP_CP = sum_TP_CP + TP_CP_L;
    sum_GT = sum_GT + GT_L;
    sum_CP = sum_CP + CP_L;

    [TP_GT_L, TP_CP_L, GT_L, CP_L] = Compare_Curve_Fragment_Maps_v2_new_score(GT_3_CEM, CP_CEM);
    sum_TP_GT = sum_TP_GT + TP_GT_L;
    sum_TP_CP = sum_TP_CP + TP_CP_L;
    sum_GT = sum_GT + GT_L;
    sum_CP = sum_CP + CP_L;

end

Recall = [sum_TP_GT/sum_GT, Recall];
Precision = [sum_TP_CP/sum_CP, Precision];
F = [2*Recall(1,f)*Precision(1,f)/(Recall(1,f)+Precision(1,f)), F];

Recall
Precision
F

PR_mat = [Recall;Precision;F];
save([prefix 'PR_mat.mat'], 'PR_mat');

figure;
plot(Recall, Precision, 'rs-');
xlabel('Recall');
ylabel('Precision');
axis square;
grid on;
axis([0 1 0 1]);

% print_pdf([prefix 'PR.pdf']);
