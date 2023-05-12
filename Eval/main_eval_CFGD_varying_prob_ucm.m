clear all; close all;
addpath (genpath('../util/'));

GT_dir_path = '../Data/CFGD/ground_truth/';
img_path = '../Data/CFGD/img/';
beta_src = '../training/';
cem_path = '../Data/CFGD/CFGD_SE/SE_ucm_CFGD/';
edge_path = '../Data/CFGD/CFGD_SE/edges/';
% cem_path_root = '/media/New_Volume/Research/Project_contour/SEL_cf_grouping/Tests/PB_TO_Kovesi/selection_test/';
% cem_path_root = '/media/New_Volume/Research/Project_contour/SEL_cf_grouping/Tests/PB_Kokkinos/selection_test/';
prefix = 'SE_ucm_CFGD_eval_vary_prob_';


% vary_num = [1 3 5 7 10 15 20 30 40 60 80 120 160 240 320 400 500 600 750];
vary_p = 0.95:-0.05:0.05;

% cem_folders = dir(cem_path_root);
% cem_folders = cem_folders(end-length(p_th_vec)+1:end);
img_files = dir([img_path '*.jpg']);

sum_TP_GT = zeros(length(img_files), length(vary_p));
sum_TP_CP = zeros(length(img_files), length(vary_p));
sum_GT = zeros(length(img_files), length(vary_p));
sum_CP = zeros(length(img_files), length(vary_p));
c_nums_all = zeros(length(img_files), length(vary_p));


for i = 1:length(img_files)
    disp([num2str(i) '/' num2str(length(img_files))]);
    img_name = img_files(i).name
    cp_file = [cem_path img_name(1:end-4) '_cfrag_maps.mat'];
    load(cp_file);
    GT_1_CEM = load_contours([GT_dir_path img_name(1:end-4) '_s1.cem']);
    GT_2_CEM = load_contours([GT_dir_path img_name(1:end-4) '_s2.cem']);
    GT_3_CEM = load_contours([GT_dir_path img_name(1:end-4) '_s3.cem']);

    
    % compute probs of all tha contours
    img = imread([img_path img_files(i).name]);


    disp('evaluate curve fragments');
    for f = 1:length(vary_p)
        
%         num_th = min(vary_num(f), length(new_cfrags))        
        p_th = vary_p(f)
        
        CP_CEM =  cell(3,1);
        CP_CEM{1} = cfrags_maps.edges{f};
        CP_CEM{2} = cfrags_maps.cfrags{f};
        c_nums_all(i,f) = length(cfrags_maps.cfrags{f});
        
        [TP_GT_L, TP_CP_L, GT_L, CP_L] = Compare_Curve_Fragment_Maps_2(GT_1_CEM, CP_CEM);
        sum_TP_GT(i,f) = sum_TP_GT(i,f) + TP_GT_L;
        sum_TP_CP(i,f) = sum_TP_CP(i,f) + TP_CP_L;
        sum_GT(i,f) = sum_GT(i,f) + GT_L;
        sum_CP(i,f) = sum_CP(i,f) + CP_L;

        [TP_GT_L, TP_CP_L, GT_L, CP_L] = Compare_Curve_Fragment_Maps_2(GT_2_CEM, CP_CEM);
        sum_TP_GT(i,f) = sum_TP_GT(i,f) + TP_GT_L;
        sum_TP_CP(i,f) = sum_TP_CP(i,f) + TP_CP_L;
        sum_GT(i,f) = sum_GT(i,f) + GT_L;
        sum_CP(i,f) = sum_CP(i,f) + CP_L;

        [TP_GT_L, TP_CP_L, GT_L, CP_L] = Compare_Curve_Fragment_Maps_2(GT_3_CEM, CP_CEM);
        sum_TP_GT(i,f) = sum_TP_GT(i,f) + TP_GT_L;
        sum_TP_CP(i,f) = sum_TP_CP(i,f) + TP_CP_L;
        sum_GT(i,f) = sum_GT(i,f) + GT_L;
        sum_CP(i,f) = sum_CP(i,f) + CP_L;

    end

end

Recall = sum(sum_TP_GT)./sum(sum_GT);
Precision = sum(sum_TP_CP)./sum(sum_CP);
F = 2*Recall.*Precision./(Recall+Precision); 
vary_num = mean(c_nums_all,1);

save([prefix 'PR_v2.mat'], 'Recall', 'Precision', 'F', 'vary_p', 'vary_num','sum_TP_GT', 'sum_TP_CP', 'sum_GT', 'sum_CP', 'c_nums_all');

figure;
plot(Recall, Precision, 'rs-');
xlabel('Recall');
ylabel('Precision');
axis square;
grid on;
axis([0 1 0 1]);

% print_pdf([prefix 'PR.pdf']);
