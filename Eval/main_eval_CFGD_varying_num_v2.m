clear all; close all;
addpath (genpath('../util/'));

GT_dir_path = '../Data/CFGD/ground_truth/';
img_path = '../Data/CFGD/img/';
beta_src = '../training/';
cem_path = '../Data/CFGD/CFGD_SE/SE_TCG_CFGD/';
edge_path = '../Data/CFGD/CFGD_SE/edges/';
% cem_path_root = '/media/New_Volume/Research/Project_contour/SEL_cf_grouping/Tests/PB_TO_Kovesi/selection_test/';
% cem_path_root = '/media/New_Volume/Research/Project_contour/SEL_cf_grouping/Tests/PB_Kokkinos/selection_test/';
prefix = 'SE_TCG_CFGD_eval_vary_num_v2_';
maxDist = 0.0075;

% load in beta
beta_prefix = 'SE_SEL_';

input_beta_2 = load([beta_src beta_prefix 'beta_of_cues_for_seletion.txt']);
fmean_2 = input_beta_2(1,:);
fstd_2 = input_beta_2(2,:);
beta_2 = input_beta_2(3,:);
beta_2 = beta_2./fstd_2;


vary_num = [1 3 5 7 10 15 20 30 40 60 80 120 160 240 320 400 500 600 750];

% cem_folders = dir(cem_path_root);
% cem_folders = cem_folders(end-length(p_th_vec)+1:end);
img_files = dir([img_path '*.jpg']);

sum_TP_GT = zeros(length(img_files), 1+length(vary_num));
sum_TP_CP = zeros(length(img_files), 1+length(vary_num));
sum_GT = zeros(length(img_files), 1+length(vary_num));
sum_CP = zeros(length(img_files), 1+length(vary_num));

for i = 1:length(img_files)
    disp([num2str(i) '/' num2str(length(img_files))]);
    img_name = img_files(i).name
    img = imread([img_path img_files(i).name]);
    [h,w,~] = size(img);
    cp_file_path = [cem_path img_name(1:end-4) '.cem'];
    CP_CEM = load_contours(cp_file_path);
    GT_1_CEM = load_contours([GT_dir_path img_name(1:end-4) '_s1.cem']);
    GT_2_CEM = load_contours([GT_dir_path img_name(1:end-4) '_s2.cem']);
    GT_3_CEM = load_contours([GT_dir_path img_name(1:end-4) '_s3.cem']);

    EdgeGroupMap_GT_1 = convert_cfrags_to_EdgeGroupMap (GT_1_CEM{2}, h, w);
    EdgeGroupMap_GT_2 = convert_cfrags_to_EdgeGroupMap (GT_2_CEM{2}, h, w);
    EdgeGroupMap_GT_3 = convert_cfrags_to_EdgeGroupMap (GT_3_CEM{2}, h, w);
    
%     % compute probs of all tha contours
%     % compute hsv space map
%     hsv_img = rgb2hsv(img);
%     [~, edgemap, thetamap] = load_edg([edge_path img_files(i).name(1:end-4) '.edg']);

    disp('rank curve fragments');
    new_cfrags = CP_CEM{2};
    P_vec = [];
    for c = 1:length(new_cfrags)
%         [bg_grad, sat_grad, hue_grad, abs_k, edge_sparsity, wigg, len, mean_conf] = curve_fragment_cues(new_cfrags{c}, hsv_img, edgemap);
% 
%         p = 1 / (1 + exp(-([1, bg_grad, sat_grad, hue_grad, abs_k, edge_sparsity, wigg, len, mean_conf]-fmean_2)*beta_2'));
        p = sum(new_cfrags{c}(:,4));
        P_vec = [P_vec p];
    end 
    
    % rank the cfrags
    [~, sort_id] = sort(P_vec, 2, 'descend');
    new_cfrags = new_cfrags(sort_id);
    EdgeGroupMap_CP = convert_cfrags_to_EdgeGroupMap (new_cfrags, h, w);

    
    disp('evaluate curve fragments');
    for f = 1:length(vary_num)
        
        num_th = min(vary_num(f), length(new_cfrags))       
        cur_EdgeGroupMap = EdgeGroupMap_CP;
        cur_EdgeGroupMap (cur_EdgeGroupMap>num_th) = 0;
        [cntP, sumP, cntR, sumR] = matching_edge_grouping_maps_One2Multi(EdgeGroupMap_GT_1, cur_EdgeGroupMap, maxDist);
        sum_TP_GT(i,f) = sum_TP_GT(i,f) + cntR;
        sum_TP_CP(i,f) = sum_TP_CP(i,f) + cntP;
        sum_GT(i,f) = sum_GT(i,f) + sumR;
        sum_CP(i,f) = sum_CP(i,f) + sumP;        

        [cntP, sumP, cntR, sumR] = matching_edge_grouping_maps_One2Multi(EdgeGroupMap_GT_2, cur_EdgeGroupMap, maxDist);
        sum_TP_GT(i,f) = sum_TP_GT(i,f) + cntR;
        sum_TP_CP(i,f) = sum_TP_CP(i,f) + cntP;
        sum_GT(i,f) = sum_GT(i,f) + sumR;
        sum_CP(i,f) = sum_CP(i,f) + sumP; 

        [cntP, sumP, cntR, sumR] = matching_edge_grouping_maps_One2Multi(EdgeGroupMap_GT_3, cur_EdgeGroupMap, maxDist);
        sum_TP_GT(i,f) = sum_TP_GT(i,f) + cntR;
        sum_TP_CP(i,f) = sum_TP_CP(i,f) + cntP;
        sum_GT(i,f) = sum_GT(i,f) + sumR;
        sum_CP(i,f) = sum_CP(i,f) + sumP; 

    end
    [cntP, sumP, cntR, sumR] = matching_edge_grouping_maps_One2Multi(EdgeGroupMap_GT_1, EdgeGroupMap_CP, maxDist);
    sum_TP_GT(i,f+1) = sum_TP_GT(i,f) + cntR;
    sum_TP_CP(i,f+1) = sum_TP_CP(i,f) + cntP;
    sum_GT(i,f+1) = sum_GT(i,f) + sumR;
    sum_CP(i,f+1) = sum_CP(i,f) + sumP;  

    [cntP, sumP, cntR, sumR] = matching_edge_grouping_maps_One2Multi(EdgeGroupMap_GT_2, EdgeGroupMap_CP, maxDist);
    sum_TP_GT(i,f+1) = sum_TP_GT(i,f) + cntR;
    sum_TP_CP(i,f+1) = sum_TP_CP(i,f) + cntP;
    sum_GT(i,f+1) = sum_GT(i,f) + sumR;
    sum_CP(i,f+1) = sum_CP(i,f) + sumP;     
    
    [cntP, sumP, cntR, sumR] = matching_edge_grouping_maps_One2Multi(EdgeGroupMap_GT_3, EdgeGroupMap_CP, maxDist);
    sum_TP_GT(i,f+1) = sum_TP_GT(i,f) + cntR;
    sum_TP_CP(i,f+1) = sum_TP_CP(i,f) + cntP;
    sum_GT(i,f+1) = sum_GT(i,f) + sumR;
    sum_CP(i,f+1) = sum_CP(i,f) + sumP;  
    
    
    Recall = sum(sum_TP_GT)./sum(sum_GT)
    Precision = sum(sum_TP_CP)./sum(sum_CP)
end

Recall = sum(sum_TP_GT)./sum(sum_GT);
Precision = sum(sum_TP_CP)./sum(sum_CP);
F = 2*Recall.*Precision./(Recall+Precision); 

save([prefix 'PR.mat'], 'Recall', 'Precision', 'F', 'vary_num', 'sum_TP_GT', 'sum_TP_CP', 'sum_GT', 'sum_CP');

figure;
plot(Recall, Precision, 'rs-');
xlabel('Recall');
ylabel('Precision');
axis square;
grid on;
axis([0 1 0 1]);

% print_pdf([prefix 'PR.pdf']);
