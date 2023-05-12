clear all; close all;
addpath (genpath('../util/'));

maxDist = 0.01;


GT_dir_path = '../Data/VOC2007_GT/';
img_path = '../Data/VOC2007_img/';
beta_src = '../training_voc/';
cem_path = '../Data/TO_Kovesi_VOC2007/';
edge_path = '../Data/TO_SEL_VOC2007/edges/';
% cem_path_root = '/media/New_Volume/Research/Project_contour/SEL_cf_grouping/Tests/PB_TO_Kovesi/selection_test/';
% cem_path_root = '/media/New_Volume/Research/Project_contour/SEL_cf_grouping/Tests/PB_Kokkinos/selection_test/';
prefix = 'TO_Kovesi_VOC2007_eval_vary_num_';



vary_num = [1 3 5 7 10 15 20 30 40 60 80 120 160 240 320];

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
    cp_file_path = [cem_path img_name(1:end-4) '.cem'];
    CP_CEM = load_contours(cp_file_path);

    load([GT_dir_path img_name(1:end-4) '.mat' ]);

    % compute probs of all tha contours
    % compute hsv space map
    img = imread([img_path img_files(i).name]);
    hsv_img = rgb2hsv(img);
    [~, edgemap, thetamap] = load_edg([edge_path img_files(i).name(1:end-4) '.edg']);

    new_cfrags = CP_CEM{2};

    disp('rank curve fragments');
    new_cfrags = CP_CEM{2};
    P_vec = [];
    for c = 1:length(new_cfrags)
        p = contour_length_mex(new_cfrags{c}');
        P_vec = [P_vec p];
    end 
    
    % rank the cfrags
    [~, sort_id] = sort(P_vec, 2, 'descend');
    new_cfrags = new_cfrags(sort_id);

    disp('evaluate curve fragments');
    for f = 1:length(vary_num)
        
        num_th = min(vary_num(f), length(new_cfrags))        
        CP_CEM{2} = new_cfrags(1:num_th);
        
        % match curve frags to each object's silhouette and compute PR
        % pixels
        for igt = 1:length(gt_edgemap)
            disp(['object: ' num2str(igt)]);

            [sum_tp_gt, sum_tp_cp, sum_gt, sum_cp] = cfrags_vs_gt_edgemap(CP_CEM{2}, edgemap, gt_edgemap{igt}, maxDist, 1);
            sum_TP_GT(i,f) = sum_TP_GT(i,f) + sum_tp_gt;
            sum_TP_CP(i,f) = sum_TP_CP(i,f) + sum_tp_cp;
            sum_GT(i,f) = sum_GT(i,f) + sum_gt;
            sum_CP(i,f) = sum_CP(i,f) + sum_cp;
        end

    end
    % all the curve frags
    CP_CEM{2} = new_cfrags;
    for igt = 1:length(gt_edgemap)
        disp(['object: ' num2str(igt)]);
        [sum_tp_gt, sum_tp_cp, sum_gt, sum_cp] = cfrags_vs_gt_edgemap(CP_CEM{2}, edgemap, gt_edgemap{igt}, maxDist,1);

        sum_TP_GT(i,f+1) = sum_TP_GT(i,f+1) + sum_tp_gt;
        sum_TP_CP(i,f+1) = sum_TP_CP(i,f+1) + sum_tp_cp;
        sum_GT(i,f+1) = sum_GT(i,f+1) + sum_gt;
        sum_CP(i,f+1) = sum_CP(i,f+1) + sum_cp;
    end
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
