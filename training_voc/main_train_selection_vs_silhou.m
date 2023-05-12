clear all; close all;
addpath(genpath('../util'));

patio_th = 0.7; % contours with GT edge over this threshold will be positive samples
len_th = 5;
maxDist = 0.01;
prefix = 'gPb_SEL_';

img_src_path = '../Data/VOC2007_img/';
edg_src_path = '../Data/gPb_SEL_VOC2007/edges/';
cem_src_path = '../Data/gPb_SEL_VOC2007/final_curves/';
gt_src_path = '../Data/VOC2007_GT/';

%% compute features and corresponding labels/probs

cem_files = dir([cem_src_path '*.cem']);
% cem_files = dir([gt_src_path '*.mat']);

pos_features = [];
neg_features = [];
pos_probs = [];
neg_probs = [];

for c = 1:length(cem_files)  
    disp([num2str(c) '/' num2str(length(cem_files))]);
    name = cem_files(c).name(1:end-4)
    input_file = [cem_src_path name '.cem'];
    [CEM, edges, cf_idx] = load_contours(input_file);
    [~, edgemap, thetamap] = load_edg([edg_src_path name '.edg']);
    
    img = imread([img_src_path name '.jpg']);    
    % compute hsv space map
    hsv_img = rgb2hsv(img);
    
    load([gt_src_path name '.mat' ]);
    
    % match silhouette of each object individually to select positive
    % samples
    disp('Selecting Pos Samples');

    gt_edgemap_union = zeros(size(gt_edgemap{1}));
    for igt = 1:length(gt_edgemap)
        disp(['object: ' num2str(igt)]);
        gt_edgemap_union = gt_edgemap_union | gt_edgemap{igt};
        
        [match1,match2] = correspondPixels(255*gt_edgemap{igt}, 255*(edgemap>0), maxDist);

    %         thetamap = thetamap.*double(match2>0) | (edgemap>0);

        final_gt_edgemap = uint8((match2>0) | (gt_edgemap{igt}>0))*255;

        [fg_contours, fg_edgemap, bg_contours, pos_probs_i] = refine_fg_cfrags_using_gt_edgemap(CEM, final_gt_edgemap, patio_th, len_th);

%         imshow(fg_edgemap)

        pos_probs = [pos_probs pos_probs_i];
        for i=1:length(fg_contours)
            [bg_grad, sat_grad, hue_grad, abs_k, edge_sparsity, wigg, len, mean_conf, area] = curve_fragment_cues(fg_contours{i}, hsv_img, edgemap);
            pos_features = [pos_features [bg_grad; sat_grad; hue_grad; abs_k; edge_sparsity; wigg; len; mean_conf; area]];
        end

        
    
    end
    
    % match to the union edge map to select negative samples
    disp('Selecting Neg Samples');
    [match1,match2] = correspondPixels(255*gt_edgemap_union, 255*(edgemap>0), maxDist);

%         thetamap = thetamap.*double(match2>0) | (edgemap>0);

    final_gt_edgemap = uint8((match2>0) | (gt_edgemap_union>0))*255;
    
    [fg_contours, fg_edgemap, bg_contours, pos_probs_i, neg_probs_i] = refine_fg_cfrags_using_gt_edgemap(CEM, final_gt_edgemap, patio_th, len_th);

    neg_probs = [neg_probs neg_probs_i];
    for i=1:length(bg_contours)
        [bg_grad, sat_grad, hue_grad, abs_k, edge_sparsity, wigg, len, mean_conf, area] = curve_fragment_cues(bg_contours{i}, hsv_img, edgemap);
        neg_features = [neg_features [bg_grad; sat_grad; hue_grad; abs_k; edge_sparsity; wigg; len; mean_conf; area]];
    end
%     keyboard;
end

save([prefix 'pos_features.mat'], 'pos_features', 'pos_probs');
save([prefix 'neg_features.mat'], 'neg_features', 'neg_probs');

% %% learn logistic regression beta (version 1, using pos, neg examples)
% 
% load([prefix 'pos_features.mat']);
% load([prefix 'neg_features.mat']);
% 
% disp('Train contour pruning model for silhouette');
% cnt_pos = size(pos_features,2);
% cnt_neg = size(neg_features,2);
% 
% if(cnt_pos < cnt_neg)
%     idx_neg = randperm(cnt_neg);
% 
%     % neg_contours_2 = neg_contours(1, idx_neg(1:cnt_pos));
%     neg_features = neg_features(:, idx_neg(1:cnt_pos));
% 
%     Y = [ones(1, cnt_pos), zeros(1, cnt_pos)];
%     % contour_samples = [pos_contours neg_contours_2];
% else
%     idx_pos = randperm(cnt_pos);
%     
%     pos_features = pos_features(:, idx_pos(1:cnt_neg));
%     Y = [ones(1, cnt_neg), zeros(1, cnt_neg)];
%     
% end
% 
% Features = [pos_features, neg_features];
% Features_1 = [ones(size(Features,2),1) Features'];  Y = Y';
% fstd = std(Features_1);
% fstd = fstd + (fstd==0);
% fmean = mean(Features_1);
% fmean(1) = 0;
% 
% Features_1 = (Features_1- repmat(fmean,size(Features_1,1),1)) ./ repmat(fstd,size(Features_1,1),1);
% 
% % fit the model
% fprintf(2,'Fitting model...\n');
% beta = logist2(Y,Features_1);
% beta = beta';
% 
% save_name = [prefix 'beta_of_cues_for_seletion.txt'];
% save(save_name,'fmean','fstd','beta','-ascii');

%% learn logistic regression beta (version 2, regressor)

load([prefix 'pos_features.mat']);
load([prefix 'neg_features.mat']);

disp('Train contour pruning model for silhouette');
cnt_pos = size(pos_features,2);
cnt_neg = size(neg_features,2);

if(cnt_pos < cnt_neg)
    idx_neg = randperm(cnt_neg);

    % neg_contours_2 = neg_contours(1, idx_neg(1:cnt_pos));
    neg_features = neg_features(:, idx_neg(1:cnt_pos));

    Y = [pos_probs(1: cnt_pos), neg_probs(1: cnt_pos)];
    % contour_samples = [pos_contours neg_contours_2];
else
    idx_pos = randperm(cnt_pos);
    
    pos_features = pos_features(:, idx_pos(1:cnt_neg));
    Y = [pos_probs(1: cnt_neg), neg_probs(1: cnt_neg)];
    
end

Features = [pos_features, neg_features];
Features_1 = [ones(size(Features,2),1) Features'];  Y = Y';
fstd = std(Features_1);
fstd = fstd + (fstd==0);
fmean = mean(Features_1);
fmean(1) = 0;

Features_1 = (Features_1- repmat(fmean,size(Features_1,1),1)) ./ repmat(fstd,size(Features_1,1),1);

% fit the model
fprintf(2,'Fitting model...\n');
beta = logist2(Y,Features_1);
beta = beta';

save_name = [prefix 'beta_of_cues_for_seletion.txt'];
save(save_name,'fmean','fstd','beta','-ascii');