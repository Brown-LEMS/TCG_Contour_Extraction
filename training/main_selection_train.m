clear all; close all;
addpath (genpath('../util/'));



%% select pos/neg contour samples
img_src_path = '../Data/CFGD_img/';
edg_src_path = '../Data/TO_SEL_CFGD/edges/';
cem_src_path = '../Data/TO_SEL_CFGD/final_curves/';
% vis_dst_path = 'Tests/TO_SEL_cfrags_CFGD2_vis/';
% mkdir(vis_dst_path);
% hist_dst_path = 'hists/TO_SEL_hists_cues_for_merging/';
% mkdir(hist_dst_path);
gt_src_path = '../Data/CFGD/';
prune_len_thresh = 0;
prefix = 'TO_SEL_';

cem_files = dir([cem_src_path '*.cem']);
pos_contours = cell(1,0);
neg_contours = cell(1,0);
pos_features = [];
neg_features = [];

for c = 1:length(cem_files)
    disp([num2str(c) '/'  num2str(length(cem_files))]);
    
    input_file = [cem_src_path cem_files(c).name];
    [CEM, edges, cf_idx] = load_contours(input_file);
    [~, edgemap, thetamap] = load_edg([edg_src_path cem_files(c).name(1:end-3) 'edg']);
    
    img = imread([img_src_path cem_files(c).name(1:end-4) '.jpg']);    
    % compute hsv space map
    hsv_img = rgb2hsv(img);
    
    gt_file_1 = [gt_src_path cem_files(c).name(1:end-4) '_s1.cem'];
    CEM_gt_1 = load_contours(gt_file_1);
    
    gt_file_2 = [gt_src_path cem_files(c).name(1:end-4) '_s2.cem'];
    CEM_gt_2 = load_contours(gt_file_2);
    
    gt_file_3 = [gt_src_path cem_files(c).name(1:end-4) '_s3.cem'];
    CEM_gt_3 = load_contours(gt_file_3);
    
%     % use the GT curve fragments confirmed by all users as pos samples
%     CEM_pos_intersect = CEM_gt_1;
%     CEM_pos_intersect{3,1} = [];
% 
%     [TP_GT_L, TP_CP_L, GT_L, CP_L, match_GT_frags, match_CP_frags] = Compare_Curve_Fragment_Maps_2(CEM_gt_1, CEM_gt_2, 3);
% 
%     CEM_pos_intersect{2,1} = match_GT_frags;
%     
%     [TP_GT_L, TP_CP_L, GT_L, CP_L, match_GT_frags, match_CP_frags] = Compare_Curve_Fragment_Maps_2(CEM_pos_intersect, CEM_gt_3, 3);
%     CEM_pos_intersect{2,1} = match_GT_frags;
    
%  it seems that using CP contours confirmed by all GT is more appropiate
%  for training

    CEM_pos_intersect = CEM;
    CEM_pos_intersect{3,1} = [];
    [TP_GT_L, TP_CP_L, GT_L, CP_L, match_GT_frags, match_CP_frags] = Compare_Curve_Fragment_Maps_2(CEM_gt_1, CEM, prune_len_thresh);

    CEM_pos_intersect{2,1} = match_CP_frags;
    
    [TP_GT_L, TP_CP_L, GT_L, CP_L, match_GT_frags, match_CP_frags] = Compare_Curve_Fragment_Maps_2(CEM_gt_2, CEM_pos_intersect, prune_len_thresh);
    CEM_pos_intersect{2,1} = match_CP_frags;    
    
    [TP_GT_L, TP_CP_L, GT_L, CP_L, match_GT_frags, match_CP_frags] = Compare_Curve_Fragment_Maps_2(CEM_gt_3, CEM_pos_intersect, prune_len_thresh);
    CEM_pos_intersect{2,1} = match_CP_frags;
    
    pos_contours = [pos_contours CEM_pos_intersect{2,1}];
    
    for i=1:length(CEM_pos_intersect{2,1})
        [bg_grad, sat_grad, hue_grad, abs_k, edge_sparsity, wigg, len, mean_conf] = curve_fragment_cues(CEM_pos_intersect{2,1}{i}, hsv_img, edgemap);
        pos_features = [pos_features [bg_grad; sat_grad; hue_grad; abs_k; edge_sparsity; wigg; len; mean_conf]];
    end
    
    % use the computed curve fragments which are not matched to any GT as
    % neg samples
    CEM_FP = CEM;
    CEM_FP{3,1} = [];
    [TP_GT_L, TP_CP_L, GT_L, CP_L, match_GT_frags, match_CP_frags, miss_GT_frags, extra_CP_frags] = Compare_Curve_Fragment_Maps_2(CEM_gt_1, CEM_FP, prune_len_thresh);
    CEM_FP{2,1} = extra_CP_frags;
    [TP_GT_L, TP_CP_L, GT_L, CP_L, match_GT_frags, match_CP_frags, miss_GT_frags, extra_CP_frags] = Compare_Curve_Fragment_Maps_2(CEM_gt_2, CEM_FP, prune_len_thresh);
    CEM_FP{2,1} = extra_CP_frags;
    [TP_GT_L, TP_CP_L, GT_L, CP_L, match_GT_frags, match_CP_frags, miss_GT_frags, extra_CP_frags] = Compare_Curve_Fragment_Maps_2(CEM_gt_3, CEM_FP, prune_len_thresh);
    CEM_FP{2,1} = extra_CP_frags;
    
    neg_contours = [neg_contours CEM_FP{2,1}];
    for i=1:length(CEM_FP{2,1})
        [bg_grad, sat_grad, hue_grad, abs_k, edge_sparsity, wigg, len, mean_conf] = curve_fragment_cues(CEM_FP{2,1}{i}, hsv_img, edgemap);
        neg_features = [neg_features [bg_grad; sat_grad; hue_grad; abs_k; edge_sparsity; wigg; len; mean_conf]];
    end
%     keyboard;
    
end

% save('PB_TO_SEL_pos_contours.mat', 'pos_contours');
% save('PB_TO_SEL_neg_contours.mat', 'neg_contours');
save([prefix 'pos_features.mat'], 'pos_features');
save([prefix 'neg_features.mat'], 'neg_features');
%% learn logistic regression beta
% load('PB_TO_SEL_pos_contours.mat');
% load('PB_TO_SEL_neg_contours.mat');
load([prefix 'pos_features.mat']);
load([prefix 'neg_features.mat']);

%%%%%%%%%%%%%%  rule out outliers, and construct classifier focusing on the
%%%%%%%%%%%%%%  confusing range
outlier_ids_pos = [];
outlier_ids_neg = [];

for d = 1:size(pos_features, 1)
   f_min = max(min(pos_features(d,:)), min(neg_features(d,:)));
   f_max = min(max(pos_features(d,:)), max(neg_features(d,:)));

   outlier_ids_pos = [outlier_ids_pos find(pos_features(d,:)<f_min)];
   outlier_ids_pos = [outlier_ids_pos find(pos_features(d,:)>f_max)];
   
   outlier_ids_neg = [outlier_ids_neg find(neg_features(d,:)<f_min)];
   outlier_ids_neg = [outlier_ids_neg find(neg_features(d,:)>f_max)];
   
end

outlier_ids_pos = unique(outlier_ids_pos);
outlier_ids_neg = unique(outlier_ids_neg);

pos_features(:, outlier_ids_pos) = [];
neg_features(:, outlier_ids_neg) = [];



%%%%%%%%%%%%%  learn the parameters of logistic regression
cnt_pos = size(pos_features,2);
cnt_neg = size(neg_features,2);

if(cnt_pos < cnt_neg)
    idx_neg = randperm(cnt_neg);

    % neg_contours_2 = neg_contours(1, idx_neg(1:cnt_pos));
    neg_features = neg_features(:, idx_neg(1:cnt_pos));

    Y = [ones(1, cnt_pos), zeros(1, cnt_pos)];
    % contour_samples = [pos_contours neg_contours_2];
else
    idx_pos = randperm(cnt_pos);
    
    pos_features = pos_features(:, idx_pos(1:cnt_neg));
    Y = [ones(1, cnt_neg), zeros(1, cnt_neg)];
    
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

P_all = 1 ./ (1 + exp(-(Features_1)*beta'));

varying_th = (1:0.5:9)/10;

error = [];
for i = 1:length(varying_th)
    estimates = zeros(size(P_all));
    estimates(find(P_all> varying_th(i))) = 1;
    
    error = [error sum(abs(estimates - Y))/length(P_all)];
    
end

plot(varying_th, error, 'r');

min(error)


%% visuliaze hists
feature_names = cell(size(Features, 1), 1);
feature_names{1} = 'BgGrad';
feature_names{2} = 'SatGrad';
feature_names{3} = 'HueGrad';
feature_names{4} = 'AbsCurvature';
feature_names{5} = 'EdgeSparsity';
feature_names{6} = 'Wigg';
feature_names{7} = 'Len';
feature_names{8} = 'AvgConf';


pos_features = Features(:, find(Y>0));
neg_features = Features(:, find(Y==0));

for d = 1:size(Features, 1)
    [pos_hist, pos_X] = hist(pos_features(d,:));
    [neg_hist, neg_X] = hist(neg_features(d,:));

    H = figure; hold on;
    plot(pos_X, pos_hist, '-r', 'LineWidth',8);
    plot(neg_X, neg_hist, '-b', 'LineWidth',8);
    hold off;
    AX = legend([ ' pos'], [ ' neg'], 'Location', 'NorthEast');
    LEG = findobj(AX,'type','text');
    set(LEG,'FontSize',25)
    TX = title(feature_names{d});
    TT = findobj(TX,'type','text');
    set(TT,'FontSize',25)
    print(H,'-dpng','-r0',['hist_selection/hist_' feature_names{d} '.png']);

end

pos_P = P_all(find(Y>0));
neg_P = P_all(find(Y==0));
[pos_hist, pos_X] = hist(pos_P);
[neg_hist, neg_X] = hist(neg_P);

H = figure; hold on;
plot(pos_X, pos_hist, '-r', 'LineWidth',8);
plot(neg_X, neg_hist, '-b', 'LineWidth',8);
hold off;
AX = legend('pos', 'neg', 'Location', 'NorthEast');
LEG = findobj(AX,'type','text');
set(LEG,'FontSize',25)
TX = title('Output Prob');
TT = findobj(TX,'type','text');
set(TT,'FontSize',25)
print(H,'-dpng','-r0',['hist_selection/hist_prob' '.png']);