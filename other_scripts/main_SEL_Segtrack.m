clear all; close all;
addpath (genpath('util/'));

%% set paths
edge_src_path = '/media/New_Volume/Research/Project_tracking/SegTrack_skel_tracking/SegTrackv2/SETO/girl/';
cfrags_dst_path = '/media/New_Volume/Research/Project_tracking/SegTrack_skel_tracking/SegTrackv2/SETO/girl/cfrags/';
mkdir(cfrags_dst_path);
img_src_path = '/media/New_Volume/Research/Project_tracking/SegTrack_skel_tracking/SegTrackv2/JPEGImages/girl/';
cem_src_path = cfrags_dst_path;
% broken_cfrag_path = 'Data/TO_SEL_VOC2007/broken_cfrags/';
% mkdir(broken_cfrag_path);
final_cem_path = '/media/New_Volume/Research/Project_tracking/SegTrack_skel_tracking/SegTrackv2/SETO/girl/final_curves/';
mkdir(final_cem_path);
% final_cem_path = './'
beta_src = 'training/';
top_K = 100;
img_files = dir([img_src_path '*.bmp']);


%% load in beta
prefix = 'SE_SEL_';

input_beta_0 = load([beta_src prefix 'beta_of_cues_for_merging.txt']);
fmean_0 = input_beta_0(1,:);
fstd_0 = input_beta_0(2,:);
beta_0 = input_beta_0(3,:);
beta_0 = beta_0./fstd_0;

input_beta_1 = load([beta_src prefix 'beta_of_geomcon_cue_for_merging.txt']);
fmean_1 = input_beta_1(1,:);
fstd_1 = input_beta_1(2,:);
beta_1 = input_beta_1(3,:);
beta_1 = beta_1./fstd_1;

input_beta_2 = load([beta_src prefix 'beta_of_cues_for_seletion.txt']);
fmean_2 = input_beta_2(1,:);
fstd_2 = input_beta_2(2,:);
beta_2 = input_beta_2(3,:);
beta_2 = beta_2./fstd_2;

params.beta_0 = beta_0;
params.fmean_0 = fmean_0;
params.beta_1 = beta_1;
params.fmean_1 = fmean_1;
params.beta_2 = beta_2;
params.fmean_2 = fmean_2;
%% set parameters

nbr_num_edges = 20; % only consider number of edges close to the connecting points for feature computation % also for matching in building ground truth
params.diag_of_train = 578.275;
params.nbr_len_th = 5; % short curve under this length will be grouped due to geometry.
params.merge_th = 0.4;
params.merge_th_geom = 0.4;
params.max_iter = 2;    


% params for computing texton
no = 6;
ss = 1;
ns = 2;
sc = sqrt(2);
el = 2;
k = 64;
fname = sprintf( ...
    'unitex_%.2g_%.2g_%.2g_%.2g_%.2g_%d.mat',no,ss,ns,sc,el,k);
textonData = load(fname); % defines fb,tex,tsim
    
    
%% Extrac Curve Fragments from edges

for c = 1:length(img_files)
    
    disp([num2str(c) '/'  num2str(length(img_files))]); 
    
    %% Compute cfrags from edges: curvelets, unambiguous cfrags, HPG cfrags
    input = [edge_src_path img_files(c).name(1:end-4) '.edg'];
    output = [cfrags_dst_path img_files(c).name(1:end-4) '.cem'];
    
    cmd = ['!util/dborl_compute_curve_frags ' input ' ' output];
    eval(cmd);
    
    
    %% Break curve fragment basing on geometric and semantic cues
    % load in computed curve fragment candidates
    input_file = [cem_src_path img_files(c).name(1:end-4) '.cem'];
    [CEM, edges, cfrags_idx] = load_contours(input_file);
    [~, edgemap, thetamap] = load_edg([edge_src_path img_files(c).name(1:end-4) '.edg']);
    
    img = imread([img_src_path img_files(c).name]);
    [h,w,~]= size(img);
    diag = sqrt(h^2+w^2);
    
    params.nbr_num_edges = max(round(nbr_num_edges*diag/params.diag_of_train), 5);
    params.diag_ratio = diag/ params.diag_of_train;

    % compute hsv space map
    hsv_img = rgb2hsv(img);
    
    % compute texton map
    tmap = assignTextons(fbRun(textonData.fb,rgb2gray(img)),textonData.tex);
    
    %%%%%%%%%%%%%%%%%%%%%% break curve fragments at conners and junctions
    tic;
    [new_cfrags, new_cfrags_idx, corner_pts, junction_pts] = contour_breaker_geom(CEM{2}, cfrags_idx, edgemap, params);
    toc;
    
    %%%%%%%%%%%%%%%%%%%%%% break curve fragments given semantic cues
    tic
    [new_cfrags, new_cfrags_idx, break_pts] = contour_breaker_semantic(new_cfrags, new_cfrags_idx, hsv_img, tmap, edgemap, params);
    toc
    
%     %%%%%%%%%%%%%%%%%%%%%%  save results of broken cfrags %%%%%%%%%%%%%%%%%%%%%%%
%     % display and save results
%     imshow(img, 'border', 'tight'); hold on;
%     draw_contours(new_cfrags);
%     if(~isempty(junction_pts))
%         plot(junction_pts(:,1)+1, junction_pts(:,2)+1, 'rx');
%     end
%     if(~isempty(corner_pts))
%         plot(corner_pts(:,1)+1, corner_pts(:,2)+1, 'gx');
%     end
%     if(~isempty(break_pts))
%         plot(break_pts(:,1)+1, break_pts(:,2)+1, 'yx');
%     end
%     hold off;
%     H = figure(1);
%     set(H, 'PaperPositionMode','auto')
%     print(H,'-dpng','-r0',[final_cem_path img_files(c).name(1:end-4) '.png']);
%     print_pdf([broken_cfrag_path img_files(c).name(1:end-4) '.pdf'])
%     tic
%     write_cem_fixed_edge_id([broken_cfrag_path img_files(c).name(1:end-4) '.cem'], new_cfrags, edges, h, w, new_cfrags_idx)
%     toc
%     keyboard;
    
    %% Merge curve fragments under graphical model
    tic
    [G, merged_cem, merged_cf_idx] = merge_cfrags_graphical_model(new_cfrags, new_cfrags_idx, edges, hsv_img, tmap, edgemap, params);
    toc
    
    
    %%%%%%%%%%%%%%%%%%  save results  %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % display and save results
    imshow(img, 'border', 'tight'); hold on;
    draw_contours(merged_cem);
    hold off;
    H = figure(1);
    set(H, 'PaperPositionMode','auto')
    print(H,'-dpng','-r0',[final_cem_path img_files(c).name(1:end-4) '.png']);
%     print_pdf([final_cem_path img_files(c).name(1:end-4) '.pdf'])
%     tic
    write_cem_fixed_edge_id([final_cem_path img_files(c).name(1:end-4) '.cem'], merged_cem, edges, h, w, merged_cf_idx)
%     toc
%     det_save_cemv([final_cem_path img_files(c).name(1:end-4) '.cemv'], new_cfrags);
        
% %% Rank of results cfrags and prune using logistic regressor        
%     disp('rank curve fragments');
%     P_vec = [];
%     for i = 1:length( merged_cem)
%         [bg_grad, sat_grad, hue_grad, abs_k, edge_sparsity, wigg, len, mean_conf] = curve_fragment_cues( merged_cem{i}, hsv_img, edgemap);
% 
%         p = 1 / (1 + exp(-([1, bg_grad, sat_grad, hue_grad, abs_k, edge_sparsity, wigg, len, mean_conf]-fmean_2)*beta_2'));
% %         p = len;
%         P_vec = [P_vec p];
%     end 
%     
%     % rank the cfrags
%     [~, sort_id] = sort(P_vec, 2, 'descend');
%      merged_cem =  merged_cem(sort_id); 
%     
%     top_cfrags =  merged_cem(1:top_K);
%     
% %     det_save_cemv([final_cem_path img_files(c).name(1:end-4) '_top_' num2str(top_K) '.cemv'], top_cfrags);
%  
%     imshow(rgb2gray(img), 'border', 'tight'); hold on;
%     draw_contours(top_cfrags);
%     hold off;
%     H = figure(1);
%     set(H, 'PaperPositionMode','auto')
% %     print(H,'-dpng','-r0',[final_cem_path img_files(c).name(1:end-4) '.png']);
%     print_pdf([best_F_vis_path img_files(c).name(1:end-4) '.pdf'])
    
end