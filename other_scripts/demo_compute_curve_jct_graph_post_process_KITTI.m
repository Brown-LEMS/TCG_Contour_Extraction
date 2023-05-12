clear all; close all;
addpath (genpath('util'))
img_src_path = '/gpfs/scratch/yg13/Datasets/KITTI/data_stereo_flow/training/colored_1/';
edge_src_path = '/gpfs/scratch/yg13/Datasets/KITTI/data_stereo_flow/training/colored_1_SE/';
frags_dst_path = '/gpfs/scratch/yg13/Datasets/KITTI/data_stereo_flow/training/colored_1_SE/';
final_cem_path = frags_dst_path;
mkdir(frags_dst_path);
beta_src = 'training/';
img_ext = '*.png';
%% load in beta
prefix = 'TO_SEL_';

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
params.merge_th = 0.3;
params.merge_th_geom = 0.3;
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

%% processing data
img_files = dir([img_src_path img_ext]);

for c = 1:length(img_files)
    
    disp([num2str(c) '/'  num2str(length(img_files))]); 
    
    %% im0:  Compute cfrags from edges: curvelets, unambiguous cfrags, HPG cfrags
    input = [edge_src_path img_files(c).name(1:end-4) '.edg'];
    output = [frags_dst_path img_files(c).name(1:end-4) '.cem'];
    
    cmd = ['!util/dborl_compute_curve_frags ' input ' ' output];
%     cmd = ['!util/dborl_compute_curve_frags_pos_uncertain05 ' input ' ' output];
    eval(cmd);
    
    
    %%%%%%%%%%%%%%%%%% Post Process introduce more junctions
    img = imread([img_src_path img_files(c).name]);
    [h,w,~]= size(img);
    diag = sqrt(h^2+w^2);
    params.nbr_num_edges = max(round(nbr_num_edges*diag/params.diag_of_train), 5);
    params.diag_ratio = diag/ params.diag_of_train;
    [~, edgemap, ~] = load_edg(input);
    
    % compute hsv space map
    hsv_img = rgb2hsv(img);
    
    % compute texton map
    tmap = assignTextons(fbRun(textonData.fb,rgb2gray(img)),textonData.tex);
   
    [CEM, edges, cfrags_idx] = load_contours(output);    
    tic;
    [new_cfrags, new_cfrags_idx, T_junctions] = contour_introduce_jct(CEM{2}, cfrags_idx, h,w, edges);
    toc;
    
    %%%%%%%%%%%%%%%%%%%%%% break curve fragments given semantic cues
    tic
    [new_cfrags, new_cfrags_idx, break_pts] = contour_breaker_semantic(new_cfrags, new_cfrags_idx, hsv_img, tmap, edgemap, params);
    toc
    
    
    %%%%%%%%%%%%%%%%%% Merge curve fragments under graphical model
    tic
    [G_test, new_cfrags, new_cfrags_idx, T_junctions] = merge_cfrags_graphical_model(new_cfrags, new_cfrags_idx, edges, hsv_img, tmap, edgemap, params);
    toc

    %%%%%%%%%%%%%%%%%  visulization %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     G_test = construct_fac_graph_from_curve_fragments (new_cfrags_idx, new_cfrags);
    H = figure(1);
    imshow(img, 'border', 'tight'); hold on;
    draw_contours(new_cfrags);
    if(~isempty(T_junctions))
        plot(T_junctions(:,1)+1, T_junctions(:,2)+1, 'rx');
    end
    for v =1:length(G_test.var)
        if(length(G_test.var(v).nbrs_fac)>=3)
            c_x = edges(G_test.var(v).actual_edge_id, 1)+1;
            c_y = edges(G_test.var(v).actual_edge_id, 2)+1;
            plot(c_x, c_y, 'gx');
        end
    end
    hold off;
    set(H, 'PaperPositionMode','auto')
    outImgFile = [frags_dst_path img_files(c).name(1:end-4) '.pdf'];
    print(H,'-dpdf','-r0', outImgFile);
    cmd = ['!pdfcrop ' outImgFile ' ' outImgFile];
    eval(cmd);
    
    %%%%%%%%%%%%%%%%%%  save results  %%%%%%%%%%%%%%%%%%%%%%%%%%%
    write_cem_fixed_edge_id(output, new_cfrags, edges, h, w, new_cfrags_idx)
    
    
    
end