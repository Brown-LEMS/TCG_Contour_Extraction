clear all; close all;
addpath (genpath('/media/guoy/Research/Project_contour/MSEL_contour_extraction_full/util'))
img_src_path = '/media/guoy/Research/Project_contour/wireframe-master/data/v1.1/test/';
edge_src_path = '/media/guoy/Research/Project_contour/MSEL_contour_extraction_full/Data/wireframe/wireframe_SE/edges/';
frags_dst_path = '/media/guoy/Research/Project_contour/MSEL_contour_extraction_full/Data/wireframe/wireframe_SE/SE_TCG_wireframe/';
final_cem_path = frags_dst_path;
mkdir(frags_dst_path);
%% set parameters
% top_K = 60;   
params.vis = 1;
params.DP_gap_range = 15;
params.DP_angle_th = pi/4;
params.DP_contrast_th = 0.1;
params.shape_gap_range = 6;
params.shape_ori_range = pi/9;
params.BP_clen_th = 15;
params.BP_merge_angle_th = pi/9;
params.BP_nbr_num_edges = 20;
params.geom_merge_angle_th = pi/6;
params.nbr_num_edges = 20;
params.corner_angle_th = pi/6;
params.noise_prob_th = 0.1;
params.noise_len_th = 5;


beta_src = 'training/';
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

params.diag_of_train = 578.275;

%%
img_files = dir([img_src_path '*.jpg']);

for i = 1:length(img_files)
    
    disp([num2str(i) '/'  num2str(length(img_files))]); 
    
    %% im0:  Compute cfrags from edges: curvelets, unambiguous cfrags, HPG cfrags
    input = [edge_src_path img_files(i).name(1:end-4) '.edg'];
    output = [frags_dst_path img_files(i).name(1:end-4) '.cem'];
    output_jct = [frags_dst_path img_files(i).name(1:end-4) '.jct'];

    cmd = ['!util/dborl_compute_curve_frags ' input ' ' output];
    eval(cmd);
    
    img = imread([img_src_path img_files(i).name(1:end-4) '.jpg']);
%     img = imresize(img, 0.5);
    % compute hsv space map
    hsv_img = rgb2hsv(img);
    
%     % compute texton map
%     tmap = assignTextons(fbRun(textonData.fb,rgb2gray(img)),textonData.tex);
    [h,w,~]= size(img);
    [~, edgemap, thetamap] = load_edg(input);
    edgemap_soft0 = imread([edge_src_path img_files(i).name(1:end-4) '_bry.png']);
    edgemap_soft0 = 1-im2double(edgemap_soft0);
%     ucm = imread([frags_dst_path '/' img_files{c} '/im0.png' ]);
%     ucm = 1-im2double(ucm);
   
    [CEM, edges, cfrags_idx] = load_contours(output);    


%     draw_contours(CEM{2},0,0);
    
    %%%%%%%%%%%%%%%%%% Break curve fragments at conners
    [new_cfrags, new_cfrags_idx, corner_pts] = contour_breaker_at_conner(CEM{2}, cfrags_idx, params, pi/18);
    if(params.vis)
        figure(1);
        imshow(img, 'border', 'tight'); 
        hold on;
        draw_contours(new_cfrags,0,0);
    end

    %%%%%%%%%%%%%%%%%% Fill large gaps
    [edgemap_soft, thetamap] = imgradient(rgb2gray(img));
    edgemap_soft = edgemap_soft/max(edgemap_soft(:));
    thetamap = thetamap/180*pi;
    thetamap = wrapToPi(-thetamap+pi/2);
    tic;
    [new_cfrags, new_cfrags_idx, edges] = contour_fill_gaps_DP(new_cfrags, new_cfrags_idx, h,w, edges, edgemap_soft, thetamap, params);
%     [new_cfrags, new_cfrags_idx, edges] = contour_fill_gaps_DP_v2(new_cfrags, new_cfrags_idx, h,w, edges, edgemap_soft, thetamap, params);
%     [new_cfrags, new_cfrags_idx, edges] = contour_fill_gaps_DP_v2(CEM{2}, cfrags_idx, h,w, edges, edgemap_soft, edgemap, thetamap);
    toc;
    
    if(params.vis)
        hold off;
        H = figure(1);
        set(H, 'PaperPositionMode','auto')
        outimg = [frags_dst_path img_files(i).name(1:end-4) '_gap.pdf'];
        print(H,'-dpdf','-r0', outimg);
        cmd = ['!pdfcrop ' outimg ' ' outimg];
        eval(cmd);
    end
    
    tic
    %%%%%%%%%%%%%%%%%% Break curve fragments at T junctions
    [new_cfrags, new_cfrags_idx, edges] = break_contours_at_T_junctions(new_cfrags, new_cfrags_idx, edges);

    %%%%%%%%%%%%%%%%%% Prune noise curves
    [new_cfrags, new_cfrags_idx] = prune_noise_curves(new_cfrags, new_cfrags_idx, edgemap_soft0, params);
    
    %%%%%%%%%%%%%%%%%% Merge curve fragments under graphical model
    [G_test, new_cfrags, new_cfrags_idx] = merge_cfrags_graphical_model_geom(new_cfrags, new_cfrags_idx, edges, params);
    
    %%%%%%%%%%%%%%%%%% Classify junction types and merge at T-junction
    [G_test, new_cfrags, new_cfrags_idx, T_junctions] = classify_junction_type_wrt_graph_BP(new_cfrags, new_cfrags_idx, edges, params);
    
%     %%%%%%%%%%%%%%%%%% Break curve fragments at conners
%     [new_cfrags, new_cfrags_idx, corner_pts] = contour_breaker_at_conner(new_cfrags, new_cfrags_idx, params);
    
    %%%%%%%%%%%%%%%%%% Prune noise curves
    [new_cfrags, new_cfrags_idx] = prune_noise_curves(new_cfrags, new_cfrags_idx, edgemap_soft0, params);

    %%%%%%%%%%%%%%%%% Merge curve fragments under graphical model
    [G_test, new_cfrags, new_cfrags_idx] = merge_cfrags_graphical_model_geom(new_cfrags, new_cfrags_idx, edges, params, pi/18);

%     %%%%%%%%%%%%%%%%%% prune contour fragments
%     P_vec = [];
%     for j = 1:length( new_cfrags)
%         [bg_grad, sat_grad, hue_grad, abs_k, edge_sparsity, wigg, len, mean_conf] = curve_fragment_cues( new_cfrags{j}, hsv_img, edgemap);
% 
%         p = 1 / (1 + exp(-([1, bg_grad, sat_grad, hue_grad, abs_k, edge_sparsity, wigg, len, mean_conf]-fmean_2)*beta_2'));
%         P_vec = [P_vec p];
%     end 
%     
%     % rank the cfrags
%     [~, sort_id] = sort(P_vec, 2, 'descend');
%     new_cfrags =  new_cfrags(sort_id); 
%     new_cfrags =  new_cfrags(1:min(top_K, length(new_cfrags)));

    toc
    
    if(params.vis)
        H = figure(2);
        imshow(img, 'border', 'tight'); 
        hold on;
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
        outimg = [frags_dst_path img_files(i).name(1:end-4) '.pdf'];
        print(H,'-dpdf','-r0', outimg);
        cmd = ['!pdfcrop ' outimg ' ' outimg];
        eval(cmd);
    end
    %%%%%%%%%%%%%%%%%%  save results  %%%%%%%%%%%%%%%%%%%%%%%%%%%
    write_cem_fixed_edge_id(output, new_cfrags, edges, h, w, new_cfrags_idx)
    G_out = construct_fac_graph_from_curve_fragments_v2 (new_cfrags_idx, new_cfrags);
    write_jct(output_jct, G_out, size(edges, 1));
%     det_save_cemv([output 'v'], new_cfrags);
    close all;
end
