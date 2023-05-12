clear all; close all;
addpath (genpath('/media/guoy/Research/Project_contour/MSEL_contour_extraction_full/util'))
img_src_path = '/media/guoy/Research/Project_contour/MSEL_contour_extraction_full/Data/Middlebury/trainingQ/';
edge_src_path = '/media/guoy/Research/Project_contour/MSEL_contour_extraction_full/Data/Middlebury/trainingQ_SE/';
frags_dst_path = '/media/guoy/Research/Project_contour/MSEL_contour_extraction_full/Data/Middlebury/trainingQ_SE/results_TCG/';
final_cem_path = frags_dst_path;
mkdir(frags_dst_path);

%% set parameters
params.vis = 0;
params.DP_gap_range = 10;
params.DP_angle_th = pi/4;
params.DP_contrast_th = 0.03;
params.shape_gap_range = 8;
params.shape_ori_range = pi/9;
params.BP_clen_th = 15;
params.BP_merge_angle_th = pi/9;
params.BP_nbr_num_edges = 30;
params.geom_merge_angle_th = pi/6;
params.nbr_num_edges = 30;
params.corner_angle_th = pi/6;
params.noise_prob_th = 0.1;
params.noise_len_th = 5;
% params.diag_of_train = 578.275; 
% load('training/SE_TCG_params.mat');

%% processing data
image_names{1} = 'Adirondack';
image_names{2} = 'ArtL';
image_names{3} = 'Jadeplant';
image_names{4} = 'Motorcycle';
image_names{5} = 'MotorcycleE';
image_names{6} = 'Piano';
image_names{7} = 'PianoL';
image_names{8} = 'Pipes';
image_names{9} = 'Playroom';
image_names{10} = 'Playtable';
image_names{11} = 'PlaytableP';
image_names{12} = 'Recycle';
image_names{13} = 'Shelves';
image_names{14} = 'Teddy';
image_names{15} = 'Vintage';

% img_files = dir([img_src_path img_ext]);

for c = 6%1:length(image_names)
    
    disp([num2str(c) '/'  num2str(length(image_names))]); 
    
    %% im0:  Compute cfrags from edges: curvelets, unambiguous cfrags, HPG cfrags
    input = [edge_src_path image_names{c} '/im0.edg'];
    output = [frags_dst_path image_names{c} 'im0.cem'];
    
    cmd = ['!util/dborl_compute_curve_frags ' input ' ' output];
    eval(cmd);
    
    img = imread([img_src_path image_names{c} '/im0.png']);
    [h,w,~]= size(img);
    diag = sqrt(h^2+w^2);
%     params.nbr_num_edges = max(round(nbr_num_edges*diag/params.diag_of_train), 5);
%     params.diag_ratio = diag/ params.diag_of_train;
    [~, edgemap, thetamap] = load_edg([edge_src_path image_names{c} '/im0.edg']);
    edgemap_soft0 = imread([edge_src_path image_names{c} '/im0_bry.png']);
    edgemap_soft0 = 1-im2double(edgemap_soft0);
%     ucm = imread([frags_dst_path '/' image_names{c} '/im0.png' ]);
%     ucm = 1-im2double(ucm);
   
    [CEM, edges, cfrags_idx] = load_contours(output);    


%     draw_contours(CEM{2},0,0);
    
    %%%%%%%%%%%%%%%%%% Break curve fragments at conners
    [new_cfrags, new_cfrags_idx, corner_pts] = contour_breaker_at_conner(CEM{2}, cfrags_idx, params, pi/15);
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
        print(H,'-dpdf','-r0',[frags_dst_path image_names{c} 'im0_gap.pdf']);
        cmd = ['!pdfcrop ' frags_dst_path image_names{c} 'im0_gap.pdf ' frags_dst_path image_names{c} 'im0_gap.pdf'];
        eval(cmd);
    end
    
    tic
    %%%%%%%%%%%%%%%%%% Break curve fragments at junctions
    [new_cfrags, new_cfrags_idx, edges] = break_contours_at_T_junctions(new_cfrags, new_cfrags_idx, edges);

    %%%%%%%%%%%%%%%%%% Prune noise curves
    [new_cfrags, new_cfrags_idx] = prune_noise_curves(new_cfrags, new_cfrags_idx, edgemap_soft0, params);
    
    %%%%%%%%%%%%%%%%%% Merge curve fragments under graphical model
    [G_test, new_cfrags, new_cfrags_idx] = merge_cfrags_graphical_model_geom(new_cfrags, new_cfrags_idx, edges, params);
    
    %%%%%%%%%%%%%%%%%% Classify junction types and merge at T-junction
    [G_test, new_cfrags, new_cfrags_idx, T_junctions] = classify_junction_type_wrt_graph_BP(new_cfrags, new_cfrags_idx, edges, params);
    
    %%%%%%%%%%%%%%%%%% Break curve fragments at conners
    [new_cfrags, new_cfrags_idx, corner_pts] = contour_breaker_at_conner(new_cfrags, new_cfrags_idx, params);
    
    %%%%%%%%%%%%%%%%%% Prune noise curves
    [new_cfrags, new_cfrags_idx] = prune_noise_curves(new_cfrags, new_cfrags_idx, edgemap_soft0, params);

    %%%%%%%%%%%%%%%%%% Merge curve fragments under graphical model
%     [G_test, new_cfrags, new_cfrags_idx, T_junctions] = merge_cfrags_graphical_model_geom(new_cfrags, new_cfrags_idx, edges, params, pi/18);
    toc
    
    if(params.vis)
        H = figure(2);
        imshow(img, 'border', 'tight'); 
        hold on;
        draw_contours(new_cfrags);
    %     if(~isempty(T_junctions))
    %         plot(T_junctions(:,1)+1, T_junctions(:,2)+1, 'rx');
    %     end
    % 
    %     for v =1:length(G_test.var)
    %         if(length(G_test.var(v).nbrs_fac)>=3)
    %             c_x = edges(G_test.var(v).actual_edge_id, 1)+1;
    %             c_y = edges(G_test.var(v).actual_edge_id, 2)+1;
    %             plot(c_x, c_y, 'gx');
    %         end
    %     end
        hold off;
        set(H, 'PaperPositionMode','auto')
        print(H,'-dpdf','-r0',[frags_dst_path image_names{c} 'im0.pdf']);
        cmd = ['!pdfcrop ' frags_dst_path image_names{c} 'im0.pdf ' frags_dst_path image_names{c} 'im0.pdf'];
        eval(cmd);
    end
    
    %%%%%%%%%%%%%%%%%%  save results  %%%%%%%%%%%%%%%%%%%%%%%%%%%
    write_cem_fixed_edge_id(output, new_cfrags, edges, h, w, new_cfrags_idx)
    %% im1: Compute cfrags from edges: curvelets, unambiguous cfrags, HPG cfrags
    input = [edge_src_path image_names{c} '/im1.edg'];
    output = [frags_dst_path image_names{c} 'im1.cem'];
    
    cmd = ['!util/dborl_compute_curve_frags ' input ' ' output];
    eval(cmd);
    
    img = imread([img_src_path image_names{c} '/im1.png']);
    [h,w,~]= size(img);
    
    [~, edgemap, ~] = load_edg([edge_src_path image_names{c} '/im1.edg']);
    edgemap_soft1 = imread([edge_src_path image_names{c} '/im1_bry.png']);
    edgemap_soft1 = 1-im2double(edgemap_soft1);    
    
    [CEM, edges, cfrags_idx] = load_contours(output);    
%     draw_contours(CEM{2},0,0);
    
    %%%%%%%%%%%%%%%%%% Break curve fragments at conners
    [new_cfrags, new_cfrags_idx, corner_pts] = contour_breaker_at_conner(CEM{2}, cfrags_idx, params, pi/15);
    if(params.vis)
        figure(1);
        imshow(img, 'border', 'tight'); 
        hold on;
        draw_contours(new_cfrags,0,0);
    end

    
    %%%%%%%%%%%%%%%%%%%%%% Fill large gaps
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
        print(H,'-dpdf','-r0',[frags_dst_path image_names{c} 'im1_gap.pdf']);
        cmd = ['!pdfcrop ' frags_dst_path image_names{c} 'im1_gap.pdf ' frags_dst_path image_names{c} 'im1_gap.pdf'];
        eval(cmd);
    end
    
    tic
    %%%%%%%%%%%%%%%%%% Break curve fragments at junctions
    [new_cfrags, new_cfrags_idx, edges] = break_contours_at_T_junctions(new_cfrags, new_cfrags_idx, edges);

    %%%%%%%%%%%%%%%%%% Prune noise curves
    [new_cfrags, new_cfrags_idx] = prune_noise_curves(new_cfrags, new_cfrags_idx, edgemap_soft1, params);
    
    %%%%%%%%%%%%%%%%%% Merge curve fragments under graphical model
    [G_test, new_cfrags, new_cfrags_idx, T_junctions] = merge_cfrags_graphical_model_geom(new_cfrags, new_cfrags_idx, edges, params);
    
    %%%%%%%%%%%%%%%%%% Classify junction types and merge at T-junction
    [G_test, new_cfrags, new_cfrags_idx, T_junctions] = classify_junction_type_wrt_graph_BP(new_cfrags, new_cfrags_idx, edges, params);
    
    %%%%%%%%%%%%%%%%%% Break curve fragments at conners
    [new_cfrags, new_cfrags_idx, corner_pts] = contour_breaker_at_conner(new_cfrags, new_cfrags_idx, params);
    
    %%%%%%%%%%%%%%%%%% Prune noise curves
    [new_cfrags, new_cfrags_idx] = prune_noise_curves(new_cfrags, new_cfrags_idx,edgemap_soft0, params);

%     %%%%%%%%%%%%%%%%%% Merge curve fragments under graphical model
%     [G_test, new_cfrags, new_cfrags_idx, T_junctions] = merge_cfrags_graphical_model_geom(new_cfrags, new_cfrags_idx, edges, params, pi/18);    
    toc
    
    if(params.vis)
        H = figure(2);
        imshow(img, 'border', 'tight'); 
        hold on;
        draw_contours(new_cfrags);
    %     if(~isempty(T_junctions))
    %         plot(T_junctions(:,1)+1, T_junctions(:,2)+1, 'rx');
    %     end
    %     for v =1:length(G_test.var)
    %         if(length(G_test.var(v).nbrs_fac)>=3)
    %             c_x = edges(G_test.var(v).actual_edge_id, 1)+1;
    %             c_y = edges(G_test.var(v).actual_edge_id, 2)+1;
    %             plot(c_x, c_y, 'gx');
    %         end
    %     end
    %     hold off;
        set(H, 'PaperPositionMode','auto')
        print(H,'-dpdf','-r0',[frags_dst_path image_names{c} 'im1.pdf']);
        cmd = ['!pdfcrop ' frags_dst_path image_names{c} 'im1.pdf ' frags_dst_path image_names{c} 'im1.pdf'];
        eval(cmd);
    end
    %%%%%%%%%%%%%%%%%%  save results  %%%%%%%%%%%%%%%%%%%%%%%%%%%
    write_cem_fixed_edge_id(output, new_cfrags, edges, h, w, new_cfrags_idx)
end
