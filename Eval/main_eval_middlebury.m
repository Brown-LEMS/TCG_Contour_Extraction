clear all; close all;
addpath(genpath('/media/guoy/Research/Project_contour/MSEL_contour_extraction_full/util'))
img_path = '/media/guoy/Research/Datasets/MiddEval/trainingQ/';
cem_path = '/media/guoy/Research/Datasets/MiddEval/trainingQ_COB/';
out_path = '/media/guoy/Research/Datasets/MiddEval/trainingQ_COB/results_MSEL/';
eval_file = '/media/guoy/Research/Datasets/MiddEval/trainingQ_COB_MSEL_eval.mat';
maxDist = 0.005; % use the ratio of img diagnal as a threshold in assignments
sumP_bry = [];
cntP_bry = [];
sumR_bry = [];
cntR_bry = [];
sumP_jct = [];
cntP_jct = [];
sumR_jct = [];
cntR_jct = [];



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

for i = 1:length(image_names)

    i
    %Load the CEM file
    [CEM0, edges0, cfrags_idx0] = load_contours([cem_path image_names{i} '/im0.cem']);
    img0 = imread([img_path image_names{i} '/im0.png']);
    img0_mask = imread([img_path image_names{i} '/mask0nocc.png']);
    img0_mask = img0_mask>128;

    [CEM1, edges1, cfrags_idx1] = load_contours([cem_path image_names{i} '/im1.cem']);
    img1 = imread([img_path image_names{i} '/im1.png']);
    img1_mask = imread([img_path image_names{i} '/mask1nocc.png']);
    img1_mask = img1_mask>128;
    [h,w,~] = size(img0);
    
    [CEM01, edges01, cfrags_idx01] = load_contours([out_path image_names{i} 'im_0to1.cem']);
    [CEM10, edges10, cfrags_idx10] = load_contours([out_path image_names{i} 'im_1to0.cem']);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    G = construct_fac_graph_from_curve_fragments (cfrags_idx0, CEM0{2});
    cfrags = CEM0{2};
    edges = edges0;
    jct_map_0 = zeros(h, w);
    for v =1:length(G.var)
        if(length(G.var(v).nbrs_fac)>=3)
            c_x = edges(G.var(v).actual_edge_id, 1)+1;
            c_y = edges(G.var(v).actual_edge_id, 2)+1;
            c_x = min(c_x, w);
            c_x = max(c_x, 1);
            c_y = min(c_y, h);
            c_y = max(c_y, 1);
            jct_map_0(round(c_y), round(c_x)) = 1;
        end
    end
%     jct_map_0 = jct_map_0.*img0_mask;

    LabelEdgemap_0 = convert_cem_to_EdgeGroupMap (cfrags, h, w);
%     LabelEdgemap_0 = LabelEdgemap_0.*img0_mask;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    G = construct_fac_graph_from_curve_fragments (cfrags_idx1, CEM1{2});
    cfrags = CEM1{2};
    edges = edges1;
    jct_map_1 = zeros(h, w);
    for v =1:length(G.var)
        if(length(G.var(v).nbrs_fac)>=3)
            c_x = edges(G.var(v).actual_edge_id, 1)+1;
            c_y = edges(G.var(v).actual_edge_id, 2)+1;
            c_x = min(c_x, w);
            c_x = max(c_x, 1);
            c_y = min(c_y, h);
            c_y = max(c_y, 1);
            jct_map_1(round(c_y), round(c_x)) = 1;
        end
    end
%     jct_map_1 = jct_map_1.*img1_mask;
    
    LabelEdgemap_1 = convert_cem_to_EdgeGroupMap (cfrags, h, w);
%     LabelEdgemap_1 = LabelEdgemap_1.*img1_mask;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    G = construct_fac_graph_from_curve_fragments (cfrags_idx01, CEM01{2});
    cfrags = CEM01{2};
    edges = edges01;
    jct_map_01 = zeros(h, w);
    for v =1:length(G.var)
        if(length(G.var(v).nbrs_fac)>=3)
            c_x = edges(G.var(v).actual_edge_id, 1)+1;
            c_y = edges(G.var(v).actual_edge_id, 2)+1;
            c_x = min(c_x, w);
            c_x = max(c_x, 1);
            c_y = min(c_y, h);
            c_y = max(c_y, 1);
            jct_map_01(round(c_y), round(c_x)) = 1;
        end
    end
%     jct_map_01 = jct_map_01.*img1_mask;
    
    LabelEdgemap_01 = convert_cem_to_EdgeGroupMap (cfrags, h, w);
%     LabelEdgemap_01 = LabelEdgemap_01.*img1_mask;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    G = construct_fac_graph_from_curve_fragments (cfrags_idx10, CEM10{2});
    cfrags = CEM10{2};
    edges = edges10;
    jct_map_10 = zeros(h, w);
    for v =1:length(G.var)
        if(length(G.var(v).nbrs_fac)>=3)
            c_x = edges(G.var(v).actual_edge_id, 1)+1;
            c_y = edges(G.var(v).actual_edge_id, 2)+1;
            c_x = min(c_x, w);
            c_x = max(c_x, 1);
            c_y = min(c_y, h);
            c_y = max(c_y, 1);
            jct_map_10(round(c_y), round(c_x)) = 1;
        end
    end
%     jct_map_10 = jct_map_10.*img0_mask;
    
    LabelEdgemap_10 = convert_cem_to_EdgeGroupMap (cfrags, h, w);
%     LabelEdgemap_10 = LabelEdgemap_10.*img0_mask;
    
    %% evaluation
    % boundary assignments 1 vs. 01
    [match1,match2] = correspondPixels(double(LabelEdgemap_1>0),double(LabelEdgemap_01>0),maxDist);
    sumP_bry = [sumP_bry sum(LabelEdgemap_1(:))];
    cntP_bry = [cntP_bry sum(match1(:)>0)];
    sumR_bry = [sumR_bry sum(LabelEdgemap_01(:))];
    cntR_bry = [cntR_bry sum(match2(:)>0)];    
    
    vis_bry_01 = zeros(h,w,3);
    vis_bry_01(:,:, 1) = LabelEdgemap_1>0;
    vis_bry_01(:,:, 2) = LabelEdgemap_01>0;
    vis_bry_01(:,:, 2) = max(vis_bry_01(:,:, 2), match1>0);
    vis_bry_01(:,:, 2) = max(vis_bry_01(:,:, 2), match2>0);
    vis_bry_01(:,:, 3) = max(vis_bry_01(:,:, 3), match1>0);
    vis_bry_01(:,:, 3) = max(vis_bry_01(:,:, 3), match2>0);
    H = figure(1);
    imshow(vis_bry_01);
    hold on;
    
    % jct assignments 1 vs. 01
    [match1,match2] = correspondPixels(double(jct_map_1>0),double(jct_map_01>0),maxDist);
    sumP_jct = [sumP_jct sum(jct_map_1(:))];
    cntP_jct = [cntP_jct sum(match1(:)>0)];
    sumR_jct = [sumR_jct sum(jct_map_01(:))];
    cntR_jct = [cntR_jct sum(match2(:)>0)];    
    
    [Y, X] = find(jct_map_1>0);
    plot(X,Y, 'rx');
    [Y, X] = find(jct_map_01>0);
    plot(X,Y, 'gx');
    hold off;
    set(H, 'PaperPositionMode','auto')
    out_name = [out_path image_names{i} 'im_0to1_match.png'];
    print(H,'-dpng','-r0', out_name);
    
    % boundary assignments 0 vs. 10
    [match1,match2] = correspondPixels(double(LabelEdgemap_0>0),double(LabelEdgemap_10>0),maxDist);
    sumP_bry = [sumP_bry sum(LabelEdgemap_0(:))];
    cntP_bry = [cntP_bry sum(match1(:)>0)];
    sumR_bry = [sumR_bry sum(LabelEdgemap_10(:))];
    cntR_bry = [cntR_bry sum(match2(:)>0)];    
    
    vis_bry_10 = zeros(h,w,3);
    vis_bry_10(:,:, 1) = LabelEdgemap_0>0;
    vis_bry_10(:,:, 2) = LabelEdgemap_10>0;
    vis_bry_10(:,:, 2) = max(vis_bry_10(:,:, 2), match1>0);
    vis_bry_10(:,:, 2) = max(vis_bry_10(:,:, 2), match2>0);
    vis_bry_10(:,:, 3) = max(vis_bry_10(:,:, 3), match1>0);
    vis_bry_10(:,:, 3) = max(vis_bry_10(:,:, 3), match2>0);
    H = figure(2);
    imshow(vis_bry_10);
    hold on;
    
    % jct assignments 0 vs. 10
    [match1,match2] = correspondPixels(double(jct_map_0>0),double(jct_map_10>0),maxDist);
    sumP_jct = [sumP_jct sum(jct_map_0(:))];
    cntP_jct = [cntP_jct sum(match1(:)>0)];
    sumR_jct = [sumR_jct sum(jct_map_10(:))];
    cntR_jct = [cntR_jct sum(match2(:)>0)];  
    
    [Y, X] = find(jct_map_0>0);
    plot(X,Y, 'rx');
    [Y, X] = find(jct_map_10>0);
    plot(X,Y, 'gx');
    hold off;
    set(H, 'PaperPositionMode','auto')
    out_name = [out_path image_names{i} 'im_1to0_match.png'];
    print(H,'-dpng','-r0', out_name);
    
end


%% save eval results
P_bry = cntP_bry./sumP_bry;
R_bry = cntR_bry./sumR_bry;
P_jct = cntP_jct./sumP_jct;
R_jct = cntR_jct./sumR_jct;
P_bry_mu = sum(cntP_bry)/sum(sumP_bry);
R_bry_mu = sum(cntR_bry)/sum(sumR_bry);
P_jct_mu = sum(cntP_jct)/sum(sumP_jct);
R_jct_mu = sum(cntR_jct)/sum(sumR_jct);
save(eval_file, 'P_bry', 'R_bry', 'P_jct', 'R_jct', 'P_bry_mu', 'R_bry_mu', 'P_jct_mu', 'R_jct_mu');

