%Load the CEM file
clear all; close all;
addpath(genpath('/media/guoy/Research/Project_contour/MSEL_contour_extraction_full/util'))
img_path = '/media/guoy/Research/Datasets/Drapery_Datasets/Cloth1/';
cem_path = '/media/guoy/Research/Datasets/Drapery_Datasets/Cloth1_GE/MSEL/';
eval_file = '/media/guoy/Research/Datasets/Drapery_Datasets/Cloth1_GE_MSEL_eval.mat';
maxDist = 0.01; % use the ratio of img diagnal as a threshold in assignments
sumP_bry = [];
cntP_bry = [];
sumR_bry = [];
cntR_bry = [];
sumP_jct = [];
cntP_jct = [];
sumR_jct = [];
cntR_jct = [];

for i = 1:23
    figure(1)
    id0 = i-1;
    id1 = i;
    [CEM0, edges, cfrags_idx] = load_contours([cem_path 'Camera' num2str(id0) '.cem']);
    img0 = imread([img_path 'Camera' num2str(id0) '.png']);
    [h,w,~] = size(img0);
    imshow(img0); hold on;
    draw_contours(CEM0{2});
    G_trans = construct_fac_graph_from_curve_fragments (cfrags_idx, CEM0{2});
    count = 0;
    for v =1:length(G_trans.var)
        if(length(G_trans.var(v).nbrs_fac)>=3)
            c_x = edges(G_trans.var(v).actual_edge_id, 1)+1;
            c_y = edges(G_trans.var(v).actual_edge_id, 2)+1;
            plot(c_x, c_y, 'gx');
            count = count+1;
        end
    end
    hold off;
    disp(['number of original junctions: ' num2str(count)]);


    %Get Corresponding Ground Truth Contours
    %Arguments: Image ID1 
    %			Image ID2
    %			CEM data just loaded
    %			The path of dataset
    %Output: A Cell with all the corresponding curves in the second image.
%     correspContours = GetGroundTruthCorrespondingCurves(id0,id1,CEM0,img_path);
    [edges_trans, cfrags_trans, cfrags_idx_trans] = transform_Contours_given_GT_depth_v2(id0, id1, cfrags_idx, edges, img_path);
    
    H=figure(2);
    img1 = imread([img_path 'Camera' num2str(id1) '.png']);
    imshow(img1); hold on;
    draw_contours(cfrags_trans);
    G_trans = construct_fac_graph_from_curve_fragments (cfrags_idx_trans, cfrags_trans);
    count = 0;
    for v =1:length(G_trans.var)
        if(length(G_trans.var(v).nbrs_fac)>=3)
            c_x = edges_trans(G_trans.var(v).actual_edge_id, 1)+1;
            c_y = edges_trans(G_trans.var(v).actual_edge_id, 2)+1;
            plot(c_x, c_y, 'gx');
            jct_map_trans(round(c_x), round(c_y)) = 1;
            count = count+1;
        end
    end
    hold off;
    disp(['number of junction transformed: ' num2str(count)]);

    set(H, 'PaperPositionMode','auto')
    out_name = [cem_path 'Camera_' num2str(id0) 'to' num2str(id1) '.pdf'];
    print(H,'-dpdf','-r0', out_name);
    cmd = ['!pdfcrop ' out_name ' ' out_name];
    eval(cmd);

    
    %% %%%%%%%%%%%%%%%%%%%% evaluation: junction PR, bry PR %%%%%%%%%%%%%%%%%%%%%%%
    
    % transformed maps
    jct_map_trans = zeros(h, w);
    for v =1:length(G_trans.var)
        if(length(G_trans.var(v).nbrs_fac)>=3)
            c_x = edges_trans(G_trans.var(v).actual_edge_id, 1)+1;
            c_y = edges_trans(G_trans.var(v).actual_edge_id, 2)+1;
            jct_map_trans(round(c_y), round(c_x)) = 1;
        end
    end
    bry_map_trans = zeros(h, w);
    for c = 1:length(cfrags_trans)
        x_coords = round(cfrags_trans{c}(:,1))+1;
        y_coords = round(cfrags_trans{c}(:,2))+1;
        x_coords = min(x_coords, w);
        x_coords = max(x_coords, 1);
        y_coords = min(y_coords, h);
        y_coords = max(y_coords, 1);
        bry_map_trans(sub2ind([h,w], y_coords, x_coords)) = 1;
    end
    
    % map 1
    [CEM1, edges_1, cfrags_idx_1] = load_contours([cem_path 'Camera' num2str(id1) '.cem']);
    cfrags_1 = CEM1{2};
    G_1 = construct_fac_graph_from_curve_fragments (cfrags_idx_1, cfrags_1);
    jct_map_1 = zeros(h, w);
    for v =1:length(G_1.var)
        if(length(G_1.var(v).nbrs_fac)>=3)
            c_x = edges_1(G_1.var(v).actual_edge_id, 1)+1;
            c_y = edges_1(G_1.var(v).actual_edge_id, 2)+1;
            jct_map_1(round(c_y), round(c_x)) = 1;
        end
    end
    bry_map_1 = zeros(h, w);
    for c = 1:length(cfrags_1)
        x_coords = round(cfrags_1{c}(:,1))+1;
        y_coords = round(cfrags_1{c}(:,2))+1;
        x_coords = min(x_coords, w);
        x_coords = max(x_coords, 1);
        y_coords = min(y_coords, h);
        y_coords = max(y_coords, 1);
        bry_map_1(sub2ind([h,w], y_coords, x_coords)) = 1;
    end
    
    % boundary assignments
    [match1,match2] = correspondPixels(double(bry_map_1>0),double(bry_map_trans>0),maxDist);
    sumP_bry = [sumP_bry sum(bry_map_1(:))];
    cntP_bry = [cntP_bry sum(match1(:)>0)];
    sumR_bry = [sumR_bry sum(bry_map_trans(:))];
    cntR_bry = [cntR_bry sum(match2(:)>0)];    
    
    % jct assignments
    [match1,match2] = correspondPixels(double(jct_map_1>0),double(jct_map_trans>0),maxDist);
    sumP_jct = [sumP_jct sum(jct_map_1(:))];
    cntP_jct = [cntP_jct sum(match1(:)>0)];
    sumR_jct = [sumR_jct sum(jct_map_trans(:))];
    cntR_jct = [cntR_jct sum(match2(:)>0)];        
%     keyboard
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

% H = figure(3); hold on;
% plot(R_bry, P_bry, 'rx');
% plot(R_bry_mu, P_bry_mu, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
% hold off;
% title('boundary evaluation')
% axis([0 1 0 1]);
% xlabel('Recall')
% ylabel('Precision')
% grid on;
% set(H, 'PaperPositionMode','auto')
% out_name = [eval_file(1:end-4) '_bry.png'];
% print(H,'-dpng','-r0', out_name);
% 
% 
% 
% H=figure(4); hold on;
% plot(R_jct, P_jct, 'rx');
% plot(R_jct_mu, P_jct_mu, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
% hold off;
% title('junction evaluation')
% axis([0 1 0 1]);
% xlabel('Recall')
% ylabel('Precision')
% grid on;
% set(H, 'PaperPositionMode','auto')
% out_name = [eval_file(1:end-4) '_jct.png'];
% print(H,'-dpng','-r0', out_name);

