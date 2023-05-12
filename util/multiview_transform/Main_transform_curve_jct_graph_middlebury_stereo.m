%Load the CEM file
clear all; close all;
addpath(genpath('/media/guoy/Research/Project_contour/MSEL_contour_extraction_full/util'))
img_path = '/media/guoy/Research/Project_contour/MSEL_contour_extraction_full/Data/Middlebury/trainingQ/';
cem_path = '/media/guoy/Research/Project_contour/MSEL_contour_extraction_full/Data/Middlebury/trainingQ_SE/results_TCG/';
prob_path = '/media/guoy/Research/Project_contour/MSEL_contour_extraction_full/Data/Middlebury/trainingQ_SE/';
out_path = '/media/guoy/Research/Project_contour/MSEL_contour_extraction_full/Data/Middlebury/trainingQ_SE/results_TCG/';
mkdir(out_path);
% is_ucm = 0;

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

for i = 6%1:length(image_names)

    i

    [CEM0, edges0, cfrags_idx0] = load_contours([cem_path image_names{i} 'im0.cem']);
    img0 = imread([img_path image_names{i} '/im0.png']);
    [CEM1, edges1, cfrags_idx1] = load_contours([cem_path image_names{i} 'im1.cem']);        
    img1 = imread([img_path image_names{i} '/im1.png']);
    [h,w,~] = size(img0);
    cfrags0 = CEM0{2};
    cfrags1 = CEM1{2};
    
    load([prob_path image_names{i} '/im0.mat']);
    E0 = E;
%         cfrags_prob0 = compute_cfrags_prob_given_probmap(cfrags0, E{1});
    cfrags_prob0 = compute_cfrags_prob_given_probmap(cfrags0, E0);
    load([prob_path image_names{i} '/im1.mat']);
    E1 = E;
%         cfrags_prob1 = compute_cfrags_prob_given_probmap (cfrags1, E{1});
    cfrags_prob1 = compute_cfrags_prob_given_probmap (cfrags1, E1);
    
    GT_disparity0 = readpfm([img_path image_names{i} '/disp0GT.pfm']);
    GT_disparity1 = readpfm([img_path image_names{i} '/disp1GT.pfm']);

    %%%%%%%%%%%%% transform contours given disparity
    [edges_trans01, cfrags_trans01, cfrags_idx_trans01, cfrags_prob01] = transform_Contours_given_GT_disparity_v2(cfrags_idx0, edges0, GT_disparity0, true, cfrags_prob0);
    [edges_trans10, cfrags_trans10, cfrags_idx_trans10, cfrags_prob10] = transform_Contours_given_GT_disparity_v2(cfrags_idx1, edges1, GT_disparity1, false, cfrags_prob1);
    %%%%%%%%%%%%% update corresponding boundary probability
    cfrags_prob01 = compute_cfrags_prob_given_probmap (cfrags_trans01, E1);
    cfrags_prob10 = compute_cfrags_prob_given_probmap (cfrags_trans10, E0);
    yy = round(edges_trans01(:,2));
    xx = round(edges_trans01(:,1));
    yy = min(yy, h);
    yy = max(yy, 1);
    xx = min(xx, w);
    xx = max(xx, 1);
    edges_trans01(:,4) = E1(sub2ind([h,w], yy, xx));
    yy = round(edges_trans10(:,2));
    xx = round(edges_trans10(:,1));
    yy = min(yy, h);
    yy = max(yy, 1);
    xx = min(xx, w);
    xx = max(xx, 1);
    edges_trans10(:,4) = E0(sub2ind([h,w], yy, xx));
    % Attention: cfrags_trans01 cfrags_trans10 has not been updated!
    
    
%     write_cem_fixed_edge_id([out_path image_names{i} 'im_0to1.cem'], cfrags_trans01, edges_trans01, h, w, cfrags_idx_trans01);
%     write_cem_fixed_edge_id([out_path image_names{i} 'im_1to0.cem'], cfrags_trans10, edges_trans10, h, w, cfrags_idx_trans10);
    save([out_path image_names{i} '_results.mat'], 'edges0', 'cfrags0', 'cfrags_idx0', 'cfrags_prob0',...
                                                        'edges1', 'cfrags1', 'cfrags_idx1', 'cfrags_prob1',...
                                                        'edges_trans01', 'cfrags_trans01', 'cfrags_idx_trans01', 'cfrags_prob01',...
                                                        'edges_trans10', 'cfrags_trans10', 'cfrags_idx_trans10', 'cfrags_prob10');
    
    %%%%%%%%%%%% visulize im1 to im0
    H=figure(1);
    imshow(img1); hold on;
    draw_contours(cfrags_trans01);
%     G_test = construct_fac_graph_from_curve_fragments (cfrags_idx_trans01, cfrags_trans01);
%     for v =1:length(G_test.var)
%         if(length(G_test.var(v).nbrs_fac)>=3)
%             c_x = round(edges_trans01(G_test.var(v).actual_edge_id, 1)+1);
%             c_y = round(edges_trans01(G_test.var(v).actual_edge_id, 2)+1);
%             plot(c_x, c_y, 'gx');
%         end
%     end
    hold off;

    set(H, 'PaperPositionMode','auto')
    out_name = [out_path image_names{i} 'im_0to1.pdf'];
    print(H,'-dpdf','-r0', out_name);
    cmd = ['!pdfcrop ' out_name ' ' out_name];
    eval(cmd);

    %%%%%%%%%%%% visulize im1 to im0
    H=figure(2);
    imshow(img0); hold on;
    draw_contours(cfrags_trans10);
%     G_test = construct_fac_graph_from_curve_fragments (cfrags_idx_trans10, cfrags_trans10);
%     for v =1:length(G_test.var)
%         if(length(G_test.var(v).nbrs_fac)>=3)
%             c_x = round(edges_trans10(G_test.var(v).actual_edge_id, 1)+1);
%             c_y = round(edges_trans10(G_test.var(v).actual_edge_id, 2)+1);
%             plot(c_x, c_y, 'gx');
%         end
%     end
    hold off;
    
    set(H, 'PaperPositionMode','auto')
    out_name = [out_path image_names{i} 'im_1to0.pdf'];
    print(H,'-dpdf','-r0', out_name);
    cmd = ['!pdfcrop ' out_name ' ' out_name];
    eval(cmd);
end