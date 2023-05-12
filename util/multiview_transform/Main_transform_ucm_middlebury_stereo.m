%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% the specialty of UCM: at each prob threshold, it make different
%%% close regions, the edge groupings are different
%%% A: project contour fragments at different levels separately 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;
addpath(genpath('/media/guoy/Research/Project_contour/MSEL_contour_extraction_full/util'))
img_path = '/media/guoy/Research/Project_contour/MSEL_contour_extraction_full/Data/Middlebury/trainingQ/';
cem_path = '/media/guoy/Research/Project_contour/MSEL_contour_extraction_full/Data/Middlebury/trainingQ_SE/results_ucm/';
prob_path = '/media/guoy/Research/Project_contour/MSEL_contour_extraction_full/Data/Middlebury/trainingQ_SE/';
out_path = '/media/guoy/Research/Project_contour/MSEL_contour_extraction_full/Data/Middlebury/trainingQ_SE/results_ucm/';
mkdir(out_path);

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
    load([cem_path image_names{i} 'im0_cfrag_maps.mat']);
    cfrags_maps0 = cfrags_maps;
    load([cem_path image_names{i} 'im1_cfrag_maps.mat']);
    cfrags_maps1 = cfrags_maps;
    
    cfrags_maps01 = cell(length(vary_p),1);
    cfrags_maps10 = cell(length(vary_p),1);
    
%     [CEM0, edges0, cfrags_idx0] = load_contours([cem_path image_names{i} '/im0_ucm.cem']);
    img0 = imread([img_path image_names{i} '/im0.png']);
%     [CEM1, edges1, cfrags_idx1] = load_contours([cem_path image_names{i} '/im1_ucm.cem']);
    img1 = imread([img_path image_names{i} '/im1.png']);
    [h,w,~] = size(img0);
    load([prob_path image_names{i} '/im0.mat']);
    ucms0 = ucms(:,:,1);
    ucms1 = ucms(:,:,1);
    E0 = E;
    load([prob_path image_names{i} '/im1.mat']);
    E1 = E;

    for p_i = 1:length(vary_p)
        disp([num2str(p_i) '/' num2str(length(vary_p))]);
        edges0 = cfrags_maps0.edges{p_i};
        cfrags0 = cfrags_maps0.cfrags{p_i};
        cfrags_idx0 = cfrags_maps0.cfrags_idx{p_i};
        edges1 = cfrags_maps1.edges{p_i};
        cfrags1 = cfrags_maps1.cfrags{p_i};
        cfrags_idx1 = cfrags_maps1.cfrags_idx{p_i};
        
        %%%%%%%% ucm is twice the size as cfrags
        cfrags_prob0 = compute_cfrags_prob_given_ucm (cfrags0, ucms0);
        cfrags_prob1 = compute_cfrags_prob_given_ucm (cfrags1, ucms1);

        GT_disparity0 = readpfm([img_path image_names{i} '/disp0GT.pfm']);
        GT_disparity1 = readpfm([img_path image_names{i} '/disp1GT.pfm']);

        %%%%%%%%%%%%% transform contours given disparity
        [edges_trans01, cfrags_trans01, cfrags_idx_trans01, cfrags_prob01] = transform_Contours_given_GT_disparity_v2(cfrags_idx0, edges0, GT_disparity0, true, cfrags_prob0);
        [edges_trans10, cfrags_trans10, cfrags_idx_trans10, cfrags_prob10] = transform_Contours_given_GT_disparity_v2(cfrags_idx1, edges1, GT_disparity1, false, cfrags_prob1);

        cfrags_maps01{p_i} = cfrags_trans01;
        cfrags_maps10{p_i} = cfrags_trans10;
%         %%%%%%%%%%%%% update corresponding boundary probability
%         cfrags_prob01 = compute_cfrags_prob_given_probmap (cfrags_trans01, E1);
%         cfrags_prob10 = compute_cfrags_prob_given_probmap (cfrags_trans10, E0);
%         yy = round(edges_trans01(:,2));
%         xx = round(edges_trans01(:,1));
%         yy = min(yy, h);
%         yy = max(yy, 1);
%         xx = min(xx, w);
%         xx = max(xx, 1);
%         edges_trans01(:,4) = E1(sub2ind([h,w], yy, xx));
%         yy = round(edges_trans10(:,2));
%         xx = round(edges_trans10(:,1));
%         yy = min(yy, h);
%         yy = max(yy, 1);
%         xx = min(xx, w);
%         xx = max(xx, 1);
%         edges_trans10(:,4) = E0(sub2ind([h,w], yy, xx));
%         % Attention: cfrags_trans01 cfrags_trans10 has not been updated!

    
    end
%     write_cem_fixed_edge_id([out_path image_names{i} 'im_0to1.cem'], cfrags_trans01, edges_trans01, h, w, cfrags_idx_trans01);
%     write_cem_fixed_edge_id([out_path image_names{i} 'im_1to0.cem'], cfrags_trans10, edges_trans10, h, w, cfrags_idx_trans10);
    
    cfrags_maps0 = cfrags_maps0.cfrags;
    cfrags_maps1 = cfrags_maps1.cfrags;
    save([out_path image_names{i} '_results.mat'], 'cfrags_maps0', 'cfrags_maps1', 'cfrags_maps01', 'cfrags_maps10');
    
   
end