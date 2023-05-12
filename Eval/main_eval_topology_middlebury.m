clear all; close all;
addpath(genpath('/media/guoy/Research/Project_contour/MSEL_contour_extraction_full/util'))
img_path = '/media/guoy/Research/Datasets/MiddEval/trainingQ/';
cem_path = '/media/guoy/Research/Datasets/MiddEval/trainingQ_SE/';
out_path = '/media/guoy/Research/Datasets/MiddEval/trainingQ_SE/results_MSEL/';
eval_file = '/media/guoy/Research/Datasets/MiddEval/trainingQ_SE_MSEL_topo_eval.mat';
maxDist = 0.075; % use the ratio of img diagnal as a threshold in assignments
sumP_bry = [];
cntP_bry = [];
sumR_bry = [];
cntR_bry = [];


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
    se = offsetstrel('ball',3,3);
    img0_mask = imdilate(img0_mask, se);
    img0_mask = (img0_mask>128);

    [CEM1, edges1, cfrags_idx1] = load_contours([cem_path image_names{i} '/im1.cem']);
    img1 = imread([img_path image_names{i} '/im1.png']);
    img1_mask = imread([img_path image_names{i} '/mask1nocc.png']);
    se = offsetstrel('ball',3,3);
    img1_mask = imdilate(img1_mask, se);
    img1_mask = (img1_mask>128);
    [h,w,~] = size(img0);
    
    [CEM01, edges01, cfrags_idx01] = load_contours([out_path image_names{i} 'im_0to1.cem']);
    [CEM10, edges10, cfrags_idx10] = load_contours([out_path image_names{i} 'im_1to0.cem']);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cfrags = CEM0{2};

    EdgeGroupMap_0 = convert_cfrags_to_EdgeGroupMap (cfrags, h, w);
    EdgeGroupMap_0 =EdgeGroupMap_0.*img0_mask;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cfrags = CEM1{2};

    EdgeGroupMap_1 = convert_cfrags_to_EdgeGroupMap (cfrags, h, w);
    EdgeGroupMap_1 =EdgeGroupMap_1.*img1_mask;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cfrags = CEM01{2};

    EdgeGroupMap_01 = convert_cfrags_to_EdgeGroupMap (cfrags, h, w);
    EdgeGroupMap_01 =EdgeGroupMap_01.*img1_mask;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cfrags = CEM10{2};

    EdgeGroupMap_10 = convert_cfrags_to_EdgeGroupMap (cfrags, h, w);
    EdgeGroupMap_10 =EdgeGroupMap_10.*img0_mask;

    %% evaluation
    % boundary assignments 1 vs. 01
    %     [match1,match2] = correspondPixels(double(EdgeGroupMap_1>0),double(EdgeGroupMap_01>0),maxDist);
    tic
    [cntP, sumP, cntR, sumR] = matching_edge_grouping_maps(EdgeGroupMap_1, EdgeGroupMap_01, maxDist);
    toc

    sumP_bry = [sumP_bry sumP];
    cntP_bry = [cntP_bry cntP];
    sumR_bry = [sumR_bry sumR];
    cntR_bry = [cntR_bry cntR];    
    
    tic
    [cntP, sumP, cntR, sumR] = matching_edge_grouping_maps(EdgeGroupMap_0, EdgeGroupMap_10, maxDist);
    toc
    
    sumP_bry = [sumP_bry sumP];
    cntP_bry = [cntP_bry cntP];
    sumR_bry = [sumR_bry sumR];
    cntR_bry = [cntR_bry cntR];  
    
end


%% save eval results
P_bry = cntP_bry./sumP_bry;
R_bry = cntR_bry./sumR_bry;
P_bry_mu = sum(cntP_bry)/sum(sumP_bry);
R_bry_mu = sum(cntR_bry)/sum(sumR_bry);

save(eval_file, 'P_bry', 'R_bry', 'P_bry_mu', 'R_bry_mu');

