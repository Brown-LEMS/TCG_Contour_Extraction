clear all; close all;
addpath(genpath('/media/guoy/Research/Project_contour/MSEL_contour_extraction_full/util'))
img_path = '/media/guoy/Research/Project_contour/MSEL_contour_extraction_full/Data/Middlebury/trainingQ/';
% cem_path = '/media/guoy/Research/Project_contour/MSEL_contour_extraction_full/Data/Middlebury/trainingQ_SE/results_Kovesi/';
cem_path = '/media/guoy/Research/Project_contour/MSEL_contour_extraction_full/Data/Middlebury/trainingQ_SE/results_Kokkinos/';
out_path = [cem_path 'vis/'];
mkdir(out_path);
% beta_src = '../training/';
% beta_prefix = 'TO_SEL_';
% input_beta_2 = load([beta_src beta_prefix 'beta_of_cues_for_seletion.txt']);
% fmean_2 = input_beta_2(1,:);
% fstd_2 = input_beta_2(2,:);
% beta_2 = input_beta_2(3,:);
% beta_2 = beta_2./fstd_2;
vary_num = [1 3 5 7 10 15 20 30 40 60 80 120 160 240 320 400 500 600 750];
num_th = 400
maxDist = 0.075; % use the ratio of img diagnal as a threshold in assignments
sumP_bry = [];
cntP_bry = [];
sumR_bry = [];
cntR_bry = [];
sumConsist_bry = [];
cntConsist_bry = [];

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

for i = 14:length(image_names)

    i
    % load ground truth edges
    im0_gt_bry = imread([img_path image_names{i} '/im_0_disp_bry.png']);
    im1_gt_bry = imread([img_path image_names{i} '/im_1_disp_bry.png']);
    % Load the CEM file
    load([cem_path image_names{i} '_results.mat'])
%     [CEM0, edges0, cfrags_idx0] = load_contours([cem_path image_names{i} '/im0_ucm.cem']);
    img0 = imread([img_path image_names{i} '/im0.png']);
    img0_mask = imread([img_path image_names{i} '/mask0nocc.png']);
    se = offsetstrel('ball',3,3);
    img0_mask = imdilate(img0_mask, se);
    img0_mask = (img0_mask>128);

%     [CEM1, edges1, cfrags_idx1] = load_contours([cem_path image_names{i} '/im1_ucm.cem']);
    img1 = imread([img_path image_names{i} '/im1.png']);
    img1_mask = imread([img_path image_names{i} '/mask1nocc.png']);
    se = offsetstrel('ball',3,3);
    img1_mask = imdilate(img1_mask, se);
    img1_mask = (img1_mask>128);
    [h,w,~] = size(img0);
    
%     [CEM01, edges01, cfrags_idx01] = load_contours([out_path image_names{i} 'im_0to1.cem']);
%     [CEM10, edges10, cfrags_idx10] = load_contours([out_path image_names{i} 'im_1to0.cem']);
    
%     % compute hsv space map
%     hsv_img0 = rgb2hsv(img0);
%     [~, edgemap0, thetamap0] = load_edg([cem_path image_names{i} '/im0.edg']);
%     hsv_img1 = rgb2hsv(img1);
%     [~, edgemap1, thetamap1] = load_edg([cem_path image_names{i} '/im1.edg']);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     cfrags0 = CEM0{2};
    P_vec = zeros(1,length(cfrags0));
    for c = 1:length(cfrags0)
        P_vec(c) = cfrags_prob0(c)*size(cfrags0{c},1);
    end 
%     P_vec = cfrags_prob0;
    P_vec(isnan(P_vec))=0;
    [P_vec, sort_id] = sort(P_vec, 2, 'descend');
    cfrags0 = cfrags0(sort_id);
    
    EdgeGroupMap_0 = convert_cfrags_to_EdgeGroupMap (cfrags0, h, w);
    EdgeGroupMap_0 =EdgeGroupMap_0.*img0_mask;

    H = figure(1);
    imshow(rgb2gray(img0));
    hold on;
    draw_contours(cfrags0(1:min(length(cfrags0), num_th)));
    hold off;
    outfile = [out_path image_names{i} '_view0.pdf'];
    set(H, 'PaperPositionMode','auto')
    print(H,'-dpdf','-r0',outfile);
    cmd = ['!pdfcrop ' outfile ' ' outfile];
    eval(cmd);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     cfrags1 = CEM1{2};
    P_vec = zeros(1,length(cfrags1));
    for c = 1:length(cfrags1)
        P_vec(c) = cfrags_prob1(c)*size(cfrags1{c},1);
    end 
%     P_vec = cfrags_prob1;
    P_vec(isnan(P_vec))=0;
    [P_vec, sort_id] = sort(P_vec, 2, 'descend');
    cfrags1 = cfrags1(sort_id);
    
    EdgeGroupMap_1 = convert_cfrags_to_EdgeGroupMap (cfrags1, h, w);
    EdgeGroupMap_1 =EdgeGroupMap_1.*img1_mask;

    H = figure(1);
    imshow(rgb2gray(img1));
    hold on;
    draw_contours(cfrags1(1:min(length(cfrags1), num_th)));
    hold off;
    outfile = [out_path image_names{i} '_view1.pdf'];
    set(H, 'PaperPositionMode','auto')
    print(H,'-dpdf','-r0',outfile);
    cmd = ['!pdfcrop ' outfile ' ' outfile];
    eval(cmd);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     cfrags_trans01 = CEM01{2};
    P_vec = zeros(1,length(cfrags_trans01));
    for c = 1:length(cfrags_trans01)
        P_vec(c) = cfrags_prob01(c)*size(cfrags_trans01{c},1);
    end 
%     P_vec = cfrags_prob01;
    P_vec(isnan(P_vec))=0;
    [P_vec, sort_id] = sort(P_vec, 2, 'descend');
    cfrags_trans01 = cfrags_trans01(sort_id);
    
    EdgeGroupMap_01 = convert_cfrags_to_EdgeGroupMap (cfrags_trans01, h, w);
    EdgeGroupMap_01 =EdgeGroupMap_01.*img1_mask;

    H = figure(1);
    imshow(rgb2gray(img1));
    hold on;
    draw_contours(cfrags_trans01(1:min(length(cfrags_trans01), num_th)));
    hold off;
    outfile = [out_path image_names{i} '_view01.pdf'];
    set(H, 'PaperPositionMode','auto')
    print(H,'-dpdf','-r0',outfile);
    cmd = ['!pdfcrop ' outfile ' ' outfile];
    eval(cmd);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     cfrags_trans10 = CEM10{2};
    P_vec = zeros(1,length(cfrags_trans10));
    for c = 1:length(cfrags_trans10)
        P_vec(c) = cfrags_prob10(c)*size(cfrags_trans10{c},1);
    end 
%     P_vec = cfrags_prob10;
    P_vec(isnan(P_vec))=0;
    [P_vec, sort_id] = sort(P_vec, 2, 'descend');
    cfrags_trans10 = cfrags_trans10(sort_id);
    
    EdgeGroupMap_10 = convert_cfrags_to_EdgeGroupMap (cfrags_trans10, h, w);
    EdgeGroupMap_10 =EdgeGroupMap_10.*img0_mask;

    H = figure(1);
    imshow(rgb2gray(img0));
    hold on;
    draw_contours(cfrags_trans10(1:min(length(cfrags_trans10), num_th)));
    hold off;
    outfile = [out_path image_names{i} '_view10.pdf'];
    set(H, 'PaperPositionMode','auto')
    print(H,'-dpdf','-r0',outfile);
    cmd = ['!pdfcrop ' outfile ' ' outfile];
    eval(cmd);
    
    %% evaluation
    sumP_vec = zeros(length(vary_num)+1, 1);
    cntP_vec = zeros(length(vary_num)+1, 1);
    sumR_vec = zeros(length(vary_num)+1, 1);
    cntR_vec = zeros(length(vary_num)+1, 1);
    sumConsist_vec = zeros(length(vary_num)+1, 2);
    cntConsist_vec = zeros(length(vary_num)+1, 2);
    
    cur_EdgeGroupMap_0 = EdgeGroupMap_0;
    cur_EdgeGroupMap_0 (cur_EdgeGroupMap_0>num_th) = 0;
    cur_EdgeGroupMap_1 = EdgeGroupMap_1;
    cur_EdgeGroupMap_1 (cur_EdgeGroupMap_1>num_th) = 0;
    cur_EdgeGroupMap_01 = EdgeGroupMap_01;
    cur_EdgeGroupMap_01 (cur_EdgeGroupMap_01>num_th) = 0;
    cur_EdgeGroupMap_10 = EdgeGroupMap_10;
    cur_EdgeGroupMap_10 (cur_EdgeGroupMap_10>num_th) = 0;


    tic
    [cntP, sumP, cntR, sumR, vis1, vis01] = matching_edge_grouping_maps(cur_EdgeGroupMap_1, cur_EdgeGroupMap_01, maxDist);
    toc
    sumConsist_vec(end,1) = sumP+sumR;
    cntConsist_vec(end,1) = cntP+cntR;
    %         sumConsist_vec(ii+1,1) = sumR;
    %         cntConsist_vec(ii+1,1) = cntR;

    tic
    [cntP, sumP, cntR, sumR, vis0, vis10] = matching_edge_grouping_maps(cur_EdgeGroupMap_0, cur_EdgeGroupMap_10, maxDist);
    toc
    sumConsist_vec(end,2) = sumP+sumR;
    cntConsist_vec(end,2) = cntP+cntR;
    %         sumConsist_vec(ii,2) = sumR;
    %         cntConsist_vec(ii,2) = cntR;
    
%     tic
%     [cntP, sumP, cntR, sumR, vis1, vis01] = matching_edge_grouping_maps(EdgeGroupMap_1, EdgeGroupMap_01, maxDist);
%     toc
%     sumConsist_vec(end,1) = sumP+sumR;
%     cntConsist_vec(end,1) = cntP+cntR;
% %     sumConsist_vec(ii+1,1) = sumR;
% %     cntConsist_vec(ii+1,1) = cntR;
% 
%     tic
%     [cntP, sumP, cntR, sumR, vis0, vis10] = matching_edge_grouping_maps(EdgeGroupMap_0, EdgeGroupMap_10, maxDist);
%     toc
%     sumConsist_vec(end,2) = sumP+sumR;
%     cntConsist_vec(end,2) = cntP+cntR;
% %     sumConsist_vec(ii+1,2) = sumR;
% %     cntConsist_vec(ii+1,2) = cntR;

    H = figure(1);
    imshow(vis0)
    outfile = [out_path image_names{i} '_match0.pdf'];
    set(H, 'PaperPositionMode','auto')
    print(H,'-dpdf','-r0',outfile);
    cmd = ['!pdfcrop ' outfile ' ' outfile];
    eval(cmd);
    
    imshow(vis10)
    outfile = [out_path image_names{i} '_match10.pdf'];
    set(H, 'PaperPositionMode','auto')
    print(H,'-dpdf','-r0',outfile);
    cmd = ['!pdfcrop ' outfile ' ' outfile];
    eval(cmd);
    
    imshow(vis1)
    outfile = [out_path image_names{i} '_match1.pdf'];
    set(H, 'PaperPositionMode','auto')
    print(H,'-dpdf','-r0',outfile);
    cmd = ['!pdfcrop ' outfile ' ' outfile];
    eval(cmd);

    imshow(vis01)
    outfile = [out_path image_names{i} '_match01.pdf'];
    set(H, 'PaperPositionMode','auto')
    print(H,'-dpdf','-r0',outfile);
    cmd = ['!pdfcrop ' outfile ' ' outfile];
    eval(cmd);
     

    C = cntConsist_vec(end,1) / sumConsist_vec(end,1)   

    
    

end


