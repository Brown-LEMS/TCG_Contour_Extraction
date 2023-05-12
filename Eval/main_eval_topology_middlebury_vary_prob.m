clear all; close all;
addpath(genpath('/media/guoy/Research/Project_contour/MSEL_contour_extraction_full/util'))
img_path = '/media/guoy/Research/Project_contour/MSEL_contour_extraction_full/Data/Middlebury/trainingQ/';
% cem_path = '/media/guoy/Research/Project_contour/MSEL_contour_extraction_full/Data/Middlebury/trainingQ_SE/results_Kovesi/';
out_path = '/media/guoy/Research/Project_contour/MSEL_contour_extraction_full/Data/Middlebury/trainingQ_SE/results_Smin/';
eval_file = '/media/guoy/Research/Project_contour/MSEL_contour_extraction_full/Data/Middlebury/trainingQ_SE_Smin_topo_eval_vary_th.mat';
% beta_src = '../training/';
% beta_prefix = 'TO_SEL_';
% input_beta_2 = load([beta_src beta_prefix 'beta_of_cues_for_seletion.txt']);
% fmean_2 = input_beta_2(1,:);
% fstd_2 = input_beta_2(2,:);
% beta_2 = input_beta_2(3,:);
% beta_2 = beta_2./fstd_2;
vary_p = 0.85:-0.05:0;


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

for i = 1:length(image_names)

    i
    %load ground truth edges
    im0_gt_bry = imread([img_path image_names{i} '/im_0_disp_bry.png']);
    im1_gt_bry = imread([img_path image_names{i} '/im_1_disp_bry.png']);
    %Load the CEM file
    load([out_path image_names{i} '_results.mat'])
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
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     cfrags0 = CEM0{2};
    P_vec = cfrags_prob0;
    P_vec(isnan(P_vec))=0;
    [cfrags_prob0, sort_id] = sort(P_vec, 2, 'descend');
    cfrags0 = cfrags0(sort_id);
    
    EdgeGroupMap_0 = convert_cfrags_to_EdgeGroupMap (cfrags0, h, w);
    EdgeGroupMap_0 =EdgeGroupMap_0.*img0_mask;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     cfrags1 = CEM1{2};
    P_vec = cfrags_prob1;
    P_vec(isnan(P_vec))=0;
    [cfrags_prob1, sort_id] = sort(P_vec, 2, 'descend');
    cfrags1 = cfrags1(sort_id);
    
    EdgeGroupMap_1 = convert_cfrags_to_EdgeGroupMap (cfrags1, h, w);
    EdgeGroupMap_1 =EdgeGroupMap_1.*img1_mask;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     cfrags01 = CEM01{2};
    cfrags01 = cfrags_trans01;
    P_vec = cfrags_prob01;
    P_vec(isnan(P_vec))=0;
    [cfrags_prob01, sort_id] = sort(P_vec, 2, 'descend');
    cfrags01 = cfrags01(sort_id);
    
    EdgeGroupMap_01 = convert_cfrags_to_EdgeGroupMap (cfrags01, h, w);
    EdgeGroupMap_01 =EdgeGroupMap_01.*img1_mask;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     cfrags10 = CEM10{2};
    cfrags10 = cfrags_trans10;
    P_vec = cfrags_prob10;
    P_vec(isnan(P_vec))=0;
    [cfrags_prob10, sort_id] = sort(P_vec, 2, 'descend');
    cfrags10 = cfrags10(sort_id);
    
    EdgeGroupMap_10 = convert_cfrags_to_EdgeGroupMap (cfrags10, h, w);
    EdgeGroupMap_10 =EdgeGroupMap_10.*img0_mask;
    
    %% evaluation
    sumP_vec = zeros(length(vary_p), 1);
    cntP_vec = zeros(length(vary_p), 1);
    sumR_vec = zeros(length(vary_p), 1);
    cntR_vec = zeros(length(vary_p), 1);
    sumConsist_vec = zeros(length(vary_p), 2);
    cntConsist_vec = zeros(length(vary_p), 2);
    
    for ii =  1:length(vary_p)
        p_th = vary_p(ii)
        num_th0 = find(cfrags_prob0 > p_th);
        if(isempty(num_th0))
            num_th0 = 0;
        else
            num_th0 = num_th0(end);
        end
        num_th1 = find(cfrags_prob1 > p_th); 
        if(isempty(num_th1))
            num_th1 = 0;
        else
            num_th1 = num_th1(end);
        end
        num_th01 = find(cfrags_prob01 > p_th);
        if(isempty(num_th01))
            num_th01 = 0;
        else
            num_th01 = num_th01(end);
        end
        num_th10 = find(cfrags_prob10 > p_th); 
        if(isempty(num_th10))
            num_th10 = 0;
        else
            num_th10 = num_th10(end);
        end
        %%%%%%%%%%%%%%%%%% compare between transformed edges
        cur_EdgeGroupMap_0 = EdgeGroupMap_0;
        cur_EdgeGroupMap_0 (cur_EdgeGroupMap_0>num_th0) = 0;
        cur_EdgeGroupMap_1 = EdgeGroupMap_1;
        cur_EdgeGroupMap_1 (cur_EdgeGroupMap_1>num_th1) = 0;
        cur_EdgeGroupMap_01 = EdgeGroupMap_01;
        cur_EdgeGroupMap_01 (cur_EdgeGroupMap_01>num_th01) = 0;
        cur_EdgeGroupMap_10 = EdgeGroupMap_10;
        cur_EdgeGroupMap_10 (cur_EdgeGroupMap_10>num_th10) = 0;
%         % visulization
%         subplot(2,2,1);hold on;
%         title('view 0')
%         imshow(cur_EdgeGroupMap_0)
%         subplot(2,2,2);hold on;
%         title('view 10')
%         imshow(cur_EdgeGroupMap_10)
%         subplot(2,2,3);hold on;
%         title('view 1')
%         imshow(cur_EdgeGroupMap_1)
%         subplot(2,2,4);hold on;
%         title('view 01')
%         imshow(cur_EdgeGroupMap_01)
        
        tic
        [cntP, sumP, cntR, sumR] = matching_edge_grouping_maps(cur_EdgeGroupMap_1, cur_EdgeGroupMap_01, maxDist);
        toc
        sumConsist_vec(ii,1) = sumP+sumR;
        cntConsist_vec(ii,1) = cntP+cntR;

        
        tic
        [cntP, sumP, cntR, sumR] = matching_edge_grouping_maps(cur_EdgeGroupMap_0, cur_EdgeGroupMap_10, maxDist);
        toc
        sumConsist_vec(ii,2) = sumP+sumR;
        cntConsist_vec(ii,2) = cntP+cntR;
        
        %%%%%%%%%%%%%%%%%% compare to ground-truth edges
        tic
        [match1,match2] = correspondPixels(double(im0_gt_bry>0),double(cur_EdgeGroupMap_0>0),maxDist);
        toc
        sumP = sum(cur_EdgeGroupMap_0(:)>0);
        cntP = sum(match2(:)>0);
        sumR = sum(im0_gt_bry(:)>0);
        cntR = sum(match1(:)>0);
        
        sumP_vec(ii,1) = sumP;
        cntP_vec(ii,1) = cntP;
        sumR_vec(ii,1) = sumR;
        cntR_vec(ii,1) = cntR;    
        
%         tic
%         [cntP, sumP, cntR, sumR] = matching_edge_grouping_maps(double(im1_gt_bry>0), double(cur_EdgeGroupMap_1>0), maxDist);
%         toc
        tic
        [match1,match2] = correspondPixels(double(im1_gt_bry>0),double(cur_EdgeGroupMap_1>0),maxDist);
        toc
        sumP = sum(cur_EdgeGroupMap_1(:)>0);
        cntP = sum(match2(:)>0);
        sumR = sum(im1_gt_bry(:)>0);
        cntR = sum(match1(:)>0);
        
        sumP_vec(ii,2) = sumP;
        cntP_vec(ii,2) = cntP;
        sumR_vec(ii,2) = sumR;
        cntR_vec(ii,2) = cntR;     
        
%         P = cntP/sumP
%         R = cntR/sumR
        C = cntConsist_vec(ii,1) / sumConsist_vec(ii,1)
        R = cntR_vec(ii,1) /  sumR_vec(ii,1);
        C_score = 2*C*R/(C+R)

    end
    
    sumP_bry = [sumP_bry sumP_vec];
    cntP_bry = [cntP_bry cntP_vec];
    sumR_bry = [sumR_bry sumR_vec];
    cntR_bry = [cntR_bry cntR_vec];   
    sumConsist_bry = [sumConsist_bry sumConsist_vec];
    cntConsist_bry = [cntConsist_bry cntConsist_vec];
end


%% save eval results
% P_bry = cntP_bry./sumP_bry;
% R_bry = cntR_bry./sumR_bry;
% P_bry_mu = sum(cntP_bry, 2)./sum(sumP_bry, 2);
% R_bry_mu = sum(cntR_bry, 2)./sum(sumR_bry, 2);

P = sum(cntP_bry, 2)./sum(sumP_bry, 2);
R = sum(cntR_bry, 2)./sum(sumR_bry, 2);
Consist = sum(cntConsist_bry, 2)./sum(sumConsist_bry, 2);

save(eval_file, 'P', 'R', 'Consist', 'vary_p', 'sumP_bry', 'cntP_bry', 'sumR_bry', 'cntR_bry', 'sumConsist_bry', 'cntConsist_bry');

