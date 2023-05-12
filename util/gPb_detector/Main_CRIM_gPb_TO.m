%% Compute globalPb and hierarchical segmentation for an example image.
addpath(genpath('../../third_order'));
addpath(fullfile(pwd,'lib'));

%% 1. compute globalPb on a BSDS image (5Gb of RAM required)
clear all; close all; clc;

img_path = 'BMS_test_retrained/';
imgs = dir([img_path '*.jpg']);
dst_path = 'BMS_test_retrained/';
mkdir(dst_path);
load('beta_gPb_BMS.txt')
load('beta_mPb_BMS.txt')

mbeta = beta_mPb_BMS(2,:)./beta_mPb_BMS(1,:);
gbeta = beta_gPb_BMS(2,:)./beta_gPb_BMS(1,:);


for i = 1:length(imgs)
    disp([num2str(i) '/' num2str(length(imgs))])
    img_name = imgs(i).name;
    imgFile = [img_path img_name];
    outFile = [dst_path img_name(1:end-4) '_gPb.mat'];
%     outimg = [img_path img_name(1:end-4) '_ucm.bmp'];
    outedg = [dst_path img_name(1:end-4) '.edg'];
%     gPb_orient = globalPb(imgFile, outFile);

    [gPb_map, TO_edge_map, gen_edge_map]=globalPb_subpixel_TO(imgFile,outFile,1.0, 0.01, gbeta, mbeta);
    
    img = imread(imgFile);
    [h,w,~] = size(img);
    
    save_edg(outedg, TO_edge_map, [w,h])
%     %% 2. compute Hierarchical Regions
% 
%     % for boundaries
%     ucm = contours2ucm(gPb_orient, 'imageSize');
%     imwrite(ucm, outimg);

end