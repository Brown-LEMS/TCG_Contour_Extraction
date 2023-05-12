%% Compute globalPb and hierarchical segmentation for an example image.
addpath ../TO-edge-detectorToolbox/;
addpath(fullfile(pwd,'lib'));

%% 1. compute globalPb on a BSDS image (5Gb of RAM required)
clear all; close all; clc;

img_path = '../../Data/CFGD_img/';
imgs = dir([img_path '*.jpg']);
dst_path = '../../Data/TO_SEL_CFGD/edges/';
mkdir(dst_path);

load('beta_mPb_CFGD.txt');
mbeta = beta_mPb_CFGD(2,:)./beta_mPb_CFGD(1,:);
load('beta_gPb_CFGD.txt');
gbeta = beta_gPb_CFGD(2,:)./beta_gPb_CFGD(1,:);

for i = 1:length(imgs)
    disp([num2str(i) '/' num2str(length(imgs))])
    img_name = imgs(i).name
    imgFile = [img_path img_name];
    outFile = [dst_path img_name(1:end-4) '_gPb.mat'];
%     outimg = [img_path img_name(1:end-4) '_ucm.bmp'];
    outedg = [dst_path img_name(1:end-4) '.edg'];
%     gPb_orient = globalPb(imgFile, outFile);

    [gPb_map, TO_edge_map, gen_edge_map]=globalPb_subpixel_TO(imgFile,outFile,1.0, 0.07);

%     [gPb_map, TO_edge_map, gen_edge_map]=globalPb_subpixel_TO(imgFile,outFile,1.0, 0.07, gbeta, mbeta);

    img = imread(imgFile);
    [h,w,~] = size(img);
    
    TO_edge_map(:,1) = TO_edge_map(:,1)-1;
    TO_edge_map(:,2) = TO_edge_map(:,2)-1;
    
    save_edg(outedg, TO_edge_map, [w,h])


end