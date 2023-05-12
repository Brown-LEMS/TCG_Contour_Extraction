%% Compute globalPb and hierarchical segmentation for an example image.
addpath ../TO-edge-detectorToolbox/;
addpath(fullfile(pwd,'lib'));

%% 1. compute globalPb on a BSDS image (5Gb of RAM required)
clear all; close all; clc;

img_path = '/media/guoy/Research/Datasets/VOCdevkit/VOC2007/JPEGImages/';
segObjGT_path = '/media/guoy/Research/Datasets/VOCdevkit/VOC2007/SegmentationObject/';
imgs = dir([segObjGT_path '*.png']);
dst_path = '/media/guoy/Research/Project_contour/MSEL_contour_extraction_full/Data/VOC2007/VOC2007_SE/edges/';
mkdir(dst_path);
th = 0.07;
r =1;
for i = 1:length(imgs)
    disp([num2str(i) '/' num2str(length(imgs))])
    img_name = imgs(i).name
    imgFile = [img_path img_name(1:end-4) '.jpg'];
    outFile = [dst_path img_name(1:end-4) '_gPb.mat'];
%     outimg = [img_path img_name(1:end-4) '_ucm.bmp'];
    outedg = [dst_path img_name(1:end-4) '.edg'];
%     gPb_orient = globalPb(imgFile, outFile);

    [gPb_map, TO_edge_map, gen_edge_map, gPb]=globalPb_subpixel_TO(imgFile,outFile, r, th);
    outImg = [dst_path img_name(1:end-4) '_bry.png'];
    imwrite(1-gPb, outImg, 'PNG');
%     [gPb_map, TO_edge_map, gen_edge_map]=globalPb_subpixel_TO(imgFile,outFile,1.0, 0.07, gbeta, mbeta);

    img = imread(imgFile);
    [h,w,~] = size(img);
    
    TO_edge_map(:,1) = TO_edge_map(:,1)-1;
    TO_edge_map(:,2) = TO_edge_map(:,2)-1;
    
    save_edg(outedg, TO_edge_map, [w,h])


end