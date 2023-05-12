%% Compute globalPb and hierarchical segmentation for an example image.

addpath(fullfile(pwd,'lib'));

%% 1. compute globalPb on a BSDS image (5Gb of RAM required)
clear all; close all; clc;

img_path = 'BMS_test/';
imgs = dir([img_path '*.jpg']);

for i = 2:length(imgs)
    img_name = imgs(i).name;
    imgFile = [img_path img_name];
    outFile = [img_path img_name(1:end-4) '_gPb.mat'];
    outimg = [img_path img_name(1:end-4) '_ucm.bmp'];

    gPb_orient = globalPb(imgFile, outFile);

    %% 2. compute Hierarchical Regions

%     % for boundaries
%     ucm = contours2ucm(gPb_orient, 'imageSize');
%     imwrite(ucm, outimg);

end