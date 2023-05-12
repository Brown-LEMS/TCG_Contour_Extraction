clear all; close all;
addpath(genpath('/media/New_Volume/Research/Project_contour/third_order/'))
addpath(genpath('/media/New_Volume/Research/Print_PDF/'))
addpath ('lib');

gPb_path = 'CFGD/';
img_path = 'CFGD/';
th = 0.07;
r =1;

mat_files = dir([gPb_path '*.mat']);

for f = 1:length(mat_files)
    fname = mat_files(f).name(1:end-8)
    load([gPb_path mat_files(f).name]);
    I = imread([img_path fname '.jpg']);
    outFile = [gPb_path fname '.edg'];
    
%     gtheta = [1.5708    1.1781    0.7854    0.3927   0    2.7489    2.3562    1.9635];
%     [gPb_max, maxo] = max(gPb_orient,[],3);
%     ori_max = gtheta(maxo);
    edginfo = subpix_TO_correction_gPb( gPb_orient,I,th,r);
    edginfo(:,1) = edginfo(:,1)-1;
    edginfo(:,2) = edginfo(:,2)-1;
%     imshow(gPb_thin);
    save_edg(outFile, edginfo, [size(I,2) size(I,1)]);
end