clear all; close all;
addpath(genpath('./'))


img_src_path = '/media/guoy/Research/Project_contour/MSEL_contour_extraction_full/Data/VOC2007/VOC2007_gPb_Maruthi/images/';
edge_src_path = '/media/guoy/Research/Project_contour/MSEL_contour_extraction_full/Data/VOC2007/VOC2007_gPb_Maruthi/edges/';

img_files = dir([img_src_path '*.png']);

for i = 1:length(img_files)
    i
    edge_file = [edge_src_path img_files(i).name(1:end-4) '_gpb.edg'];
    [edges, edgemap, thetamap] = load_edg(edge_file);
    I = imread([img_src_path img_files(i).name]);
    
    edges(:, 1:2) = edges(:, 1:2) -1;
    save_edg(edge_file, edges, [size(I,2) size(I,1)]);

end