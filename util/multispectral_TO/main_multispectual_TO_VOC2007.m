clear all; close all
addpath toolbox/

img_src_path = '/home/guoy/Desktop/Chris/silhouettes_frog/img/';
edg_dst_path = '/home/guoy/Desktop/Chris/silhouettes_frog/edges/';

img_files = dir([img_src_path '*.png']);

interp=1;
sigma=2;
T0=2;

for i = 1:length(img_files)
    i
    img=imread([img_src_path img_files(i).name]);
    [h,w,~]=size(img);
    dim = [w,h];
    T0_edge=multi_spect_TO_edge_detector(img, interp, sigma, T0, 'Lab', 0);
    
    saveedgname = [edg_dst_path img_files(i).name(1:end-4) '.edg'];
    save_edg(saveedgname, T0_edge, dim);
end