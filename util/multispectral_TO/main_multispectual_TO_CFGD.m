clear all; close all
addpath toolbox/

img_src_path = '../../Data/CFGD_img/';
edg_dst_path = '../../Data/TO_SEL_CFGD/edges/';

img_files = dir([img_src_path '*.jpg']);

interp=1;
sigma=1;
grad_th=1;

for i = 1:length(img_files)
    i
    img=imread([img_src_path img_files(i).name]);
    [h,w,~]=size(img);
    dim = [w,h];
    T0_edge=multi_spect_TO_edge_detector(img, interp, sigma, grad_th, 'Lab', 0);
    
    saveedgname = [edg_dst_path img_files(i).name(1:end-4) '.edg'];
    save_edg(saveedgname, T0_edge, dim);
end