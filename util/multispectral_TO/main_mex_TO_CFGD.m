clear all; close all
addpath ../mex_unix64/
addpath toolbox/


img_src_path = '/media/guoy/Research/Datasets/Drapery_Datasets/Cloth1/';
edg_dst_path = '/media/guoy/Research/Datasets/Drapery_Datasets/Cloth1/GE/edges/';
mkdir(edg_dst_path);

img_files = dir([img_src_path '*.png']);

interp=1;
sigma=1.5;
grad_th=1.5;

for i = 1:length(img_files)
    i
    img=imread([img_src_path img_files(i).name]);
    [h,w,~]=size(img);
    dim = [w,h];
%     img = repmat(img, [1,1,3]);
%     TO_edge=multi_spect_TO_edge_detector(img, interp, sigma, grad_th, 'Lab', 0);
    if(size(img,3)==3)
        TO_edge = mex_third_order_color_edge_detector(img,1, sigma, grad_th, interp, w, h);
    else
        TO_edge = mex_third_order_edge_detector(img, sigma, grad_th, interp);
    end
    
    % convert to c++ coordinates
    TO_edge(:, 1:2) = TO_edge(:, 1:2)-1;
    saveedgname = [edg_dst_path img_files(i).name(1:end-4) '.edg'];    
    save_edg(saveedgname, TO_edge, dim);
end