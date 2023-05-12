clear all; close all
addpath ../mex_unix64/
addpath toolbox/


img_src_path = '/media/guoy/Research/Datasets/Junction_curve_structure_Dataset/Data/BSDS500_img/test/';
edg_dst_path = '/media/guoy/Research/Datasets/Junction_curve_structure_Dataset/Data/BSDS500_eval/test_GETO/edges/';
mkdir(edg_dst_path);

img_files = dir([img_src_path '*.jpg']);

interp=1;
sigma=1.5;
grad_th=1;

for i = 1:length(img_files)
    i
    imgfile=[img_src_path img_files(i).name];
    edgfile=[edg_dst_path img_files(i).name(1:end-4) '.edg'];
    
    % this is only for grey image
    cmd = ['!./third_order_edge_detector ' imgfile ' ' edgfile ' ' ...
        num2str(sigma) ' ' num2str(grad_th) ' ' num2str(interp)];
    
    eval(cmd)


end