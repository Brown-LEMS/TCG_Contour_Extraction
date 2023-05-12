clear all; close all
addpath ../mex_unix64/
addpath toolbox/


img_src_path = '/media/guoy/Research/Datasets/MiddEval/trainingQ/';
edg_dst_path = '/media/guoy/Research/Datasets/MiddEval/trainingQ_GE/';
mkdir(edg_dst_path);

image_names{1} = 'Adirondack';
image_names{2} = 'ArtL';
image_names{3} = 'Jadeplant';
image_names{4} = 'Motorcycle';
image_names{5} = 'MotorcycleE';
image_names{6} = 'Piano';
image_names{7} = 'PianoL';
image_names{8} = 'Pipes';
image_names{9} = 'Playroom';
image_names{10} = 'Playtable';
image_names{11} = 'PlaytableP';
image_names{12} = 'Recycle';
image_names{13} = 'Shelves';
image_names{14} = 'Teddy';
image_names{15} = 'Vintage';

interp=1;
sigma=1.5;
grad_th=1;

for i = 1:length(image_names)
    i
    img=imread([img_src_path image_names{i} '/im0.png']);
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
    saveedgname = [edg_dst_path image_names{i} '/im0.edg'];    
    mkdir([edg_dst_path image_names{i}]);
    save_edg(saveedgname, TO_edge, dim);
    
    img=imread([img_src_path image_names{i} '/im1.png']);
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
    saveedgname = [edg_dst_path image_names{i} '/im1.edg'];    
    mkdir([edg_dst_path image_names{i}]);
    save_edg(saveedgname, TO_edge, dim);
    
end