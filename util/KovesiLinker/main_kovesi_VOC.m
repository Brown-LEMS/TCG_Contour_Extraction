clear all; close all;
% addpath(genpath('/media/New_Volume/Research/Project_contour/third_order/pb/'));
addpath(genpath('../'));
imgPath  = '/media/guoy/Research/Datasets/VOCdevkit/VOC2007/JPEGImages/';
seg_src_path = '/media/guoy/Research/Datasets/VOCdevkit/VOC2007/SegmentationObject/';
edge_path = '/media/guoy/Research/Project_contour/MSEL_contour_extraction_full/Data/VOC2007/VOC2007_SE/edges/';
% edgePath = '../../data/CFGD/Pb/';
outPath = '/media/guoy/Research/Project_contour/MSEL_contour_extraction_full/Data/VOC2007/VOC2007_SE/SE_Kovesi_VOC2007/';
if ~exist(outPath,'dir')
    mkdir(outPath);
end
top_K = 60;   

imgList= dir([seg_src_path '*.png']);
numImage = length(imgList);
for iI= 1:numImage
    tic
    im = imread([imgPath '/' imgList(iI).name(1:end-4) '.jpg']);
%     if (size(im,3)==3)
%         im= rgb2gray(im);
%     end
    img_name = imgList(iI).name(1:end-4)
    
    [~, edgeim, thetamap] = load_edg([edge_path img_name '.edg']);


    % Link edge pixels together into lists of sequential edge points, one
    % list for each edge contour. A contour/edgelist starts/stops at an 
    % ending or a junction with another contour/edgelist.
    % Here we discard contours less than 18 pixels long.

    [edgelist, labelededgeim] = edgelink((edgeim>0), 5);
    edgelist2 = edgelist;
    
    for c = 1:length(edgelist)
        cur_c = zeros(size(edgelist{c},1), 4);
        cur_c(:,1) = edgelist{c}(:,2);
        cur_c(:,2) = edgelist{c}(:,1);
        for e = 1:size(cur_c,1)
            cur_c(e,3) = thetamap(cur_c(e,2),cur_c(e,1));
            cur_c(e,4) = edgeim(cur_c(e,2),cur_c(e,1));
        end
        edgelist2{c} = cur_c;
    end
    
    
    disp('rank curve fragments');
    new_cfrags = edgelist2;
    P_vec = [];
    for c = 1:length(new_cfrags)
        p = contour_length_mex(new_cfrags{c}');
        P_vec = [P_vec p];
    end 
    
    % rank the cfrags
    [~, sort_id] = sort(P_vec, 2, 'descend');
    new_cfrags = new_cfrags(sort_id);

    % here edgelist is matlab coordinates, need to change to cxx before
%     % write_cem
    cem_file = [outPath '/' imgList(iI).name(1:end-4)  '.cem'];
    write_cem(cem_file, new_cfrags, size(im,1), size(im,2));
    det_save_cemv([cem_file 'v'], new_cfrags(1:top_K));

%     smooth_cem(cem_file, cem_file);
    % Display the edgelists with random colours for each distinct edge 
    % in figure 2
    imshow(im, 'border', 'tight');
    hold on;
    draw_contours(new_cfrags(1:top_K))
    hold off;
    filename = [outPath '/' imgList(iI).name(1:end-4)  '.png'];
%     saveas(gcf,filename);
    print(filename, '-dpng');
%     print_pdf(filename);
    toc
    close all
end
