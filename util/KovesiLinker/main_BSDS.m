clear all; close all;
% addpath(genpath('/media/New_Volume/Research/Project_contour/third_order/pb/'));
% addpath /media/New_Volume/Research//Print_PDF/
addpath '/media/guoy/Research/Datasets/Junction_curve_structure_Dataset/util/io/'
imgPath  = '/media/guoy/Research/Datasets/BSDS500/data/images/test/';
edge_path = '/media/guoy/Research/Datasets/Junction_curve_structure_Dataset/Data/BSDS500_eval/test_GETO/edges/';
% edgePath = '../../data/CFGD/Pb/';
outPath = '/media/guoy/Research/Datasets/Junction_curve_structure_Dataset/Data/BSDS500_eval/test_GETO/Kovesi/';
if ~exist(outPath,'dir')
    mkdir(outPath);
end
imgList= dir([imgPath '*.jpg']);
numImage = length(imgList);
for iI= 1:numImage
    tic
    im = imread([imgPath '/' imgList(iI).name]);
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
    
    [yy,xx,mag] = find(edgeim>0);
    Angle = zeros(length(xx),1); 
    for k=1:length(xx)
        Angle(k) = thetamap(yy(k),xx(k));
    end
    edginfo  = [xx yy Angle mag];
    
    % construct edge idx map
    idx_map = zeros(size(edgeim));
    for eid = 1:size(edginfo, 1)
        idx_map(edginfo(eid, 2), edginfo(eid, 1)) = eid;
    end
    % construct cfrags_idx, contours recording edge idx
    cfrags_idx = cell(size(edgelist2));
    for c = 1:length(edgelist2)
        cfrags_idx{c} = zeros(1, size(edgelist2{c}, 1));
        for ii=1:size(edgelist2{c}, 1)
            cfrags_idx{c}(ii) = idx_map(edgelist2{c}(ii,2), edgelist2{c}(ii,1));
        end

    end
    
    
    % here edgelist is matlab coordinates, need to change to cxx before
    % write_cem
    % change edges to cxx coodinates
    edginfo(:,1) = edginfo(:,1)-1;
    edginfo(:,2) = edginfo(:,2)-1;
    cem_file = [outPath '/' imgList(iI).name(1:end-4)  '.cem'];
%     write_cem(cem_file, edgelist2, size(im,1), size(im,2));
    write_cem_fixed_edge_id(cem_file, edgelist2, edginfo,  size(im,1), size(im,2), cfrags_idx)
%     smooth_cem(cem_file, cem_file);
    % Display the edgelists with random colours for each distinct edge 
%     % in figure 2
%     imshow(im, 'border', 'tight');
%     hold on;
%      draw_contours(edgelist2, 0, 1)
%     hold off;
%     filename = [outPath '/' imgList(iI).name(1:end-4)  '.png'];
% %     saveas(gcf,filename);
%     print(filename, '-dpng');
% %     print_pdf(filename);
%     toc
%     close all
end
