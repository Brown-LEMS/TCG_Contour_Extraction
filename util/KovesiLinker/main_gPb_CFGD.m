clear all; close all;
% addpath(genpath('/media/New_Volume/Research/Project_contour/third_order/pb/'));
addpath /media/New_Volume/Research/Print_PDF/
imgPath  = '../../data/CFGD/';
gPb_path = '/media/New_Volume/Research/Project_contour/BSR/grouping/CFGD/';
% edgePath = '../../data/CFGD/Pb/';
outPath = '../../data/CFGD/gPb/';
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
    img_name = imgList(iI).name(1:end-4);
    
    [~, edgeim, thetamap] = load_edg([gPb_path img_name '.edg']);

%     load([gPb_path img_name '_gPb.mat'])
%     edgeim = gPb_thin;
%     norient=8;
%     gtheta = [1.5708    1.1781    0.7854    0.3927   0    2.7489    2.3562    1.9635];
% 
%     [unused, maxo] = max(gPb_orient,[],3);
%     thetamap = gtheta(maxo); %orientation of maxima

    % Link edge pixels together into lists of sequential edge points, one
    % list for each edge contour. A contour/edgelist starts/stops at an 
    % ending or a junction with another contour/edgelist.
    % Here we discard contours less than 18 pixels long.

    [edgelist, labelededgeim] = edgelink((edgeim>0.07), 5);
    edgelist2 = edgelist;
    
    for c = 1:length(edgelist)
        cur_c = zeros(size(edgelist{c},1), 4);
        cur_c(:,1) = edgelist{c}(:,2)-1;
        cur_c(:,2) = edgelist{c}(:,1)-1;
        for e = 1:size(cur_c,1)
            cur_c(e,3) = thetamap(cur_c(e,2)+1,cur_c(e,1)+1);
            cur_c(e,4) = edgeim(cur_c(e,2)+1,cur_c(e,1)+1);
        end
        edgelist2{c} = cur_c;
    end
    
    % here edgelist is matlab coordinates, need to change to cxx before
    % write_cem
    cem_file = [outPath '/' imgList(iI).name(1:end-4)  '.cem'];
    write_cem(cem_file, edgelist2, size(im,1), size(im,2));
    smooth_cem(cem_file, cem_file);
    % Display the edgelists with random colours for each distinct edge 
    % in figure 2
    imshow(im, 'border', 'tight');
    hold on;
     draw_contours(edgelist2, 0, 1, [0 0 0], -1)
    hold off;
    filename = [outPath '/' imgList(iI).name(1:end-4)  '.png'];
%     saveas(gcf,filename);
    print(filename, '-dpng');
%     print_pdf(filename);
    toc
    close all
end
