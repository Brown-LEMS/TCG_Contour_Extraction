clear all; close all;
addpath(genpath('/media/New_Volume/Research/Project_contour/third_order/pb/'));
addpath ~/Desktop/Print_PDF/
imgPath  = '../../data/video_capital/';
% edgePath = '../../data/CFGD/Pb/';
outPath = '../../data/video_capital/Pb/';
if ~exist(outPath,'dir')
    mkdir(outPath);
end
imgList= dir([imgPath '*.png']);
numImage = length(imgList);
for iI= 1:numImage
    tic
    im = imread([imgPath '/' imgList(iI).name]);
%     if (size(im,3)==3)
%         im= rgb2gray(im);
%     end
    img_name = imgList(iI).name(1:end-4);
    
    % Find edges using the Canny operator with default
    % edgeim = edge(im,'canny',[],5);
% 	edgeim = edge(im,'canny');
    
%     edge_file = [edgePath img_name '.edg'];
    [edgeim,thetamap] = pbBGTG(im2double(im));
%     [edg, edgeim, thetamap] = load_edg(edge_file);
    
    
%     figure(1), imshow(1-(edgeim>0.1));
    %set(gcf,'PaperPositionMode','auto');
    %filename = [outPath '/' imgList(iI).name(1:end-4)  '_canny.png'];
    %saveas(gcf,filename);
    
    % Link edge pixels together into lists of sequential edge points, one
    % list for each edge contour. A contour/edgelist starts/stops at an 
    % ending or a junction with another contour/edgelist.
    % Here we discard contours less than 18 pixels long.

    [edgelist, labelededgeim] = edgelink((edgeim>0.1), 5);
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
    
    cem_file = [outPath '/' imgList(iI).name(1:end-4)  '.cem'];
    write_cem(cem_file, edgelist2, size(im,1), size(im,2));
    smooth_cem(cem_file, cem_file);
    % Display the edgelists with random colours for each distinct edge 
    % in figure 2
%     imshow(im, 'border', 'tight');
%     hold on;
%     h=drawedgelist(edgelist, size(im), 2, 'rand', 5); axis off 
%     hold off;
%     filename = [outPath '/' imgList(iI).name(1:end-4)  '.pdf'];
%     saveas(gcf,filename);
%     print_pdf(filename);
    toc
%     close all
end
