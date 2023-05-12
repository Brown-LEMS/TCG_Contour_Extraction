clear all; close all;
% addpath /media/guoy/Research/Print_PDF/
addpath(genpath('/media/guoy/Research/Project_contour/MSEL_contour_extraction_full/util/'))
imgPath  = '/media/guoy/Research/Project_contour/MSEL_contour_extraction_full/Data/CFGD_img/';
edgePath = '/media/guoy/Research/Project_contour/MSEL_contour_extraction_full/Data/TO_SEL_CFGD/edges/';
outPath = '/media/guoy/Research/Project_contour/MSEL_contour_extraction_full/Data/TO_Kovesi_CFGD/';
if ~exist(outPath,'dir')
    mkdir(outPath);
end
imgList= dir([imgPath '*.jpg']);
numImage = length(imgList);
for iI= 1:numImage
    iI
    tic
    im = imread([imgPath '/' imgList(iI).name]);
%     if (size(im,3)==3)
%         im= rgb2gray(im);
%     end
    img_name = imgList(iI).name(1:end-4);
    
    % Find edges using the Canny operator with default
    % edgeim = edge(im,'canny',[],5);
% 	edgeim = edge(im,'canny');
    
    edge_file = [edgePath img_name '.edg'];
    [edg, edgeim, thetamap] = load_edg(edge_file);

    figure(1), imshow(1-edgeim);
    %set(gcf,'PaperPositionMode','auto');
    %filename = [outPath '/' imgList(iI).name(1:end-4)  '_canny.png'];
    %saveas(gcf,filename);
    
    % Link edge pixels together into lists of sequential edge points, one
    % list for each edge contour. A contour/edgelist starts/stops at an 
    % ending or a junction with another contour/edgelist.
    % Here we discard contours less than 18 pixels long.
    fprintf('construct a complex matrix holding subpixel positions as complex\n');    
    location = zeros(size(edgeim));
    for j = 1:length(edg)
        x = edg(j,1)+1;
        y = edg(j,2)+1;
        row = int32(y);
        col = int32(x);
        location(row,col) = complex(y, x);
    end
    
    
%     [edgelist, labelededgeim] = edgelink(edgeim, 5);
    [edgelist, labelededgeim] = edgelink(edgeim, 5, location);
    edgelist2 = edgelist;
    
    for c = 1:length(edgelist)
        cur_c = zeros(size(edgelist{c},1), 4);
        cur_c(:,1) = edgelist{c}(:,2)-1;
        cur_c(:,2) = edgelist{c}(:,1)-1;
        for e = 1:size(cur_c,1)
            cur_c(e,3) = thetamap(round(cur_c(e,2))+1,round(cur_c(e,1))+1);
            cur_c(e,4) = edgeim(round(cur_c(e,2))+1,round(cur_c(e,1))+1);
        end
        edgelist2{c} = cur_c;
    end
    
        % here edgelist is matlab coordinates, need to change to cxx before
    % write_cem
    cem_file = [outPath '/' imgList(iI).name(1:end-4)  '.cem'];
    write_cem(cem_file, edgelist2, size(im,1), size(im,2));
%     smooth_cem(cem_file, cem_file);
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
    pause(.1)
end
