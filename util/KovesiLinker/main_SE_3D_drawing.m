clear all; close all;
% addpath(genpath('/media/New_Volume/Research/Project_contour/third_order/pb/'));
addpath(genpath('../'))
imgPath  = '/media/guoy/Research/Datasets/3D_Drawing/pavilion-midday-mcs-work/';
SE_path = '/media/guoy/Research/Datasets/3D_Drawing/pavilion-midday-mcs-work/SE_edges/';
% edgePath = '../../data/CFGD/Pb/';
outPath = '/media/guoy/Research/Datasets/3D_Drawing/pavilion-midday-mcs-work/SE_Kovesi';
if ~exist(outPath,'dir')
    mkdir(outPath);
end
prune_len = 5;

imgList= dir([imgPath '*.png']);
numImage = length(imgList);
for iI= 1:numImage
    tic
    im = imread([imgPath '/' imgList(iI).name]);
%     if (size(im,3)==3)
%         im= rgb2gray(im);
%     end
    img_name = imgList(iI).name(1:end-4)
    
    [~, edgeim, thetamap] = load_edg([SE_path img_name '.edg']);


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
    P_vec = [];
    for i = 1:length( edgelist2)
        p = contour_length_mex(edgelist2{i}');
        P_vec = [P_vec p];
    end 
    top_cfrags =  edgelist2(P_vec>prune_len);
    % here edgelist is matlab coordinates, need to change to cxx before
    % write_cem
    cem_file = [outPath '/' imgList(iI).name(1:end-4)  '.cemv'];
%     write_cem(cem_file, edgelist2, size(im,1), size(im,2));
    det_save_cemv(cem_file, top_cfrags);

    % Display the edgelists with random colours for each distinct edge 
    % in figure 2
    
    H = figure(1);
    imshow(im, 'border', 'tight');
    hold on;
     draw_contours(top_cfrags)
    hold off;
    set(H, 'PaperPositionMode','auto')
    filename = [outPath '/' imgList(iI).name(1:end-4)  '.png'];
%     saveas(gcf,filename);
    print(H,'-dpng','-r0', filename);
%     print_pdf(filename);
    toc
    close all
end
