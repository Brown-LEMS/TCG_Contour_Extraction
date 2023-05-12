clear all; close all;
% addpath(genpath('/media/New_Volume/Research/Project_contour/third_order/pb/'));
addpath('/media/guoy/Research/Project_contour/MSEL_contour_extraction_full/util/io');
imgPath  = '/media/guoy/Research/Project_contour/MSEL_contour_extraction_full/Data/CFGD/img/';
src_path = '/media/guoy/Research/Project_contour/MSEL_contour_extraction_full/Data/CFGD/CFGD_SE/edges/'
% edgePath = '../../data/CFGD/Pb/';
outPath = '/media/guoy/Research/Project_contour/MSEL_contour_extraction_full/Data/CFGD/CFGD_SE/SE_ucm_CFGD/';
if ~exist(outPath,'dir')
    mkdir(outPath);
end
imgList= dir([imgPath '*.jpg']);
numImage = length(imgList);
th = 0.05;

for iI= 1:numImage
    tic
    im = imread([imgPath '/' imgList(iI).name]);
    [h,w,~] = size(im);
%     if (size(im,3)==3)
%         im= rgb2gray(im);
%     end
    img_name = imgList(iI).name(1:end-4)
    load([src_path img_name '.mat']);
    
%     E =  imresize(ucms(:,:,1), [h,w])>th;
    E = ucms(:,:,1);
    O = wrapToPi(O.angle+pi/2);
    O = imresize(O, size(E));
    % Link edge pixels together into lists of sequential edge points, one
    % list for each edge contour. A contour/edgelist starts/stops at an 
    % ending or a junction with another contour/edgelist.
    % Here we discard contours less than 18 pixels long.

    [edgelist, labelededgeim] = edgelink((E>th), 5);
    edgelist2 = edgelist;
    
    for c = 1:length(edgelist)
        cur_c = zeros(size(edgelist{c},1), 4);
        cur_c(:,1) = edgelist{c}(:,2);
        cur_c(:,2) = edgelist{c}(:,1);
        % use the contours to compute edge orientation
        dx = diff(cur_c(:,1));
        dy = diff(cur_c(:,2));
        theta = atan2(dy, dx);
        
        for e = 1:size(cur_c,1)
            if(e~=1 && e~=size(cur_c,1))
                cur_c(e,3) =  theta(e-1);
                O(cur_c(e,2),cur_c(e,1)) = theta(e-1);
            else
                cur_c(e,3) = O(cur_c(e,2),cur_c(e,1));
            end
            cur_c(e,4) = E(cur_c(e,2),cur_c(e,1));
        end
        % use the contours to compute edge orientation
        
        edgelist2{c} = cur_c;
    end
    
    % here edgelist is matlab coordinates, need to change to cxx before
    %% write_cem
    %%%%%% save GT edges  
    [yy,xx] = find(E~=0);
    mag = zeros(size(xx));
    Angle = zeros(length(xx),1); 
    for k=1:length(xx)
        Angle(k) = O(yy(k),xx(k));
        mag(k) = E(yy(k),xx(k));
    end
    edginfo  = [xx yy Angle mag];

    % construct edge idx map
    idx_map = zeros(size(E));
%     idx_map(sub2ind([h,w], edginfo(:,2), edginfo(:,1))) = 1:size(edginfo, 1);
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
        
        % change coordinates to original size
        edgelist2{c}(:,1)= floor(edgelist2{c}(:,1)/2);
        edgelist2{c}(:,2)= floor(edgelist2{c}(:,2)/2);

    end

    % change coordinates to original size
    edginfo(:,1) = floor(edginfo(:,1)/2);
    edginfo(:,2) = floor(edginfo(:,2)/2);
       
    % save cem with idx consistent with edge map
    cem_file = [outPath '/' imgList(iI).name(1:end-4)  '.cem'];
    write_cem_fixed_edge_id(cem_file, edgelist2, edginfo,  size(im,1), size(im,2), cfrags_idx)

    
%     %% Display the edgelists with random colours for each distinct edge 
%     [CEM, edges,cfrags_idx] = load_contours(cem_file);
%     H = figure(1);
%     imshow(im, 'border', 'tight');
%     hold on;
%     draw_contours(edgelist2);
% %     draw_contours(edgelist2, 0, 1, [0 0 0], -1)
%     G_test = construct_fac_graph_from_curve_fragments (cfrags_idx, edgelist2);
%     for v =1:length(G_test.var)
%         if(length(G_test.var(v).nbrs_fac)>=3)
%             c_x = round(edginfo(G_test.var(v).actual_edge_id, 1)+1);
%             c_y = round(edginfo(G_test.var(v).actual_edge_id, 2)+1);
%             plot(c_x, c_y, 'gx');
%         end
%     end
%     
%     hold off;
%     out_file = [outPath '/' imgList(iI).name(1:end-4)  '.pdf'];
%     set(H, 'PaperPositionMode','auto')
%     print(H,'-dpdf','-r0',out_file);
%     cmd = ['!pdfcrop ' out_file ' ' out_file];
%     eval(cmd);
end
