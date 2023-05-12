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
p_min = 0.05; % smaller than this is meaningless to link
vary_p = 0.95:-0.05:0.05;

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

    cfrags_maps.edges = [];
    cfrags_maps.cfrags = [];
    cfrags_maps.cfrags_idx = [];
    
    for p_i = 1:length(vary_p)
        p_th = vary_p(p_i)
        
        if(p_th< p_min) % just duplicate the last sets
            cfrags_maps.edges = cat( 1, cfrags_maps.edges, {edges});
            cfrags_maps.cfrags = cat(1, cfrags_maps.cfrags, {cfrags});
            cfrags_maps.cfrags_idx = cat(1, cfrags_maps.cfrags_idx, {cfrags_idx});
            continue;
        end
        
        E = ucms(:,:,1)>p_th;
        % Link edge pixels together into lists of sequential edge points, one
        % list for each edge contour. A contour/edgelist starts/stops at an 
        % ending or a junction with another contour/edgelist.
        % Here we discard contours less than 18 pixels long.

        [edgelist, labelededgeim] = edgelink((E), 5);
        cfrags = edgelist;

        for c = 1:length(edgelist)
            cur_c = zeros(size(edgelist{c},1), 4);
            cur_c(:,1) = edgelist{c}(:,2);
            cur_c(:,2) = edgelist{c}(:,1);
            for e = 1:size(cur_c,1)
                cur_c(e,3) = O(cur_c(e,2),cur_c(e,1));
                cur_c(e,4) = E(cur_c(e,2),cur_c(e,1));
            end
            cfrags{c} = cur_c;
        end

        % here edgelist is matlab coordinates, need to change to cxx before
        %%%%%%%%%%%%%% write results
        [yy,xx,mag] = find(E~=0);
        Angle = zeros(length(xx),1); 
        for k=1:length(xx)
            Angle(k) = O(yy(k),xx(k));
        end
        edges  = [xx yy Angle mag];

        % construct edge idx map
        idx_map = zeros(size(E));
        for eid = 1:size(edges, 1)
            idx_map(edges(eid, 2), edges(eid, 1)) = eid;
        end
        % construct cfrags_idx, contours recording edge idx
        cfrags_idx = cell(size(cfrags));
        for c = 1:length(cfrags)
            cfrags_idx{c} = zeros(1, size(cfrags{c}, 1));
            for ii=1:size(cfrags{c}, 1)
                cfrags_idx{c}(ii) = idx_map(cfrags{c}(ii,2), cfrags{c}(ii,1));
            end

            % change coordinates to original size
            cfrags{c}(:,1)= (cfrags{c}(:,1)/2)-1;
            cfrags{c}(:,2)= (cfrags{c}(:,2)/2)-1;

        end

        % change coordinates to original size
        edges(:,1) = (edges(:,1)/2)-1;
        edges(:,2) = (edges(:,2)/2)-1;
    %        
    %     % save cem with idx consistent with edge map
    %     cem_file = [outPath image_names{iI}  'im0.cem'];
    %     write_cem_fixed_edge_id(cem_file, cfrags, edges,  size(im,1), size(im,2), cfrags_idx)
        cfrags_maps.edges = cat( 1, cfrags_maps.edges, {edges});
        cfrags_maps.cfrags = cat(1, cfrags_maps.cfrags, {cfrags});
        cfrags_maps.cfrags_idx = cat(1, cfrags_maps.cfrags_idx, {cfrags_idx});
    
    end
    
    out_file = [outPath img_name  '_cfrag_maps.mat'];
    save(out_file, 'cfrags_maps', 'vary_p');
end
