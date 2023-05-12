clear all; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% the specialty of UCM: at each prob threshold, it make different
%%% close regions, the edge groupings are different
%%% A: cfrags can be extracted at lowest threshold, the high threshold keep
%%% pruning out curves, and grouping remaining ones: might erronous due to
%%% linking at lowest level
%%% B: just run linking seperately at each threshold (current)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p_min = 0.1; % smaller than this is meaningless to link
vary_p = 0.85:-0.05:0;

% addpath(genpath('/media/New_Volume/Research/Project_contour/third_order/pb/'));
addpath('/media/guoy/Research/Project_contour/MSEL_contour_extraction_full/util/io');
imgPath  = '/media/guoy/Research/Project_contour/MSEL_contour_extraction_full/Data/Middlebury/trainingQ/';
src_path = '/media/guoy/Research/Project_contour/MSEL_contour_extraction_full/Data/Middlebury/trainingQ_SE/';
% edgePath = '../../data/CFGD/Pb/';
outPath = '/media/guoy/Research/Project_contour/MSEL_contour_extraction_full/Data/Middlebury/trainingQ_SE/results_ucm/';
if ~exist(outPath,'dir')
    mkdir(outPath);
end

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

for iI= 1:length(image_names);
    tic
    iI
    %% %%%%%%%%%%%%%%%%% im0 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    im = imread([imgPath image_names{iI} '/im0.png']);
    [h,w,~] = size(im);
    load([src_path image_names{iI} '/im0.mat']);
    
    
    %     E =  imresize(ucms(:,:,1), [h,w])>th;
    E = ucms(:,:,1);
    O = wrapToPi(O.angle+pi/2);
    O = imresize(O, size(E));
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
    
    out_file = [outPath image_names{iI}  'im0_cfrag_maps.mat'];
    save(out_file, 'cfrags_maps', 'vary_p');
    
    %% %%%%%%%%%%%%%%%%% im1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    im = imread([imgPath image_names{iI} '/im1.png']);
    [h,w,~] = size(im);
    load([src_path image_names{iI} '/im1.mat']);
    
    
    %     E =  imresize(ucms(:,:,1), [h,w])>th;
    E = ucms(:,:,1);
    O = wrapToPi(O.angle+pi/2);
    O = imresize(O, size(E));
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
           
    %     % save cem with idx consistent with edge map
    %     cem_file = [outPath image_names{iI}  'im0.cem'];
    %     write_cem_fixed_edge_id(cem_file, cfrags, edges,  size(im,1), size(im,2), cfrags_idx)
        cfrags_maps.edges = cat( 1, cfrags_maps.edges, {edges});
        cfrags_maps.cfrags = cat(1, cfrags_maps.cfrags, {cfrags});
        cfrags_maps.cfrags_idx = cat(1, cfrags_maps.cfrags_idx, {cfrags_idx});
    
    end
    
    out_file = [outPath image_names{iI}  'im1_cfrag_maps.mat'];
    save(out_file, 'cfrags_maps', 'vary_p');
    
end
