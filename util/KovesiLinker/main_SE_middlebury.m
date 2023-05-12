clear all; close all;
% addpath(genpath('/media/New_Volume/Research/Project_contour/third_order/pb/'));
addpath('/media/guoy/Research/Project_contour/MSEL_contour_extraction_full/util/io');
imgPath  = '/media/guoy/Research/Project_contour/MSEL_contour_extraction_full/Data/Middlebury/trainingQ/';
src_path = '/media/guoy/Research/Project_contour/MSEL_contour_extraction_full/Data/Middlebury/trainingQ_SE/';
outPath = '/media/guoy/Research/Project_contour/MSEL_contour_extraction_full/Data/Middlebury/trainingQ_SE/results_Kovesi/';
if ~exist(outPath,'dir')
    mkdir(outPath);
end
imgList= dir([imgPath '*.png']);
th = 0.1;

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
    iI
    %% %%%%%%%%%%%%%%%%% im0 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    im = imread([imgPath image_names{iI} '/im0.png']);
    [h,w,~] = size(im);
    [~, E, O] = load_edg([src_path image_names{iI} '/im0.edg']);

    % Link edge pixels together into lists of sequential edge points, one
    % list for each edge contour. A contour/edgelist starts/stops at an 
    % ending or a junction with another contour/edgelist.
    % Here we discard contours less than 18 pixels long.

    [edgelist, labelededgeim] = edgelink((E>0), 5);
    edgelist2 = edgelist;
    
    for c = 1:length(edgelist)
        cur_c = zeros(size(edgelist{c},1), 4);
        cur_c(:,1) = edgelist{c}(:,2);
        cur_c(:,2) = edgelist{c}(:,1);
        for e = 1:size(cur_c,1)
            cur_c(e,3) = O(cur_c(e,2),cur_c(e,1));
            cur_c(e,4) = E(cur_c(e,2),cur_c(e,1));
        end
        edgelist2{c} = cur_c;
    end
    
    % here edgelist is matlab coordinates, need to change to cxx before
    %% write_cem
    %%%%%% save GT edges  
    [yy,xx,mag] = find(E~=0);
    Angle = zeros(length(xx),1); 
    for k=1:length(xx)
        Angle(k) = O(yy(k),xx(k));
    end
    edginfo  = [xx yy Angle mag];

    % construct edge idx map
    idx_map = zeros(size(E));
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
        edgelist2{c}(:,1:2)=  edgelist2{c}(:,1:2)-1;

    end

    % change coordinates to original size
    edginfo(:,1:2) = edginfo(:,1:2)-1;

       
    % save cem with idx consistent with edge map
    cem_file = [outPath image_names{iI}  'im0.cem'];
    write_cem_fixed_edge_id(cem_file, edgelist2, edginfo,  size(im,1), size(im,2), cfrags_idx)

    
    %% Display the edgelists with random colours for each distinct edge 
    [CEM, edges,cfrags_idx] = load_contours(cem_file);
    H = figure(1);
    imshow(im, 'border', 'tight');
    hold on;
    draw_contours(edgelist2);
%     draw_contours(edgelist2, 0, 1, [0 0 0], -1)
%     G_test = construct_fac_graph_from_curve_fragments (cfrags_idx, edgelist2);
%     for v =1:length(G_test.var)
%         if(length(G_test.var(v).nbrs_fac)>=3)
%             c_x = (edginfo(G_test.var(v).actual_edge_id, 1)+1);
%             c_y = (edginfo(G_test.var(v).actual_edge_id, 2)+1);
%             plot(c_x, c_y, 'gx');
%         end
%     end
    
    hold off;
    out_file = [outPath image_names{iI}  'im0.pdf'];
    set(H, 'PaperPositionMode','auto')
    print(H,'-dpdf','-r0',out_file);
    cmd = ['!pdfcrop ' out_file ' ' out_file];
    eval(cmd);
    
    
    
     %% %%%%%%%%%%%%%%%%% im1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    im = imread([imgPath image_names{iI} '/im1.png']);
    [h,w,~] = size(im);
    [~, E, O] = load_edg([src_path image_names{iI} '/im1.edg']);

    % Link edge pixels together into lists of sequential edge points, one
    % list for each edge contour. A contour/edgelist starts/stops at an 
    % ending or a junction with another contour/edgelist.
    % Here we discard contours less than 18 pixels long.

    [edgelist, labelededgeim] = edgelink((E>0), 5);
    edgelist2 = edgelist;
    
    for c = 1:length(edgelist)
        cur_c = zeros(size(edgelist{c},1), 4);
        cur_c(:,1) = edgelist{c}(:,2);
        cur_c(:,2) = edgelist{c}(:,1);
        for e = 1:size(cur_c,1)
            cur_c(e,3) = O(cur_c(e,2),cur_c(e,1));
            cur_c(e,4) = E(cur_c(e,2),cur_c(e,1));
        end
        edgelist2{c} = cur_c;
    end
    
    % here edgelist is matlab coordinates, need to change to cxx before
    %% write_cem
    %%%%%% save GT edges  
    [yy,xx,mag] = find(E~=0);
    Angle = zeros(length(xx),1); 
    for k=1:length(xx)
        Angle(k) = O(yy(k),xx(k));
    end
    edginfo  = [xx yy Angle mag];

    % construct edge idx map
    idx_map = zeros(size(E));
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
        edgelist2{c}(:,1:2)=  edgelist2{c}(:,1:2)-1;

    end

    % change coordinates to original size
    edginfo(:,1:2) = edginfo(:,1:2)-1;
       
    % save cem with idx consistent with edge map
    cem_file = [outPath image_names{iI}  'im1.cem'];
    write_cem_fixed_edge_id(cem_file, edgelist2, edginfo,  size(im,1), size(im,2), cfrags_idx)

    
    %% Display the edgelists with random colours for each distinct edge 
    [CEM, edges,cfrags_idx] = load_contours(cem_file);
    H = figure(1);
    imshow(im, 'border', 'tight');
    hold on;
    draw_contours(edgelist2);
%     draw_contours(edgelist2, 0, 1, [0 0 0], -1)
%     G_test = construct_fac_graph_from_curve_fragments (cfrags_idx, edgelist2);
%     for v =1:length(G_test.var)
%         if(length(G_test.var(v).nbrs_fac)>=3)
%             c_x = (edginfo(G_test.var(v).actual_edge_id, 1)+1);
%             c_y = (edginfo(G_test.var(v).actual_edge_id, 2)+1);
%             plot(c_x, c_y, 'gx');
%         end
%     end
    
    hold off;
    out_file = [outPath image_names{iI}  'im1.pdf'];
    set(H, 'PaperPositionMode','auto')
    print(H,'-dpdf','-r0',out_file);
    cmd = ['!pdfcrop ' out_file ' ' out_file];
    eval(cmd);
end
