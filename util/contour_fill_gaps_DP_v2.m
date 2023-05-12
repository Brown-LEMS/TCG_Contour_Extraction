function [cfrags, cfrags_idx, edges] = contour_fill_gaps_DP_v2(cfrags, cfrags_idx, h, w, edges, edgemap_soft, O_map, params)
% Last Modified: 02/27/2018
% DP-augmented curves extend existing curves.
% max_iter =2;
r_range = params.DP_gap_range;
theta_range = params.DP_angle_th;
search_nbr_range = params.shape_gap_range; % search range from an end of curve fragment 
search_ori_range = params.shape_ori_range; % search angle at an end point of cfrag
cost_th = 2;
nbr_xx_coordinates = repmat(-search_nbr_range:search_nbr_range, [2*search_nbr_range+1, 1]);
nbr_yy_coordinates = nbr_xx_coordinates';
nbr_rr = sqrt(nbr_xx_coordinates.^2 + nbr_yy_coordinates.^2);
nbr_theta = atan2(nbr_yy_coordinates, nbr_xx_coordinates);

%% preparation of reference maps
G_test = construct_fac_graph_from_curve_fragments (cfrags_idx, cfrags);
end_point_eid_map = zeros(h,w);
end_point_vid_map = zeros(h,w);
end_eids = [G_test.var(:).actual_edge_id];
end_Y = round(edges(end_eids', 2))+1;
end_X = round(edges(end_eids', 1))+1;
end_point_eid_map(sub2ind([h,w], end_Y, end_X))=end_eids;
end_point_vid_map(sub2ind([h,w], end_Y, end_X))=[G_test.var(:).id];

%%%%%%%%% construct contour id map
% ATTENTION: update EdgeGroupMap requires to use interpolated edges (inefficient)
EdgeGroupMap = convert_cfrags_to_EdgeGroupMap (cfrags, h, w);

%%%%%%%%%% construct junction map (This can be removed is T-junction is known)
EDGEIM = EdgeGroupMap ~= 0;           % make sure image is binary.
EDGEIM = bwmorph(EDGEIM,'clean');     % Remove isolated pixels
EDGEIM = bwmorph(EDGEIM,'skel',Inf);  % and make sure edges are thinned. I
                                      % think using 'skel' is better than 'thin'    
% Find endings and junctions in edge data
jct_map = zeros(h,w);
[Jct_Y, Jct_X, re, ce] = findendsjunctions(EDGEIM);
jct_map(sub2ind([h,w], Jct_Y, Jct_X))=1;
jct_map = conv2(jct_map,ones(3,3),'same');


%% %%%%%%%%%%% iterate through free ends and look for direction completion
% G_test = construct_fac_graph_from_curve_fragments (cfrags_idx, cfrags);
% prev_last_eid = size(edges, 1);
for vid =1:length(G_test.var)
   
    % extact the contour end
    % TODO: consider corner contours
    if(length(G_test.var(vid).nbrs_fac)~=1)
        edges(G_test.var(vid).actual_edge_id, 5) = length(G_test.var(vid).nbrs_fac);
        continue;
    end
    edges(G_test.var(vid).actual_edge_id, 5) = max(edges(G_test.var(vid).actual_edge_id, 5), 1);
    
    % skip marked T-junction
    if(edges(G_test.var(vid).actual_edge_id, 5)>1)
        continue;
    end
    
    cur_x = edges(G_test.var(vid).actual_edge_id,1);
    cur_y = edges(G_test.var(vid).actual_edge_id,2);
    cur_x = round(cur_x)+1;
    cur_y = round(cur_y)+1;
    
    if(jct_map(cur_y, cur_x))
        edges(G_test.var(vid).actual_edge_id, 5) = 3;
        continue;
    end
    % compute the orientation at free end
    cid = G_test.var(vid).nbrs_fac;
    cur_c = cfrags{cid};
    c_len =  contour_length_mex(cur_c');
    e_size = size(cur_c,1);
    
    % skip short contours
    if(e_size<10)
        continue;
    end
    
    is_start = 1;
    if(G_test.fac(cid).nbrs_var(1)==vid)
        end_dx =  cur_c(1,1) - cur_c(min(e_size,5),1);
        end_dy =  cur_c(1,2) - cur_c(min(e_size,5),2);
    else
        is_start = 0;
        end_dx =  cur_c(end,1) - cur_c(max(e_size-5+1, 1),1);
        end_dy =  cur_c(end,2) - cur_c(max(e_size-5+1, 1),2);
    end
        
 %%%%%%%% case 2.1 search a fan area to reach a contour end point
    add_end_edge_idx1=0;
    angle = atan2(end_dy, end_dx);
    angle_max = angle+search_ori_range;
    angle_min = angle-search_ori_range;
    if(angle_max>pi)
        qualify_idx = ((nbr_theta(:) >= angle_min & nbr_theta(:) <= pi) | ...
                        (nbr_theta(:) > -pi & nbr_theta(:) <= angle_max-2*pi)) ...
                        & nbr_rr(:)<=search_nbr_range;
    elseif(angle_min<-pi)
        qualify_idx = ((nbr_theta(:) >= angle_min +2*pi & nbr_theta(:) <= pi) | ...
                        (nbr_theta(:) > -pi & nbr_theta(:) <= angle_max))...
                        & nbr_rr(:)<=search_nbr_range;
    else
        qualify_idx = nbr_theta(:) >= angle_min & nbr_theta(:) <= angle_max & nbr_rr(:)<=search_nbr_range;
    end

    % allow larger angle in close range
    angle_max = angle+pi/3;
    angle_min = angle-pi/3;
    if(angle_max>pi)
        qualify_idx = qualify_idx | ((nbr_theta(:) >= angle_min & nbr_theta(:) <= pi) | ...
                        (nbr_theta(:) > -pi & nbr_theta(:) <= angle_max-2*pi)) ...
                        & nbr_rr(:)<=2;
    elseif(angle_min<-pi)
        qualify_idx = qualify_idx | ((nbr_theta(:) >= angle_min +2*pi & nbr_theta(:) <= pi) | ...
                        (nbr_theta(:) > -pi & nbr_theta(:) <= angle_max))...
                        & nbr_rr(:)<=2;
    else
        qualify_idx = qualify_idx | (nbr_theta(:) >= angle_min & nbr_theta(:) <= angle_max & nbr_rr(:)<=2);
    end

    x_coords_end = nbr_xx_coordinates(qualify_idx)+cur_x;
    y_coords_end = nbr_yy_coordinates(qualify_idx)+cur_y;
    x_coords_end = max(x_coords_end, ones(size(x_coords_end)));
    y_coords_end = max(y_coords_end, ones(size(y_coords_end)));
    x_coords_end = min(x_coords_end, w*ones(size(x_coords_end)));
    y_coords_end = min(y_coords_end, h*ones(size(y_coords_end)));

    reach_vids = end_point_vid_map(sub2ind([h,w], y_coords_end, x_coords_end));
    reach_vids(reach_vids<=0) =[];
    reach_vids(reach_vids==vid) =[];
    if(~isempty(reach_vids))
        add_end_edge_idx1 = G_test.var(reach_vids(1)).actual_edge_id;
        dist_1 = sqrt((cur_x-edges(add_end_edge_idx1,1)-1)^2 + (cur_y-edges(add_end_edge_idx1,2)-1)^2);
        if(dist_1 > c_len)
            add_end_edge_idx1 = 0;
        end
    end


    %%%%%%%% case 2.2 search on the extending line for a T_junction
    add_end_edge_idx2 = 0;
    % only consider half range for T-junction
    line_end_x = cur_x + cos(angle)*search_nbr_range*0.5;
    line_end_y = cur_y + sin(angle)*search_nbr_range*0.5;
    line_X = linspace(cur_x, line_end_x, search_nbr_range);
    line_Y = linspace(cur_y, line_end_y, search_nbr_range);
    line_X = max(line_X, 1);
    line_X = min(line_X, w);
    line_Y = max(line_Y, 1);
    line_Y = min(line_Y, h);            

    % mark T-junction
    reach_cids = EdgeGroupMap(sub2ind([h,w], round(line_Y), round(line_X)));
    reach_cids(reach_cids==cid)=0;
    inter_idx = find(reach_cids);
    if(~isempty(inter_idx))
        inter_idx = inter_idx(1);
        inter_x = line_X(inter_idx);
        inter_y = line_Y(inter_idx);
        reach_cid = reach_cids(inter_idx);

        attach_cfrag = cfrags{reach_cid}(:, 1:2)+1;
        attach_cfrag_idx = cfrags_idx{reach_cid};
        c_len = length(attach_cfrag_idx);
        dist2cfrag = attach_cfrag - repmat([inter_x, inter_y], c_len, 1);
        dist2cfrag = sqrt(dist2cfrag(:,1).^2 + dist2cfrag(:,2).^2);
        [temp_min_c_dist, r_idx] = min(dist2cfrag);

        add_end_edge_idx2 = attach_cfrag_idx(r_idx); 
        dist_2 = sqrt((cur_x-edges(add_end_edge_idx2,1)-1)^2 + (cur_y-edges(add_end_edge_idx2,2)-1)^2);
        if(dist_2 > c_len)
            add_end_edge_idx2 = 0;
        end
    end
    
    % compare two choices and choose the shorter one
    c_type =0;
    if(add_end_edge_idx1 && ~add_end_edge_idx2)
        c_type = 1;
        add_end_edge_idx = add_end_edge_idx1;
    elseif(~add_end_edge_idx1 && add_end_edge_idx2)
        c_type = 2;
        add_end_edge_idx = add_end_edge_idx2;
    elseif(add_end_edge_idx1 && add_end_edge_idx2)
        if(dist_1 < dist_2)
            c_type = 1;
            add_end_edge_idx = add_end_edge_idx1;
        else
            c_type = 2;
            add_end_edge_idx = add_end_edge_idx2;
        end
    end
    
    
    if(c_type==1)
        if(is_start)
            cfrags{cid} = [edges(add_end_edge_idx,:); cfrags{cid}];
            cfrags_idx{cid} = [add_end_edge_idx cfrags_idx{cid}];
            G_test.fac(cid).nbrs_var(1)=reach_vids(1);
        else
            cfrags{cid} = [cfrags{cid}; edges(add_end_edge_idx,:)];
            cfrags_idx{cid} = [cfrags_idx{cid} add_end_edge_idx];            
            G_test.fac(cid).nbrs_var(2)=reach_vids(1);
        end
        G_test.var(reach_vids(1)).nbrs_fac = [G_test.var(reach_vids(1)).nbrs_fac; cid];
        G_test.var(vid).nbrs_fac(G_test.var(vid).nbrs_fac==cid) = [];
        % mark connected edge
        edges(add_end_edge_idx,5) = 2;
        end_point_eid_map(cur_y, cur_x) = 0;
        plot([cur_x, edges(add_end_edge_idx,1)+1], [cur_y, edges(add_end_edge_idx,2)+1], 'g-');
    elseif(c_type==2)
        if(is_start)
            cfrags{cid} = [edges(add_end_edge_idx,:); cfrags{cid}];
            cfrags_idx{cid} = [add_end_edge_idx cfrags_idx{cid}];
        else
            cfrags{cid} = [cfrags{cid}; edges(add_end_edge_idx,:)];
            cfrags_idx{cid} = [cfrags_idx{cid} add_end_edge_idx];            
        end

        % Mark T-junction
        edges(add_end_edge_idx,5) = 3;
        jct_map(round(edges(add_end_edge_idx,2)+1), round(edges(add_end_edge_idx,1))+1) = 1; 
        end_point_eid_map(cur_y, cur_x) = 0;
        plot([cur_x, edges(add_end_edge_idx,1)+1], [cur_y, edges(add_end_edge_idx,2)+1], 'g-');
    end
    
    
end


%% Dynamic programming at each end point and find the optimal path within a range
% DT of binary map
DT_map = bwdist(EdgeGroupMap>0,'euclidean');
DT_map = double(exp(-DT_map.^2/2^2/2));
% for iter = 1:max_iter
prev_last_eid = size(edges, 1);
for vid =1:length(G_test.var)

    %%%%%%%%% Prepare input maps to DP 
    
    % extact the contour end
    % TODO: consider corner contours
    if(length(G_test.var(vid).nbrs_fac)~=1)
        edges(G_test.var(vid).actual_edge_id, 5) = length(G_test.var(vid).nbrs_fac);
        continue;
    end
    
    % skip marked junctions
    if(edges(G_test.var(vid).actual_edge_id, 5)>1)
        continue;
    end
    
    cur_x = edges(G_test.var(vid).actual_edge_id,1);
    cur_y = edges(G_test.var(vid).actual_edge_id,2);
    cur_x = round(cur_x)+1;
    cur_y = round(cur_y)+1;
    % skip detected junctions
    if(jct_map(cur_y, cur_x))
        edges(G_test.var(vid).actual_edge_id, 5) = 3;
        continue;
    end
    % skip completed end previous
    if(end_point_eid_map(cur_y, cur_x)==0)
        continue;
    end
    % compute the orientation at free end
    cid = G_test.var(vid).nbrs_fac;
    cur_c = cfrags{cid};
    e_size = size(cur_c,1);
    
    is_start = 1;
    if(G_test.fac(cid).nbrs_var(1)==vid)
        end_dx =  cur_c(1,1) - cur_c(min(e_size,4),1);
        end_dy =  cur_c(1,2) - cur_c(min(e_size,4),2);
    else
        end_dx =  cur_c(end,1) - cur_c(max(e_size-3, 1),1);
        end_dy =  cur_c(end,2) - cur_c(max(e_size-3, 1),2);
        is_start = 0;
    end

    start_point = [cur_x, cur_y, atan2(end_dy, end_dx)];

    x_min = max(cur_x-r_range,1);
    x_max = min(cur_x+r_range,w);
    y_min = max(cur_y-r_range,1);
    y_max = min(cur_y+r_range,h);    
    
    start_point_ref = start_point; % this is in c++ coordinates
    start_point_ref(1) = start_point_ref(1) - x_min;
    start_point_ref(2) = start_point_ref(2) - y_min;
    h_ref = y_max-y_min+1;
    w_ref = x_max-x_min+1;
    
    EdgeGroupMap_ref = EdgeGroupMap(y_min:y_max, x_min:x_max);
    edgemap_ref = DT_map(y_min:y_max, x_min:x_max);
    edgemap_ref(EdgeGroupMap_ref>0) = 1;

    %%%%%%%%%% DP: make sure all the matrix input in double format
    [cost_mat, back_p_mat, len_mat] = DP_gap_cpt(edgemap_ref, edgemap_soft(y_min:y_max, x_min:x_max), ...
                                            O_map(y_min:y_max, x_min:x_max), start_point_ref, h_ref, w_ref, theta_range);
                                                                      
    cost_mat_norm = cost_mat./len_mat; % + ones(size(len_mat))./len_mat;
    
    %%%%%%%%%% A: do not extend longer than the connected contour
    cost_mat_norm(len_mat>r_range) = 1000;

    %%%%%%%%%% B: rule out paths ending on existing junctions
    jct_map_ref = jct_map(y_min:y_max, x_min:x_max);
    cost_mat_norm(jct_map_ref>0) = 1000;
    
    [opt_s, opt_ind] = min(cost_mat_norm(:));
    if(opt_s>cost_th);
        continue;
    end

    %%%%%%%%%% Backtrack to extract the optimal path
    [end_point_y, end_point_x] = ind2sub([h_ref, w_ref], opt_ind); 
    
    cur_ind = sub2ind([h_ref, w_ref], end_point_y,end_point_x);
    opt_path = cur_ind;
    while(back_p_mat(cur_ind))
        opt_path = [opt_path back_p_mat(cur_ind)+1];
        cur_ind = back_p_mat(cur_ind)+1; % covert to matlab coordinates
    end
    [opt_path_y, opt_path_x] = ind2sub([h_ref,w_ref], opt_path);
    
    opt_path_y = opt_path_y+y_min-1;
    opt_path_x = opt_path_x+x_min-1;
    
    
    %%%%%%%%%% C: rule out those gap paths that are longer than its connected contour and are very low avg contrast
    path_prob = edgemap_soft(sub2ind([h,w], opt_path_y(1:end-1), opt_path_x(1:end-1)));
    if((mean(path_prob)<0.05 && length(path_prob)> 5)) %|| (mean(path_prob)<0.05 && length(path_prob)<= 5))
        continue;
    end

    
    %%%%%%%%%% case 1 reach an existing contour 
    path_end_x = opt_path_x(1);
    path_end_y = opt_path_y(1);
    replace_end_edge_idx = 0;
    if(EdgeGroupMap(path_end_y, path_end_x) >0)
        reach_cid = EdgeGroupMap(path_end_y, path_end_x);
        attach_cfrag = cfrags{reach_cid}(:, 1:2)+1;
        attach_cfrag_idx = cfrags_idx{reach_cid};
        
        c_len = length(attach_cfrag_idx);
        dist2cfrag = attach_cfrag - repmat([path_end_x, path_end_y], c_len, 1);
        dist2cfrag = dist2cfrag(:,1).^2 + dist2cfrag(:,2).^2;
        [~, r_idx] = min(dist2cfrag);
        opt_path_x(1) = attach_cfrag(r_idx, 1);
        opt_path_y(1) = attach_cfrag(r_idx, 2);
        replace_end_edge_idx = attach_cfrag_idx(r_idx);
        
        if(r_idx~=1 && r_idx~=c_len) % mark T-junction
            edges(replace_end_edge_idx,5) = 3;
            jct_map(path_end_y, path_end_x) =1;
            
        else % mark connected edge
            edges(replace_end_edge_idx,5) = 2;
        end 
    end
    
    plot(opt_path_x, opt_path_y, 'r-');
    end_point_eid_map(cur_y, cur_x) = 0;
    x2update = round(opt_path_x(1));
    y2update = round(opt_path_y(1));
    
    %%%%%%%%%% extend current contour
    if(is_start) % extend the starting edge of a contour
        opt_path_dx = diff(opt_path_x);
        opt_path_dy = diff(opt_path_y);
        opt_path_theta = atan2(opt_path_dy, opt_path_dx);
        
        add_edges = [opt_path_x(1:end-1)'-1, opt_path_y(1:end-1)'-1, opt_path_theta', path_prob', zeros(size(path_prob'))];
        cfrags{cid} = [add_edges; cfrags{cid}];
        
        if(replace_end_edge_idx )
            add_edges(1,:) = [];
            cfrags{cid}(1,:) = edges(replace_end_edge_idx,:);
        end
        edges = [edges; add_edges];
        cur_last_eid = size(edges, 1);
        cfrags_idx{cid} = [prev_last_eid+1:cur_last_eid cfrags_idx{cid}];
        eid2update = prev_last_eid+1;
        if(replace_end_edge_idx)
            cfrags_idx{cid} = [replace_end_edge_idx cfrags_idx{cid}];
%             eid2update = replace_end_edge_idx;
            eid2update = 0;
        end
        
        prev_last_eid = cur_last_eid;
    else % extend the ending edge of a contour
        opt_path_x = fliplr(opt_path_x);
        opt_path_y = fliplr(opt_path_y);
        path_prob = fliplr(path_prob);
        
        opt_path_dx = diff(opt_path_x);
        opt_path_dy = diff(opt_path_y);
        opt_path_theta = atan2(opt_path_dy, opt_path_dx);
        
        add_edges = [opt_path_x(2:end)'-1, opt_path_y(2:end)'-1, opt_path_theta', path_prob', zeros(size(path_prob'))];
        cfrags{cid} = [cfrags{cid}; add_edges];      
        
        if(replace_end_edge_idx)
            add_edges(end,:) = [];
            cfrags{cid}(end,:) = edges(replace_end_edge_idx,:);
        end
        
        edges = [edges; add_edges];
        cur_last_eid = size(edges, 1);
        cfrags_idx{cid} = [cfrags_idx{cid} prev_last_eid+1:cur_last_eid];
        eid2update = cur_last_eid;
        if(replace_end_edge_idx)
            cfrags_idx{cid} = [cfrags_idx{cid} replace_end_edge_idx];
%             eid2update = replace_end_edge_idx;
            eid2update = 0;
        end
        
        prev_last_eid = cur_last_eid;
    end
        
    end_point_eid_map(y2update, x2update) = eid2update;
    EdgeGroupMap(sub2ind([h,w], round(opt_path_y), round(opt_path_x))) =  cid;

end

% end
%%%%%%%%%%%%% iterate through free ends and look for direction completion
G_test = construct_fac_graph_from_curve_fragments (cfrags_idx, cfrags);
end_point_eid_map = zeros(h,w);
end_point_vid_map = zeros(h,w);
end_eids = [G_test.var(:).actual_edge_id];
end_Y = round(edges(end_eids', 2))+1;
end_X = round(edges(end_eids', 1))+1;
end_point_eid_map(sub2ind([h,w], end_Y, end_X))=end_eids;
end_point_vid_map(sub2ind([h,w], end_Y, end_X))=[G_test.var(:).id];
for vid =1:length(G_test.var)
   
    % extact the contour end
    % TODO: consider corner contours
    if(length(G_test.var(vid).nbrs_fac)~=1)
        edges(G_test.var(vid).actual_edge_id, 5) = length(G_test.var(vid).nbrs_fac);
        continue;
    end
    edges(G_test.var(vid).actual_edge_id, 5) = max(edges(G_test.var(vid).actual_edge_id, 5), 1);
    
    % skip marked T-junction
    if(edges(G_test.var(vid).actual_edge_id, 5)>1)
        continue;
    end
    
    cur_x = edges(G_test.var(vid).actual_edge_id,1);
    cur_y = edges(G_test.var(vid).actual_edge_id,2);
    cur_x = round(cur_x)+1;
    cur_y = round(cur_y)+1;
        
    if(jct_map(cur_y, cur_x))
        edges(G_test.var(vid).actual_edge_id, 5) = 1;
        continue;
    end
    % compute the orientation at free end
    cid = G_test.var(vid).nbrs_fac;
    cur_c = cfrags{cid};
    c_len =  contour_length_mex(cur_c');
    e_size = size(cur_c,1);
    
    % skip short contours
    if(e_size<5)
        continue;
    end
    
    is_start = 1;
    if(G_test.fac(cid).nbrs_var(1)==vid)
        end_dx =  cur_c(1,1) - cur_c(min(e_size,5),1);
        end_dy =  cur_c(1,2) - cur_c(min(e_size,5),2);
    else
        is_start = 0;
        end_dx =  cur_c(end,1) - cur_c(max(e_size-5+1, 1),1);
        end_dy =  cur_c(end,2) - cur_c(max(e_size-5+1, 1),2);
    end
        
    %%%%%%%% case 2.1 search a fan area to reach a contour end point
    add_end_edge_idx1=0;
    angle = atan2(end_dy, end_dx);
    angle_max = angle+search_ori_range;
    angle_min = angle-search_ori_range;
    if(angle_max>pi)
        qualify_idx = ((nbr_theta(:) >= angle_min & nbr_theta(:) <= pi) | ...
                        (nbr_theta(:) > -pi & nbr_theta(:) <= angle_max-2*pi)) ...
                        & nbr_rr(:)<=search_nbr_range;
    elseif(angle_min<-pi)
        qualify_idx = ((nbr_theta(:) >= angle_min +2*pi & nbr_theta(:) <= pi) | ...
                        (nbr_theta(:) > -pi & nbr_theta(:) <= angle_max))...
                        & nbr_rr(:)<=search_nbr_range;
    else
        qualify_idx = nbr_theta(:) >= angle_min & nbr_theta(:) <= angle_max & nbr_rr(:)<=search_nbr_range;
    end

    % allow larger angle in close range
    angle_max = angle+pi/3;
    angle_min = angle-pi/3;
    if(angle_max>pi)
        qualify_idx = qualify_idx | ((nbr_theta(:) >= angle_min & nbr_theta(:) <= pi) | ...
                        (nbr_theta(:) > -pi & nbr_theta(:) <= angle_max-2*pi)) ...
                        & nbr_rr(:)<=2;
    elseif(angle_min<-pi)
        qualify_idx = qualify_idx | ((nbr_theta(:) >= angle_min +2*pi & nbr_theta(:) <= pi) | ...
                        (nbr_theta(:) > -pi & nbr_theta(:) <= angle_max))...
                        & nbr_rr(:)<=2;
    else
        qualify_idx = qualify_idx | (nbr_theta(:) >= angle_min & nbr_theta(:) <= angle_max & nbr_rr(:)<=2);
    end

    x_coords_end = nbr_xx_coordinates(qualify_idx)+cur_x;
    y_coords_end = nbr_yy_coordinates(qualify_idx)+cur_y;
    x_coords_end = max(x_coords_end, ones(size(x_coords_end)));
    y_coords_end = max(y_coords_end, ones(size(y_coords_end)));
    x_coords_end = min(x_coords_end, w*ones(size(x_coords_end)));
    y_coords_end = min(y_coords_end, h*ones(size(y_coords_end)));

    reach_vids = end_point_vid_map(sub2ind([h,w], y_coords_end, x_coords_end));
    reach_vids(reach_vids<=0) =[];
    reach_vids(reach_vids==vid) =[];
    if(~isempty(reach_vids))
        add_end_edge_idx1 = G_test.var(reach_vids(1)).actual_edge_id;
        dist_1 = sqrt((cur_x-edges(add_end_edge_idx1,1)-1)^2 + (cur_y-edges(add_end_edge_idx1,2)-1)^2);
        if(dist_1 > 2*c_len)
            add_end_edge_idx1 = 0;
        end
    end


    %%%%%%%% case 2.2 search on the extending line for a T_junction
    add_end_edge_idx2 = 0;
    % only consider half range for T-junction
    line_end_x = cur_x + cos(angle)*search_nbr_range*0.5;
    line_end_y = cur_y + sin(angle)*search_nbr_range*0.5;
    line_X = linspace(cur_x, line_end_x, search_nbr_range);
    line_Y = linspace(cur_y, line_end_y, search_nbr_range);
    line_X = max(line_X, 1);
    line_X = min(line_X, w);
    line_Y = max(line_Y, 1);
    line_Y = min(line_Y, h);            

    % mark T-junction
    reach_cids = EdgeGroupMap(sub2ind([h,w], round(line_Y), round(line_X)));
    reach_cids(reach_cids==cid)=0;
    inter_idx = find(reach_cids);
    if(~isempty(inter_idx))
        inter_idx = inter_idx(1);
        inter_x = line_X(inter_idx);
        inter_y = line_Y(inter_idx);
        reach_cid = reach_cids(inter_idx);

        attach_cfrag = cfrags{reach_cid}(:, 1:2)+1;
        attach_cfrag_idx = cfrags_idx{reach_cid};
        c_len = length(attach_cfrag_idx);
        dist2cfrag = attach_cfrag - repmat([inter_x, inter_y], c_len, 1);
        dist2cfrag = sqrt(dist2cfrag(:,1).^2 + dist2cfrag(:,2).^2);
        [temp_min_c_dist, r_idx] = min(dist2cfrag);

        add_end_edge_idx2 = attach_cfrag_idx(r_idx); 
        dist_2 = sqrt((cur_x-edges(add_end_edge_idx2,1)-1)^2 + (cur_y-edges(add_end_edge_idx2,2)-1)^2);
        
        if(dist_2 > 2*c_len)
            add_end_edge_idx2 = 0;
        end
    end
    
    % compare two choices and choose the shorter one
    c_type =0;
    if(add_end_edge_idx1 && ~add_end_edge_idx2)
        c_type = 1;
        add_end_edge_idx = add_end_edge_idx1;
    elseif(~add_end_edge_idx1 && add_end_edge_idx2)
        c_type = 2;
        add_end_edge_idx = add_end_edge_idx2;
    elseif(add_end_edge_idx1 && add_end_edge_idx2)
        if(dist_1 < dist_2)
            c_type = 1;
            add_end_edge_idx = add_end_edge_idx1;
        else
            c_type = 2;
            add_end_edge_idx = add_end_edge_idx2;
        end
    end
    
    
    if(c_type==1)
        if(is_start)
            cfrags{cid} = [edges(add_end_edge_idx,:); cfrags{cid}];
            cfrags_idx{cid} = [add_end_edge_idx cfrags_idx{cid}];
            G_test.fac(cid).nbrs_var(1)=reach_vids(1);
        else
            cfrags{cid} = [cfrags{cid}; edges(add_end_edge_idx,:)];
            cfrags_idx{cid} = [cfrags_idx{cid} add_end_edge_idx];            
            G_test.fac(cid).nbrs_var(2)=reach_vids(1);
        end
        G_test.var(reach_vids(1)).nbrs_fac = [G_test.var(reach_vids(1)).nbrs_fac; cid];
        G_test.var(vid).nbrs_fac(G_test.var(vid).nbrs_fac==cid) = [];
        % mark connected edge
        edges(add_end_edge_idx,5) = 2;
%         end_point_eid_map(cur_y, cur_x) = 0;
        plot([cur_x, edges(add_end_edge_idx,1)+1], [cur_y, edges(add_end_edge_idx,2)+1], 'g-');
    elseif(c_type==2)
        if(is_start)
            cfrags{cid} = [edges(add_end_edge_idx,:); cfrags{cid}];
            cfrags_idx{cid} = [add_end_edge_idx cfrags_idx{cid}];
        else
            cfrags{cid} = [cfrags{cid}; edges(add_end_edge_idx,:)];
            cfrags_idx{cid} = [cfrags_idx{cid} add_end_edge_idx];            
        end

        % Mark T-junction
        edges(add_end_edge_idx,5) = 3;
        jct_map(round(edges(add_end_edge_idx,2)+1), round(edges(add_end_edge_idx,1))+1) = 1; 
%         end_point_eid_map(cur_y, cur_x) = 0;
        plot([cur_x, edges(add_end_edge_idx,1)+1], [cur_y, edges(add_end_edge_idx,2)+1], 'g-');
    end

    
end


cfrags_idx = cfrags_idx(~cellfun('isempty',cfrags_idx));
cfrags = cfrags(~cellfun('isempty',cfrags)); 

