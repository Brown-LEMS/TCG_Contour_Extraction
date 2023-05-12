function  [G, merged_cfrags, merged_cfrags_idx, T_junctions] = classify_junction_type_wrt_graph2(new_cfrags, new_cfrags_idx, edges, params)    
% this version consider more general case of junction with >=3 connected
% curve fragments
angle_diff_th = pi/12;
nbr_num_edges = params.nbr_num_edges;
% nbr_num_edges = 20;
disp('classify junction type wrt topological contour graph...');    
T_junctions = [];
    
% construct factor graph based on curve fragment candidates
G = construct_fac_graph_from_curve_fragments (new_cfrags_idx, new_cfrags);

merged_cfrags = new_cfrags;
merged_cfrags_idx = new_cfrags_idx;
  
%%  decide dim 3 nodes merging
for vid = 1:length(G.var)

    nbrs_fac = G.var(vid).nbrs_fac;
    v_degree = length(nbrs_fac);
    if( v_degree<3) % dim 3 nodes
        continue;
    end

    % skip the case including cycles
    if(length(unique(nbrs_fac))<v_degree)
        continue;
    end
%     actual_edge = edges(G.var(vid).actual_edge_id,:);


    % only consider portion of the curve fragment within
    % neighbourhood of the connecting node
    % always make directions c1 -> node <- c2
    %                                ^
    %                                |
    %                               c3
    %
    % length of contours taken into account must be equal
    
    cfrags_at_v = cell(1, v_degree);

%     cur_nbr = nbr_num_edges;
%     for i = 1:length(nbrs_fac)
%         cur_nbr =  min(cur_nbr, length(merged_cfrags_idx{nbrs_fac(i)}));
%     end

    for i = 1:length(nbrs_fac)

        cid = nbrs_fac(i);
        c1_ids = merged_cfrags_idx{cid};
        c1 = merged_cfrags{cid};
        cur_nbr = min(nbr_num_edges, length(c1_ids));
        
        if(c1_ids(1) == G.var(vid).actual_edge_id)
            c1 = c1(1:cur_nbr, :);
            % reverse the order of edges
            c1 = fliplr(c1');
            c1 = c1';
%             c1_ids = fliplr(c1_ids);
        elseif (c1_ids(end) == G.var(vid).actual_edge_id)
            c1 = c1(end-cur_nbr+1:end, :);
        end  
        cfrags_at_v{i} = c1;
%         cfrags_idx_at_v{i} = c1_ids;
    end

    %%%%%%%%%%%%%%%  compute the pairwise diff and merge probs
    angle_diff_mat = ones(v_degree, v_degree) * pi;

    for i = 1:v_degree
        for j= i+1:v_degree
            angle_diff_mat(i,j) = co_circular_cost(cfrags_at_v{i}, cfrags_at_v{j});
        end

    end
    % merge decision refer to this vector, saved order: (c1,c2), (c1,c3), (c2,c3)
    probs = exp(-angle_diff_mat);
    G.var(vid).p = probs;
    G.var(vid).angle_diff_mat = angle_diff_mat;
    
end


% Another iteration to merge, this get rid of the problem of merging order
cid_to_remove = [];
for vid = 1:length(G.var)

    nbrs_fac = G.var(vid).nbrs_fac;
    v_degree = length(nbrs_fac);
    if( v_degree<3) % dim 3 nodes
        continue;
    end

    % skip the case including cycles
    if(length(unique(nbrs_fac))<v_degree)
        continue;
    end
    actual_edge = edges(G.var(vid).actual_edge_id,:);

    cfrags_idx_at_v = cell(1, v_degree);
    for i = 1:length(nbrs_fac)

        cid = nbrs_fac(i);
        c1_ids = merged_cfrags_idx{cid};

        if(c1_ids(1) == G.var(vid).actual_edge_id)
            c1_ids = fliplr(c1_ids);
        end  

        cfrags_idx_at_v{i} = c1_ids;
    end
    
    angle_diff_mat = G.var(vid).angle_diff_mat;
    is_merge = 0;
    min_angle_diff = min(angle_diff_mat(:));
    [min_i, min_j] = find(angle_diff_mat==min_angle_diff);
    if(length(min_i)>1) % pass the case with ambiguity
        continue;
    end
    if(min_angle_diff < angle_diff_th)
        is_merge = 1;
    end

    % merge the curve fragments and update G
    if(is_merge)
        T_junctions = [T_junctions; actual_edge];
        cid_1 = nbrs_fac(min_i);
        cid_2 = nbrs_fac(min_j);

        % merge the curve fragments
        c1_ids = cfrags_idx_at_v{min_i};
        c2_ids = cfrags_idx_at_v{min_j};
        c2_ids = flip(c2_ids);
        cfrag_idx = [c1_ids, c2_ids(2:end)]; 
        cfrag = edges(cfrag_idx,:);
        merged_cfrags{cid_1} = cfrag;
        merged_cfrags_idx{cid_1} = cfrag_idx;

        % update G: change
        % v1 ---c1----- v ----- c2------v2
        %               to
        % v1 -----------c1--------------v2
        % 
        % detach c1, c2 from v
        G.var(vid).nbrs_fac([min_i, min_j]) = [];
        % locate v2
        vid_2 = G.fac(cid_2).nbrs_var( G.fac(cid_2).nbrs_var ~= vid);
        % link c1 to v2
        G.fac(cid_1).nbrs_var( G.fac(cid_1).nbrs_var == vid) = vid_2;
        G.var(vid_2).nbrs_fac(G.var(vid_2).nbrs_fac==cid_2) = cid_1;
        % right now, c2 is isolated, do not need to worry, wait to be
        % deleted at last
        cid_to_remove = [cid_to_remove cid_2];

    end    
end

merged_cfrags(cid_to_remove) = [];
merged_cfrags_idx(cid_to_remove) = [];
% only keep the non-empty contours
merged_cfrags = merged_cfrags(~cellfun('isempty',merged_cfrags));
merged_cfrags_idx = merged_cfrags_idx(~cellfun('isempty',merged_cfrags_idx)); 
