function  [G, merged_cfrags, merged_cfrags_idx, T_junctions] = classify_junction_type_wrt_graph_BP(new_cfrags, new_cfrags_idx, edges, params)    
% this version consider more general case of junction with >=3 connected
% curve fragments
w0 = 1; % preference of global reasoning: 1 means no preference
angle_diff_th = params.BP_merge_angle_th;
nbr_num_edges = params.BP_nbr_num_edges;
clen_th = params.BP_clen_th;
disp('classify junction type wrt topological contour graph...');    
T_junctions = [];
    
% construct factor graph based on curve fragment candidates
G = construct_fac_graph_from_curve_fragments (new_cfrags_idx, new_cfrags);

merged_cfrags = new_cfrags;
merged_cfrags_idx = new_cfrags_idx;
  
%%  compute merge probability at each junction node
for vid = 1:length(G.var)

%     actual_edge = edges(G.var(vid).actual_edge_id,:);
%     if(round(actual_edge(1))==142 && round(actual_edge(2))==164)
%         keyboard;
%     end
%     if(round(actual_edge(1))==424 && round(actual_edge(2))==421)
%         keyboard;
%     end
    nbrs_fac = G.var(vid).nbrs_fac;
    v_degree = length(nbrs_fac);
    if( v_degree<3) % dim 3 nodes
        continue;
    end

    % skip the case including cycles
    if(length(unique(nbrs_fac))<v_degree)
        continue;
    end

    
    % only consider portion of the curve fragment within
    % neighbourhood of the connecting node
    % always make directions c1 -> node <- c2
    %                                ^
    %                                |
    %                               c3
    %
    cfrags_at_v = cell(1, v_degree);

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
%             angle_diff_mat(j,i) = angle_diff_mat(i,j);
        end

    end
    % merge decision refer to this vector, saved order: (c1,c2), (c1,c3), (c2,c3)
    probs = exp(-angle_diff_mat);
    G.var(vid).p = probs;
    G.var(vid).angle_diff_mat = angle_diff_mat;
    G.var(vid).cfrags_at_v = cfrags_at_v;
    
end

%%%%%%%%%%%%%%% One iteration of belief propagation, a contour receives msg
%%%%%%%%%%%%%%% from both connected junctions and send back update msg
for cid = 1:length(G.fac)
   
    % only consider the case when both ends are junctions
    vid1 = G.fac(cid).nbrs_var(1);
    vid2 = G.fac(cid).nbrs_var(2);
    nbrs_fac1 = G.var(vid1).nbrs_fac;
    v_degree1 = length(nbrs_fac1);
    if( v_degree1<3)
        continue;
    end
    % skip the case including cycles
    if(length(unique(nbrs_fac1))<v_degree1)
        continue;
    end
    nbrs_fac2 = G.var(vid2).nbrs_fac;
    v_degree2 = length(nbrs_fac2);
    if( v_degree2<3)
        continue;
    end
    % skip the case including cycles
    if(length(unique(nbrs_fac2))<v_degree2)
        continue;
    end
    % skip double link case
    if(length(intersect(nbrs_fac1,nbrs_fac2))>1)
        continue;
    end
    
    % skip long curves, no msg could pass long curve
    % use the actual cfrag for the connection
    c0 = merged_cfrags{cid}; 
    c_len = contour_length_mex(c0');    
    if(c_len > clen_th)
        continue;
    end
    
    % receive msg from both junctions, and send msg back, 
    % update merge probablities at each junction
    %           |           |
    %           c1          c2
    %           |           |
    %           v           v
    % ---c1--> J1 ---c0---> J2 <--c2---
    
    % update merg prob for J1
    ind_c = find(nbrs_fac1 == cid);
    for i = 1:v_degree1
        cid1 = nbrs_fac1(i);
        if(cid1 == cid)
            continue;
        end
        min_angle_diff = pi;
        % use the previous saved cfrag
        c1 = G.var(vid1).cfrags_at_v{i};
        % consider the co-circularity between c1, and merged c0,c2
        for j = 1:v_degree2
            cid2 = nbrs_fac2(j);
            if(cid2 == cid)
                continue;
            end
            % use the previous saved cfrag
            c2 = G.var(vid2).cfrags_at_v{j};
            c_merge = [c2(1:end-1,:); flip(c0,1)]; 
            % input need to be ---c--> v <---c----           
            angle_diff = co_circular_cost(c1, c_merge)*w0; % add preference
            min_angle_diff = min(min_angle_diff, angle_diff);
        end
        
        % update minimum angle diff
        if(i<ind_c)
            G.var(vid1).angle_diff_mat(i,ind_c) = min(G.var(vid1).angle_diff_mat(i,ind_c), min_angle_diff);
        else
            G.var(vid1).angle_diff_mat(ind_c,i) = min(G.var(vid1).angle_diff_mat(ind_c,i), min_angle_diff);
        end
    end
    
    % update merg prob for J2
    ind_c = find(nbrs_fac2 == cid);
    for j = 1:v_degree2
        cid2 = nbrs_fac2(j);
        if(cid2 == cid)
            continue;
        end
        min_angle_diff = pi;
        % use the previous saved cfrag
        c2 = G.var(vid2).cfrags_at_v{j};
        % consider the co-circularity between c1, and merged c0,c2
        for i = 1:v_degree1
            cid1 = nbrs_fac1(i);
            if(cid1 == cid)
                continue;
            end
            % use the previous saved cfrag
            c1 = G.var(vid1).cfrags_at_v{i};
            c_merge = [c1(1:end-1, :); c0]; 
            % input need to be ---c--> v <---c----
            angle_diff = co_circular_cost(c2, c_merge)*w0; % add preference
            min_angle_diff = min(min_angle_diff, angle_diff);
        end
        
        % update minimum angle diff
        if(j<ind_c)
            G.var(vid2).angle_diff_mat(j,ind_c) = min(G.var(vid2).angle_diff_mat(j,ind_c), min_angle_diff);
        else
            G.var(vid2).angle_diff_mat(ind_c,j) = min(G.var(vid2).angle_diff_mat(ind_c,j), min_angle_diff);
        end
    end
    
end

%%%%%%%%%%%%%% Decision of T-Junctions based on updated merging probs
% Another iteration to merge, this get rid of the problem of merging order
cid_to_remove = [];
for vid = 1:length(G.var)

    nbrs_fac = G.var(vid).nbrs_fac;
    v_degree = length(nbrs_fac);
    if( v_degree~=3) % dim 3 nodes
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
%     if(sum(angle_diff_mat(:)<angle_diff_th)>1)
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
% update factor graph based on curve fragment candidates
G = construct_fac_graph_from_curve_fragments (merged_cfrags_idx, merged_cfrags);

