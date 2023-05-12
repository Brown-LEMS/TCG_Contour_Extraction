function [G, merged_cem, merged_cf_idx] = merge_at_degree_3_node(G, merged_cem, merged_cf_idx, v, actual_edge, c1_ids, c2_ids, c3_ids)

%%%%             directions c1 -> node <- c2
%%%%                                ^
%%%%                                |
%%%%                               c3
%%%%
%%%%%%%%%%%%%%%%  merge the 2 cfrags with highest merge prob at degree 3 node
%%%%%%%%%%%%%%%%  input c1_ids, c2_ids, c3_ids are NOT CUT portion


    [max_p, max_p_id] = max(G.var(v).p);

    switch max_p_id
        case 1
                initial_cid = G.var(v).nbrs_fac(1);
                merged_cid = G.var(v).nbrs_fac(2);
                untouch_cid = G.var(v).nbrs_fac(3);
               c_ids_merged = [c1_ids c2_ids(end-1:-1:1)];
        case 2
                initial_cid = G.var(v).nbrs_fac(1);
                merged_cid = G.var(v).nbrs_fac(3);
                untouch_cid = G.var(v).nbrs_fac(2);
                c_ids_merged = [c1_ids c3_ids(end-1:-1:1)];       
        case 3
                initial_cid = G.var(v).nbrs_fac(2);
                merged_cid = G.var(v).nbrs_fac(3);
                untouch_cid = G.var(v).nbrs_fac(1);
                c_ids_merged = [c2_ids c3_ids(end-1:-1:1)];
    end

%%%% retrieve the orignal version of cfrags
    c1_2merge = merged_cem{initial_cid};
    c2_2merge = merged_cem{merged_cid};
    c_merged = merge_two_curve_fragments(c1_2merge, c2_2merge, actual_edge);
    % update c1 in G
    merged_cem{initial_cid} = c_merged;
    merged_cf_idx{initial_cid} = c_ids_merged;
    % remove c2 in G
    merged_cem{merged_cid} = [];
    merged_cf_idx{merged_cid} = [];
    % update merged node in G
    G.var(v).nbrs_fac = untouch_cid;
    G.var(v).merged = 1;

    % n1 -c1 -> n2 - c2 -> n3,  n2 = G.var(v)
    % update nbr_var in c 1, make c2 removed
    nbr_var_id_1 = G.fac(initial_cid).nbrs_var;
    n1_id = nbr_var_id_1(find(nbr_var_id_1~=G.var(v).id));

    nbr_var_id_2 = G.fac(merged_cid).nbrs_var;
    n3_id = nbr_var_id_2(find(nbr_var_id_2~=G.var(v).id));

%     if(isempty(n1_id)||isempty(n3_id))
%         keyboard;
%     end
    G.fac(initial_cid).nbrs_var = [n1_id n3_id];
    G.fac(merged_cid).removed = 1;
    % update nbr_fac of n3, replace c2_id with c2_id
    G.var(n3_id).nbrs_fac(find(G.var(n3_id).nbrs_fac == merged_cid)) = initial_cid;


end