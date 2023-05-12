function [G, merged_cem, merged_cf_idx] = merge_at_degree_2_node(G, merged_cem, merged_cf_idx, v, actual_edge, nbrs_fac, c1_ids, c2_ids)

%%%%%%%%%%%%%%%%%  input c1_ids, c2_ids are NOT CUT portion

    c1_2merge = merged_cem{nbrs_fac(1)};
    c2_2merge = merged_cem{nbrs_fac(2)};
    c_merged = merge_two_curve_fragments(c1_2merge, c2_2merge, actual_edge);
    c_ids_merged = [c1_ids c2_ids(2:end)];
    % update c1 in G
    merged_cem{nbrs_fac(1)} = c_merged;
    merged_cf_idx{nbrs_fac(1)} = c_ids_merged;
    % remove c2 in G
    merged_cem{nbrs_fac(2)} = [];
    merged_cf_idx{nbrs_fac(2)} = [];
    % update merged node in G
    G.var(v).nbrs_fac = [];
    G.var(v).merged = 1;

    % n1 -c1 -> n2 - c2 -> n3,  n2 = G.var(v)
    % update nbr_var in c 1, make c2 removed
    nbr_var_id_1 = G.fac(nbrs_fac(1)).nbrs_var;
    n1_id = nbr_var_id_1(find(nbr_var_id_1~=G.var(v).id));

    nbr_var_id_2 = G.fac(nbrs_fac(2)).nbrs_var;
    n3_id = nbr_var_id_2(find(nbr_var_id_2~=G.var(v).id));

%                 if(isempty(n1_id)||isempty(n3_id))
%                     keyboard;
%                 end
    G.fac(nbrs_fac(1)).nbrs_var = [n1_id n3_id];
    G.fac(nbrs_fac(2)).removed = 1;
    % update nbr_fac of n3, replace c2_id with c2_id
    G.var(n3_id).nbrs_fac(find(G.var(n3_id).nbrs_fac == nbrs_fac(2))) = nbrs_fac(1);

end