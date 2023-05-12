function G = construct_fac_graph_from_curve_fragments_v2 (cf_idx, contours)
% construct a graph G: graph node represents connecting curve fragments end
% points, graph factor represents each curve fragments

G = init_graph();
cf_idx = cf_idx(~cellfun('isempty',cf_idx));

used_edge_ids = [];
for c = 1:length(cf_idx)
        
    search_1 = find(used_edge_ids==cf_idx{c}(1));

    if(isempty(search_1))
        [G, id_1] = add_varnode(G, cf_idx{c}(1));
        used_edge_ids = [used_edge_ids cf_idx{c}(1)];
    else
        id_1 = search_1(1);
    end
    
    search_2 = find(used_edge_ids==cf_idx{c}(end));
    if(isempty(search_2))
        [G, id_2]= add_varnode(G, cf_idx{c}(end));
        used_edge_ids = [used_edge_ids cf_idx{c}(end)];
    else
        id_2 = search_2(1);
    end
    
    G = add_facnode_v2(G, cf_idx{c}, contours{c}, [id_1 id_2]);
end

end