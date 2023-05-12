function T_junctions = locate_T_jct_from_graph(G, total_num_edges)
% input: 
%   G: constructed from construct_fac_graph_from_curve_fragments_v2.m
% output: 
%   T_junctions: each row record a T-junction in
%                [edge_id, stem contour id, reached contour id]
T_junctions = [];
map_eid_cid = zeros(1, total_num_edges);

%% Iterate each graph node, locate end point (degre==1), associated its contour id with edge id
for v = 1:length(G.var)
    jct_id = G.var(v).nbrs_fac;
    if(length(jct_id) ~= 1)
        continue;
    end
    map_eid_cid(G.var(v).actual_edge_id) = jct_id;
end


%% Iterate each contour, look for thoses edge id in the middle belonging to located T junctions 
for f = 1:length(G.fac)
   cfrag_idx = G.fac(f).cfrag_idx;
   for i = 2:length(cfrag_idx)-1
       if ( map_eid_cid(cfrag_idx(i)) > 0 )
           T_junction = [cfrag_idx(i), map_eid_cid(cfrag_idx(i)), f];
           T_junctions = [T_junctions; T_junction];
       end
   end
end

