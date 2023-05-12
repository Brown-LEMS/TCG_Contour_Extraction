function adjust_cood_cem (cem_file_in, cem_file_out, adjust)
addpath io/

[CEM, edges, cf_idx] = load_contours(cem_file_in);
h = CEM{1}(2);
w = CEM{1}(1);

% adjust edges
edges(:,1) = edges(:,1) + adjust;
edges(:,2) = edges(:,2) + adjust;

new_cfrags = CEM{2};

% adjsut each curve fragment
for c = 1:length(new_cfrags)
    cfrag = new_cfrags{c};
    cfrag(:,1) = cfrag(:,1) + adjust;
    cfrag(:,2) = cfrag(:,2) + adjust;
    new_cfrags{c} = cfrag;
end

write_cem_fixed_edge_id(cem_file_out, new_cfrags, edges, h, w)


end
