% Function to determin if combo_2 overlap a considerable potion of combo_1
function b = overlap_combos(combo_1, combo_2, CFG, local_dist)

ind_1 = combo_1.ind;
ind_2 = combo_2.ind;
frags_1 = CFG(ind_1);
frags_2 = CFG(ind_2);

all_edges_1 = [];
all_edges_2 = [];

% build up all edges list
for i = 1:size(frags_1, 2)
    all_edges_1 = [all_edges_1; frags_1{1,i}];
end

% build up all edges list
for i = 1:size(frags_2, 2)
    all_edges_2 = [all_edges_2; frags_2{1,i}];
end

b = overlap_combos_mex(all_edges_1', all_edges_2', local_dist);

end