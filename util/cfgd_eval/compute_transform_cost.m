% compute the transformation cost from a group of curve fragment to the
% other
function dist = compute_transform_cost(combo_1, combo_2, CFG_1, CFG_2, local_dist)
ind_1 = combo_1.ind;
ind_2 = combo_2.ind;
frags_1 = CFG_1(ind_1);
frags_2 = CFG_2(ind_2);
dist = 0;

combo_length_1 = 0;
combo_length_2 = 0;
all_edges_1 = [];
all_edges_2 = [];

% compute the total length of contours in combo_1
for i = 1:size(frags_1, 2)
    combo_length_1 = combo_length_1 + contour_length_mex(frags_1{1,i}');
    all_edges_1 = [all_edges_1; frags_1{1,i}];
end

% compute the total length of contours in combo_2
for i = 1:size(frags_2, 2)
    combo_length_2 = combo_length_2 + contour_length_mex(frags_2{1,i}');
    all_edges_2 = [all_edges_2; frags_2{1,i}];
end

dist = compute_transform_cost_mex(all_edges_1', all_edges_2', local_dist);
% Add a factor of the number of the frags and perform normalization
dist = (dist + 2*(size(frags_1, 2) + 2*size(frags_2, 2)-2))/(combo_length_1+combo_length_2);
% dist = dist/(combo_length_1+combo_length_2);
end