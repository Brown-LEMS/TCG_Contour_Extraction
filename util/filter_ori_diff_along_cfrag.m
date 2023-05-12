function ori_diff_vec = filter_ori_diff_along_cfrag(cfrag, nbr_range_th)

nbr_range_th = min(nbr_range_th, 5); % does not make sense for too long range
nbr_range_vec = 2:2:nbr_range_th;
num_range = length(nbr_range_vec);
%%%%%%%% This filter keep the original sample density


len = size(cfrag,1);
ori_diff_vec = ones(len,1);

V = cfrag(:, 1:2);
V = V +1; % change CXX to the matlab coordinates

%% compute prob along each curve
%%%%%%%%%%%%%%%%%%  TODO: can be converted into matrix
%%%%%%%%%%%%%%%%%%  computation to get rid of loop
for i = nbr_range_th +1: len-nbr_range_th
    
    %%%%%%%%%%%% compute ori diff at the connecting points
%     [geom_diff, texture_diff] = degree_2_node_cues(cfrag(idx_c1, :),cfrag(idx_c2, :), tmap);
    ori_1 = V(repmat(i, [1, num_range]), 1:2) - V(repmat(i, [1, num_range]) - nbr_range_vec, 1:2);
    ori_2 = V(repmat(i, [1, num_range]) + nbr_range_vec, 1:2) - V(repmat(i, [1, num_range]), 1:2);

    geom_diff = diag(ori_1*ori_2')./sqrt(diag(ori_1*ori_1'))./sqrt(diag(ori_2*ori_2'));
    ori_diff_vec(i) = mean(geom_diff);
end

ori_diff_vec = acos(ori_diff_vec);
end


