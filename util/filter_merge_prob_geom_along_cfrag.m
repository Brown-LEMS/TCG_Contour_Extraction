function merge_prob_vec = filter_merge_prob_geom_along_cfrag(cfrag, nbr_range_th, beta_1, fmean_1)

nbr_range_th = min(nbr_range_th, 5); % does not make sense for too long range
%%%%%%%% This filter keep the original sample density


len = size(cfrag,1);
merge_prob_vec = ones(len,1);

V = cfrag(:, 1:2);
V = V +1; % change CXX to the matlab coordinates

%% compute prob along each curve
%%%%%%%%%%%%%%%%%%  TODO: can be converted into matrix
%%%%%%%%%%%%%%%%%%  computation to get rid of loop
for i = nbr_range_th +1: len-nbr_range_th
    
    %%%%%%%%%%%% compute ori diff at the connecting points
%     [geom_diff, texture_diff] = degree_2_node_cues(cfrag(idx_c1, :),cfrag(idx_c2, :), tmap);
    ori_1 = V(i, 1:2) - V(i-nbr_range_th, 1:2);
    ori_2 = V(i+nbr_range_th, 1:2) - V(i, 1:2);

    geom_diff = ori_1*ori_2'/sqrt(ori_1*ori_1')/sqrt(ori_2*ori_2');
    
    
    %%%%%%%%%%% compute the probability given difference between cues
    p = 1 / (1 + exp(-([1 geom_diff] - fmean_1)*beta_1'));
     merge_prob_vec(i) = p;
end

end


