function [cues_diff, probs] = degree_3_node_cues(c1,c2,c3, texton_map, hsv_img, edgemap, fmean_0, beta_0, fmean_1, beta_1, nbr_len_th, nbr_range_th)
%             directions c1 -> node <- c2
%                                ^
%                                |
%                               c3
%
%  all the diff are vectors in the order of  1 vs 2,  1 vs 3,  2 vs 3

geom_diff = zeros(3,1);
texture_diff = zeros(3,1);
other_cues_diff = zeros(3,6);
num_pairs = 3;


texton_nbr_range = 3;

V_1 = c1(:, 1:2);
V_1 = V_1 +1; % change CXX to the matlab coordinates
K_1 = LineCurvature2D(V_1);
N_1 = LineNormals2D(V_1);


V_2 = c2(:, 1:2);
V_2 = V_2 +1; % change CXX to the matlab coordinates
K_2 = LineCurvature2D(V_2);
N_2 = LineNormals2D(V_2);

V_3 = c3(:, 1:2);
V_3 = V_3 +1; % change CXX to the matlab coordinates
K_3 = LineCurvature2D(V_3);
N_3 = LineNormals2D(V_3);

% tan direction of cfs at the node
nbr_range_th = min([nbr_range_th size(V_1,1) size(V_2,1) size(V_3,1)]);
ori_1 = V_1(end, 1:2) - V_1(end-nbr_range_th+1, 1:2);
ori_2 = V_2(end, 1:2) - V_2(end-nbr_range_th+1, 1:2);
ori_3 = V_3(end, 1:2) - V_3(end-nbr_range_th+1, 1:2);

%%%%%%%%%%%%%%%%%%%%%%%%% compute the orientation diffs at the node
%%%%%% the direction matters
geom_diff(1) = -ori_1*ori_2'/sqrt(ori_1*ori_1')/sqrt(ori_2*ori_2');
geom_diff(2) = -ori_1*ori_3'/sqrt(ori_1*ori_1')/sqrt(ori_3*ori_3');
geom_diff(3) = -ori_2*ori_3'/sqrt(ori_2*ori_2')/sqrt(ori_3*ori_3');



%%%%%%%%%%%%%%%%%%%%%% compute the left and right nbr texton hist of each curve 
[h w ~] = size(texton_map);

left_texton_vec_1 = [];
right_texton_vec_1 = [];

for local_dist = 1:texton_nbr_range
    [left_x, left_y, right_x, right_y] = get_left_right_coordinates(V_1, N_1, local_dist, w, h);
    left_textons = diag(texton_map(left_y, left_x, 1));
    right_textons = diag(texton_map(right_y, right_x, 1));
    left_texton_vec_1 = [left_texton_vec_1; left_textons];
    right_texton_vec_1 = [right_texton_vec_1; right_textons];
end

left_texton_vec_2 = [];
right_texton_vec_2 = [];

for local_dist = 1:texton_nbr_range
    [left_x, left_y, right_x, right_y] = get_left_right_coordinates(V_2, N_2, local_dist, w, h);
    left_textons = diag(texton_map(left_y, left_x, 1));
    right_textons = diag(texton_map(right_y, right_x, 1));
    left_texton_vec_2 = [left_texton_vec_2; left_textons];
    right_texton_vec_2 = [right_texton_vec_2; right_textons];
end

left_texton_vec_3 = [];
right_texton_vec_3 = [];

for local_dist = 1:texton_nbr_range
    [left_x, left_y, right_x, right_y] = get_left_right_coordinates(V_3, N_3, local_dist, w, h);
    left_textons = diag(texton_map(left_y, left_x, 1));
    right_textons = diag(texton_map(right_y, right_x, 1));
    left_texton_vec_3 = [left_texton_vec_3; left_textons];
    right_texton_vec_3 = [right_texton_vec_3; right_textons];
end

%%%%%%%%%%%%%%%%%%%% dist between left/right histograms
left_hist_1=hist(left_texton_vec_1, 1:64)/length(left_texton_vec_1);
right_hist_1=hist(right_texton_vec_1, 1:64)/length(right_texton_vec_1);
left_hist_2=hist(left_texton_vec_2, 1:64)/length(left_texton_vec_2);
right_hist_2=hist(right_texton_vec_2, 1:64)/length(right_texton_vec_2);
left_hist_3=hist(left_texton_vec_3, 1:64)/length(left_texton_vec_3);
right_hist_3=hist(right_texton_vec_3, 1:64)/length(right_texton_vec_3);


texture_diff(1) = pdist2(left_hist_1,right_hist_2,'chisq') + pdist2(right_hist_1,left_hist_2,'chisq');
texture_diff(2) = pdist2(left_hist_1,right_hist_3,'chisq') + pdist2(right_hist_1,left_hist_3,'chisq');
texture_diff(3) = pdist2(left_hist_2,right_hist_3,'chisq') + pdist2(right_hist_2,left_hist_3,'chisq');

%%%%%%%%%%%%%%%% compute the cues diff between portion of conneting curve fragments
[bg_grad, sat_grad, hue_grad, abs_k, edge_sparsity, wigg, c_len1] = curve_fragment_cues(c1, hsv_img, edgemap);
c1_cues = [bg_grad sat_grad hue_grad abs_k edge_sparsity wigg];
[bg_grad, sat_grad, hue_grad, abs_k, edge_sparsity, wigg, c_len2] = curve_fragment_cues(c2, hsv_img, edgemap);
c2_cues = [bg_grad sat_grad hue_grad abs_k edge_sparsity wigg];
[bg_grad, sat_grad, hue_grad, abs_k, edge_sparsity, wigg, c_len3] = curve_fragment_cues(c3, hsv_img, edgemap);
c3_cues = [bg_grad sat_grad hue_grad abs_k edge_sparsity wigg];

other_cues_diff(1,:) = abs(c1_cues-c2_cues.*[-1 -1 -1 1 1 1]);
other_cues_diff(2,:) = abs(c1_cues-c3_cues.*[-1 -1 -1 1 1 1]);
other_cues_diff(3,:) = abs(c2_cues-c3_cues.*[-1 -1 -1 1 1 1]);

cues_diff = [other_cues_diff geom_diff texture_diff];

probs = ones(num_pairs, 1) ./ (ones(num_pairs, 1) + exp(-([ones(num_pairs, 1) cues_diff]- repmat(fmean_0, [num_pairs, 1]))*beta_0'));

% if(c_len1< nbr_len_th )
%     probs(1) = min(probs(1), 1 / (1 + exp(-([1 geom_diff(1)] - fmean_1)*beta_1')));
%     probs(2) = min(probs(2), 1 / (1 + exp(-([1 geom_diff(2)] - fmean_1)*beta_1')));
% end
% 
% if(c_len2< nbr_len_th )
%     probs(1) = min(probs(1), 1 / (1 + exp(-([1 geom_diff(1)] - fmean_1)*beta_1')));
%     probs(3) = min(probs(3), 1 / (1 + exp(-([1 geom_diff(3)] - fmean_1)*beta_1')));
% end
% 
% if(c_len3< nbr_len_th )
%     probs(2) = min(probs(2), 1 / (1 + exp(-([1 geom_diff(2)] - fmean_1)*beta_1')));
%     probs(3) = min(probs(3), 1 / (1 + exp(-([1 geom_diff(3)] - fmean_1)*beta_1')));
% end