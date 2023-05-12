function [geom_diff, texture_diff] = degree_2_node_cues(c1,c2, texton_map, nbr_range_th)
if nargin<4
   nbr_range_th = 5; 
end
nbr_range_th_max = min([nbr_range_th size(c1,1)-1 size(c2,1)-1]);
nbr_range_th_min = min([nbr_range_th_max 2]);


nbr_range_vec = nbr_range_th_min:2:nbr_range_th_max;
num_range = length(nbr_range_vec);
num_edges = size(c1,1);

texton_nbr_range = 3;

V_1 = c1(:, 1:2);
V_1 = V_1 +1; % change CXX to the matlab coordinates
K_1 = LineCurvature2D(V_1);
N_1 = LineNormals2D(V_1);

V_2 = c2(:, 1:2);
V_2 = V_2 +1; % change CXX to the matlab coordinates
K_2 = LineCurvature2D(V_2);
N_2 = LineNormals2D(V_2);

% % tan direction of cfs at the node
% ori_1 = V_1(end, 1:2) - V_1(end-nbr_range_th+1, 1:2);
% ori_2 = V_2(nbr_range_th, 1:2) - V_2(1, 1:2);
% 
% geom_diff = ori_1*ori_2'/sqrt(ori_1*ori_1')/sqrt(ori_2*ori_2');

ori_1 = V_1(repmat(num_edges, [1, num_range]), 1:2) - V_1(repmat(num_edges, [1, num_range]) - nbr_range_vec, 1:2);
ori_2 = V_2(repmat(1, [1, num_range]) + nbr_range_vec, 1:2) - V_2(repmat(1, [1, num_range]), 1:2);

geom_diff = diag(ori_1*ori_2')./sqrt(diag(ori_1*ori_1'))./sqrt(diag(ori_2*ori_2'));
geom_diff = mean(geom_diff);

% compute the left and right nbr texton hist of each curve 
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

%%%%%%%%%%%%%%%%%%%% dist between left histograms
[hist_1,n1]=hist(left_texton_vec_1, 1:64);
hist_1=hist_1/length(left_texton_vec_1);

[hist_2,n2]=hist(left_texton_vec_2, 1:64);
hist_2=hist_2/length(left_texton_vec_2);
% keyboard;

left_dist = pdist2(hist_1,hist_2,'chisq');

%%%%%%%%%%%%%%%%%%%% dist between right histograms
[hist_1,n1]=hist(right_texton_vec_1, 1:64);
hist_1=hist_1/length(right_texton_vec_1);

[hist_2,n2]=hist(right_texton_vec_2, 1:64);
hist_2=hist_2/length(right_texton_vec_2);

right_dist = pdist2(hist_1,hist_2,'chisq');

texture_diff = left_dist + right_dist;
% keyboard;

end