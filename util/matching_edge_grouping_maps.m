function [cntP, sumP, cntR, sumR, match_vis1, match_vis2] = matching_edge_grouping_maps(EdgeGroupMap1, EdgeGroupMap2, maxDist)
% taking EdgeGroupMap1 as ground-truth, EdgeGroupMap2 as test
% match_vis1 = zeros([size(EdgeGroupMap1),3]);
% match_vis1(:,:, 1) = EdgeGroupMap1>0;
% match_vis2 = zeros([size(EdgeGroupMap2),3]);
% match_vis2(:,:, 1) = EdgeGroupMap2>0;
match_vis1 = repmat(EdgeGroupMap1>0, [1,1,3]);
match_vis1 =  1 - match_vis1;
match_vis2 = repmat(EdgeGroupMap2>0, [1,1,3]);
match_vis2 =  1 - match_vis2;

labels_1 = unique(EdgeGroupMap1);
labels_2 = unique(EdgeGroupMap2);
[h, w] = size(EdgeGroupMap1);

sumR = sum(EdgeGroupMap1(:)>0);
sumP = sum(EdgeGroupMap2(:)>0);
% sumP = 0;
% sumR = 0;
cntP = 0;
cntR = 0;

%%%%%%%%%% process for precision
Edgemap1 = EdgeGroupMap1>0;
for li = 2:length(labels_2)
    l = labels_2(li);
    [Y, X] = find(EdgeGroupMap2==l);
    y_min = max(min(Y)-5, 1);
    y_max = min(max(Y)+5, h);
    x_min = max(min(X)-5, 1);
    x_max = min(max(X)+5, w);
    Edgemap2_crop = (EdgeGroupMap2(y_min:y_max, x_min:x_max)==l);
    Edgemap1_crop = Edgemap1(y_min:y_max, x_min:x_max);
    EdgeGroupMap1_crop = EdgeGroupMap1(y_min:y_max, x_min:x_max);
    
    [match1,match2] = correspondPixels(double(Edgemap1_crop>0),double(Edgemap2_crop>0),maxDist);
    
%     sumP = sumP + sum(Edgemap2(:)>0);
    
    match_coords = match2(match2(:)>0);
    match_group_labels = EdgeGroupMap1_crop(match_coords);
    most_repeat_label = mode(match_group_labels);
    
    cnt_match = sum(match_group_labels==most_repeat_label(1));
    cntP = cntP + cnt_match;
%     cntR = cntR + cnt_match;
%     sumR = sumR + sum(EdgeGroupMap1(:)==most_repeat_label(1));
%     match_coords_sub = match_coords(match_group_labels==most_repeat_label(1));
%     [h_c,w_c] = size(match2);
%     [match_y, match_x] =  ind2sub([h_c,w_c], match_coords_sub);
    [match_y, match_x] = find(match2>0);
    match_y = match_y(match_group_labels==most_repeat_label(1));
    match_x = match_x(match_group_labels==most_repeat_label(1));
    match_y = match_y + y_min -1;
    match_x = match_x + x_min -1;
    match_vis2(sub2ind([h,w,3], match_y, match_x, 2*ones(size(match_x))))=1;
    
%     keyboard;
end

%%%%%%%%%% reverve process for recall
Edgemap2 = EdgeGroupMap2>0;
for li = 2:length(labels_1)
    l = labels_1(li);
    [Y, X] = find(EdgeGroupMap1==l);
    
%     if(sum(ismember([Y, X],[86, 178],'rows')))
%         keyboard;
%     end
    y_min = max(min(Y)-5, 1);
    y_max = min(max(Y)+5, h);
    x_min = max(min(X)-5, 1);
    x_max = min(max(X)+5, w);
    Edgemap1_crop = (EdgeGroupMap1(y_min:y_max, x_min:x_max)==l);
    Edgemap2_crop = Edgemap2(y_min:y_max, x_min:x_max);
    EdgeGroupMap2_crop = EdgeGroupMap2(y_min:y_max, x_min:x_max);
    
    [match1,match2] = correspondPixels(double(Edgemap1_crop>0),double(Edgemap2_crop>0),maxDist);
    
%     sumP = sumP + sum(Edgemap2(:)>0);
    
    match_coords = match1(match1(:)>0);
    match_group_labels = EdgeGroupMap2_crop(match_coords);
    most_repeat_label = mode(match_group_labels);
    
    cnt_match = sum(match_group_labels==most_repeat_label(1));
%     cntP = cntP + cnt_match;
    cntR = cntR + cnt_match;
%     sumR = sumR + sum(EdgeGroupMap1(:)==most_repeat_label(1));
    
    [match_y, match_x] = find(match1>0);
    match_y = match_y(match_group_labels==most_repeat_label(1));
    match_x = match_x(match_group_labels==most_repeat_label(1));
    match_y = match_y + y_min -1;
    match_x = match_x + x_min -1;
    match_vis1(sub2ind([h,w,3], match_y, match_x, 2*ones(size(match_x))))=1;
end

