function [cntP, sumP, cntR, sumR] = matching_edge_grouping_maps_One2Multi(EdgeGroupMap_gt, EdgeGroupMap_cp, maxDist)
% match_vis1 = repmat(EdgeGroupMap_gt>0, [1,1,3]);
% match_vis1 =  1 - match_vis1;
% match_vis2 = repmat(EdgeGroupMap_cp>0, [1,1,3]);
% match_vis2 =  1 - match_vis2;

labels_1 = unique(EdgeGroupMap_gt);
labels_2 = unique(EdgeGroupMap_cp);
[h, w] = size(EdgeGroupMap_gt);

Edgemap_gt = EdgeGroupMap_gt>0;
Edgemap_cp = EdgeGroupMap_cp>0;
[match1,match2] = correspondPixels(double(Edgemap_gt>0),double(Edgemap_cp>0),maxDist);

sumR = sum(EdgeGroupMap_gt(:)>0);
cntR = sum(match1(:)>0);
% sumP = sum(EdgeGroupMap_cp(:)>0);
sumP = 0;
cntP = 0;


%%%%%%%%%% process for precision
for li = 2:length(labels_2)
    l = labels_2(li);
    [Y, X] = find(EdgeGroupMap_cp==l);
    y_min = max(min(Y)-5, 1);
    y_max = min(max(Y)+5, h);
    x_min = max(min(X)-5, 1);
    x_max = min(max(X)+5, w);
    Edgemap_cp_crop = (EdgeGroupMap_cp(y_min:y_max, x_min:x_max)==l);
    Edgemap_gt_crop = Edgemap_gt(y_min:y_max, x_min:x_max);
    EdgeGroupMap_gt_crop = EdgeGroupMap_gt(y_min:y_max, x_min:x_max);
    
    [match1,match2] = correspondPixels(double(Edgemap_gt_crop>0),double(Edgemap_cp_crop>0),maxDist);
        
    match_coords = match2(match2(:)>0);
    match_group_labels = EdgeGroupMap_gt_crop(match_coords);
    unique_labels = unique(match_group_labels);
    
    % measure the consistency as intersection/union
    sumP = sumP + sum(Edgemap_cp_crop(:)>0) - sum(match2(:)>0);
    cntP = cntP + sum(match2(:)>0);
    for ll = 1:length(unique_labels)
        cnt_gt_pixels = sum(EdgeGroupMap_gt_crop(:)==unique_labels(ll));
%         cnt_match_pixels = sum(match_group_labels(:)==unique_labels(ll));

        sumP = sumP + cnt_gt_pixels;
%         [match_y, match_x] = find(match2>0);
%         match_y = match_y(match_group_labels==unique_labels(ll));
%         match_x = match_x(match_group_labels==unique_labels(ll));
%         match_y = match_y + y_min -1;
%         match_x = match_x + x_min -1;
%         match_vis2(sub2ind([h,w,3], match_y, match_x, 2*ones(size(match_x))))=1;
    end
%     keyboard;
end

%%%%%%%%%% reverve process for recall


% for li = 2:length(labels_1)
%     l = labels_1(li);
%     [Y, X] = find(EdgeGroupMap_gt==l);
%     
% %     if(sum(ismember([Y, X],[86, 178],'rows')))
% %         keyboard;
% %     end
%     y_min = max(min(Y)-5, 1);
%     y_max = min(max(Y)+5, h);
%     x_min = max(min(X)-5, 1);
%     x_max = min(max(X)+5, w);
%     Edgemap_gt_crop = (EdgeGroupMap_gt(y_min:y_max, x_min:x_max)==l);
%     Edgemap_cp_crop = Edgemap_cp(y_min:y_max, x_min:x_max);
%     EdgeGroupMap_cp_crop = EdgeGroupMap_cp(y_min:y_max, x_min:x_max);
%     
%     [match1,match2] = correspondPixels(double(Edgemap_gt_crop>0),double(Edgemap_cp_crop>0),maxDist);
%     match_coords = match1(match1(:)>0);
%     match_group_labels = EdgeGroupMap_cp_crop(match_coords);
%     unique_labels = unique(match_group_labels);
%     cntR = cntR + sum(match1(:)>0);
%     for ll = 1:length(unique_labels)
% 
%         cnt_match = sum(EdgeGroupMap_cp_crop(:)==unique_labels(ll));
%         sumR = sumR + cnt_match;
% 
% %         [match_y, match_x] = find(match1>0);
% %         match_y = match_y(match_group_labels==unique_labels(ll));
% %         match_x = match_x(match_group_labels==unique_labels(ll));
% %         match_y = match_y + y_min -1;
% %         match_x = match_x + x_min -1;
% %         match_vis1(sub2ind([h,w,3], match_y, match_x, 2*ones(size(match_x))))=1;
%     end
% end

