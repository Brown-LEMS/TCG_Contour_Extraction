function [TP_GT_L, TP_CP_L, GT_L, CP_L, match_GT_frags, match_CP_frags, miss_GT_frags, extra_CP_frags, match_select_frags] = Compare_Curve_Fragment_Maps_v3_new_score(CEM_1, CEM_2, prune_len_thresh, display_image, bgroup_only, cost_thresh, edit_thresh, local_dist)
% outputs: matched length of GT, matched length of CP, total length of
% GT, total lenght of CP
% inputs: groud_truth_path, computed_map_path, prune length threshold for
% computed fragments

% Parameters
if nargin<4
    display_image = 0; % choose to display the visualiztion of the matching of curve fragments
end
if nargin<5
    bgroup_only = 0; % Stop After Grouping, no edit distance process within each group
end
if nargin<6
    cost_thresh = 3; % For pruning out unrelated curve fragments and grouping process 
end
if nargin<7
    edit_thresh = 2; % Edit Distance threshold pruning out unqualified matching candidats
end
if nargin<8
    local_dist = 3; % Allowed localization Error
end

theta_diff = deg2rad(60);
maxDist = 0.01;

match_GT_frags = [];
match_CP_frags = [];
miss_GT_frags = [];
extra_CP_frags = [];
match_select_frags = [];

if(nargin < 3)
    prune_len_thresh = 0; % the pruned length_for computed maps, no prune out for 
end
fprintf('Prune Length Thresh: %d\n', prune_len_thresh);

% % Load in Curve Fragment Files
% fprintf('Loading contour maps...\n');
% fprintf([GT_map_path '\n']);
% fprintf([CP_map_path '\n']);
% % CEM_1=load_contours('../2018_s1.cem');
% % CEM_2=load_contours('../2018.cem');
% if(~isempty(strfind(GT_map_path, '.cem')))
%     CEM_1=load_contours(GT_map_path);
% elseif(~isempty(strfind(GT_map_path, '.mat')))
%     load(GT_map_path);
%     CEM_1 = CEM;
% else
%     fprintf('First input file has wrong format');
%     return;
% end
% 
% 
% if(~isempty(strfind(CP_map_path, '.cem')))
%     CEM_2=load_contours(CP_map_path);
% elseif(~isempty(strfind(CP_map_path, '.mat')))
%     load(CP_map_path);
%     CEM_2 = CEM;
% else
%     fprintf('Second input file has wrong format');
%     return;
% end

%%
% Curve Fragment Graph (CFG)
CFG_1 = CEM_1{2,1};
CFG_2 = CEM_2{2,1};
% gt_total_length = 0;
% cp_total_length = 0;

% % display the original two curve fragment sets
% if(display_image>0)
%     imgBW = zeros(CEM_1{1,1}(1,2), CEM_1{1,1}(1,1),3);
%     h = figure();
%     imshow(imgBW, 'border', 'tight');hold on;
%     draw_red_contours(CFG_1);
%     draw_green_contours(CFG_2);
%     hold off;
% end

% build up the cost_mat to remove unrelated 
cost_mat = zeros(size(CFG_1,2), size(CFG_2,2));
% Compute the total contour length of each CFG
% at the same time, prune out some computed fragments under some length
% threshold by setting the corresponding column in cost_mat as 1000
fprintf('compute total length of contour in each map...\n');


% [binary_gt, orient_gt] = convert_contours(CEM_1,3);
% [edgemap, thetamap] = convert_contours(CEM_2,3);
% 
% gt_total_length = sum(binary_gt(:));
% cp_total_length = sum(edgemap(:));
        
gt_total_length = 0;
cp_total_length = 0;
% for l = 1:size(CFG_1, 2)
%     c_length = contour_length_mex(CFG_1{1,l}');
%     gt_total_length = gt_total_length + c_length;
% end
if(prune_len_thresh >0)
    for l = 1:size(CFG_2, 2)
        c_length = contour_length_mex(CFG_2{1,l}');
        if(c_length < prune_len_thresh)
           cost_mat(:,l) = cost_mat(:,l)+1000; 
%            continue;
        end       
    %     cp_total_length = cp_total_length + c_length;
    end   
end

fprintf('compute cost matrix...\n');
for i = 1:size(CFG_1,2)
    for j = 1:size(CFG_2,2)
%         disp(sprintf('compute cost matrix %d, %d', i, j));
        if (cost_mat(i,j)==1000)
            continue;
        end
        cost_mat(i,j) = min(compute_contours_cost_mex(CFG_1{1,i}', CFG_2{1,j}', cost_thresh),...
            compute_contours_cost_mex(CFG_2{1,j}', CFG_1{1,i}', cost_thresh));
    end
end

% get the indexes of related contours
[r, c] = find(cost_mat<1000);
r = unique(r);
c = unique(c);


% for special case
if(size(r, 1)==1 && size(r, 2)>1)
    r = r';
end


if(size(c, 1)==1 && size(c, 2)>1)
    c = c';
end

% push the miss and extra in lists
tp_rm_ind = setdiff(1:size(cost_mat, 1), r);
cp_rm_ind = setdiff(1:size(cost_mat, 2), c);

miss_GT_frags = [miss_GT_frags CFG_1(1, tp_rm_ind)];
extra_CP_frags = [extra_CP_frags CFG_2(1,cp_rm_ind)];

fprintf('number of refined contours in CFG_1: %d, CFG_2: %d \n', length(r), length(c));
New_CFG_1 = cell(1,length(r));
New_CFG_2 = cell(1,length(c));
for i = 1:length(r)
    New_CFG_1{1,i} = CFG_1{1, r(i,1)};
end

for j = 1:length(c)
    New_CFG_2{1,j} = CFG_2{1, c(j,1)};
end

% CEM_1{2,1} = New_CFG_1;
% CEM_2{2,1} = New_CFG_2;

% sort the two refined curve fragment maps by the length of curve fragment
% in descent order
[Cnrows,Cncols] = cellfun(@size, New_CFG_1);
[vals,I] = sort(Cnrows, 'descend');
New_CFG_1 = New_CFG_1(I);

[Cnrows,Cncols] = cellfun(@size, New_CFG_2);
[vals,I] = sort(Cnrows, 'descend');
New_CFG_2 = New_CFG_2(I);

% % display the refine two curve fragment sets
if(display_image>0)
    imgBW = zeros(CEM_1{1,1}(1,2), CEM_1{1,1}(1,1),3);
    h = figure(1);
    imshow(imgBW, 'border', 'tight');hold on;
    draw_red_contours(New_CFG_1);
    draw_green_contours(New_CFG_2);
    hold off;
end

%% build up the cost_mat_gt_cp to remove to hold cost from Ground Truth to Computed fragments 
row = size(New_CFG_1,2);
col = size(New_CFG_2,2);

cost_gt2cp = zeros(row, col);
cost_cp2gt = zeros(row, col);

fprintf('compute GT to CP cost matrix...\n');
for i = 1:row
    for j = 1:col
%         disp(sprintf('compute cost matrix %d, %d', i, j));
        cost_gt2cp(i,j) = compute_contours_cost_mex(New_CFG_1{1,i}', New_CFG_2{1,j}', cost_thresh);
    end
end

fprintf('compute CP to GT cost matrix...\n');
for i = 1:row
    for j = 1:col
%         disp(sprintf('compute cost matrix %d, %d', i, j));
        cost_cp2gt(i,j) = compute_contours_cost_mex(New_CFG_2{1,j}', New_CFG_1{1,i}', cost_thresh);
    end
end

%% Group Raletive Curve Fragments between GT map and CP map
fprintf('Grouping and Edit Distance Processing...\n');
gt_match_length = 0;
cp_match_length = 0;
gt_index_group = [];
cp_index_group = [];

gt_matched_ind = [];
cp_matched_ind = [];

row_checked = [];
column_checked = [];
set(0,'RecursionLimit',2000);
% group_gt_id = zeros(1, row);
% group_cp_id = zeros(1, col);

if(display_image > 0)
    imgBW = zeros(CEM_1{1,1}(1,2), CEM_1{1,1}(1,1),3);
    h = figure(2);
    imshow(imgBW, 'border', 'tight');hold on;
end

% start checking columns first
for j = 1:col
    if(~isempty(find(column_checked==j, 1)))
        continue;
    end
    group_gt_id = zeros(1, row);
    group_cp_id = zeros(1, col);
    column2check = [];
    group_cp_id(1,j) = 1;
    column_checked = [column_checked, j];
    column2check = [column2check, j];
    
    col_v = cost_gt2cp(:, j);
    
    if(min(col_v)==1000)
        column_checked = column_checked(find(column_checked~=j));
        group_cp_id(1,j) = 0;
        continue;
    end
    

    [new_row_checked, new_column_checked, new_group_gt_id, new_group_cp_id] = check_column(cost_gt2cp, cost_cp2gt, column2check, row_checked, column_checked, group_gt_id, group_cp_id, cost_thresh);
    row_checked = new_row_checked;
    column_checked = new_column_checked;
    group_gt_id = new_group_gt_id;
    group_cp_id = new_group_cp_id;
    
    if(max(group_gt_id)==0 || max(group_cp_id)==0)
        continue;
    else
        %   display the grouped fragments
        [ones, gt_index] = find(new_group_gt_id==1);
        [ones, cp_index] = find(new_group_cp_id==1);
        if(bgroup_only>0)
%             group_gt_fragments = New_CFG_1(gt_index);
%             group_cp_fragments = New_CFG_2(cp_index);
%             gt_matched_ind = [gt_matched_ind gt_index];
%             cp_matched_ind = [cp_matched_ind cp_index];
%             match_GT_frags = [match_GT_frags group_gt_fragments];
%             match_CP_frags = [match_CP_frags group_cp_fragments];
% 
%             for l = 1:size(group_gt_fragments, 2)
%                 gt_match_length = gt_match_length + contour_length_mex(group_gt_fragments{1,l}');
%             end            
%             for l = 1:size(group_cp_fragments, 2)
%                 cp_match_length = cp_match_length + contour_length_mex(group_cp_fragments{1,l}');
%             end         
%             if(display_image>0)
%                 draw_contour_pairs(group_gt_fragments, group_cp_fragments); 
%             end
        else
            [gt_groups_ind, cp_groups_ind] = edit_distance_process(gt_index, cp_index, New_CFG_1, New_CFG_2, local_dist, edit_thresh);
            
            % this part count for the fragments pruned during edit distance
            % process
            A = unique(gt_groups_ind);
            B = unique(cp_groups_ind);
            if(~isempty(find(A==0)))
                A = A(2:end);
            end
            if(~isempty(find(B==0)))
                B = B(2:end);
            end            
            
            gt_matched_ind = [gt_matched_ind; A(:)];
            cp_matched_ind = [cp_matched_ind; B(:)];
%             miss_ind = setdiff(A, gt_index);
%             if(~isempty(miss_ind))
%                 miss_GT_frags = [miss_GT_frags New_CFG_1(miss_ind)];                
%             end
%             
%             extra_ind = setdiff(B, cp_index);
%             if(~isempty(extra_ind))
%                 extra_CP_frags = [extra_CP_frags New_CFG_2(extra_ind)];                
%             end
            
            
            for k = 1: size(gt_groups_ind, 1)
                gt_ind = gt_groups_ind(k,:);
                gt_ind = gt_ind(gt_ind~=0);
                cp_ind = cp_groups_ind(k,:);
                cp_ind = cp_ind(cp_ind~=0);
                group_gt_fragments = New_CFG_1(gt_ind);
                                
%                 for l = 1:size(group_gt_fragments, 2)
%                     gt_match_length = gt_match_length + contour_length_mex(group_gt_fragments{1,l}');
%                 end
                group_cp_fragments = New_CFG_2(cp_ind);
%                 for l = 1:size(group_cp_fragments, 2)
%                     cp_match_length = cp_match_length + contour_length_mex(group_cp_fragments{1,l}');
%                 end  
                match_GT_frags = [match_GT_frags group_gt_fragments];
                match_CP_frags = [match_CP_frags group_cp_fragments];
                % compute new version gt_match_length and cp_match_length,
                % only consider the matched poition rather than the whold
                % lengt of matched curve fragments
                CEM_1_1 = CEM_1;
                CEM_1_1{2} = group_gt_fragments;
                
                CEM_2_1 = CEM_2;
                CEM_2_1{2} = group_cp_fragments;
                
                [binary_1, orient_1] = convert_contours(CEM_1_1,3);
                [binary_2, orient_2] = convert_contours(CEM_2_1,3);
 
                        % edge evaluation with orientation
                [match1,match2] = correspondPixels_theta(double(binary_2),double(binary_1),orient_2,orient_1,theta_diff,maxDist);
                
%                 keyboard;
                accP = zeros(size(binary_1));
                accR = zeros(size(binary_2));

                % accumulate machine matches
                accP = accP | match1;
                % accumulate ground truth matches
                accR = accR | match2;

                % save the counts of matched pixels
                cntR = sum(accR(:));
                cntP = sum(accP(:)); 
                
                gt_match_length = gt_match_length + cntR;
                cp_match_length = cp_match_length + cntP;
                
                gt_total_length = gt_total_length + sum(binary_1(:));
                cp_total_length = cp_total_length + sum(binary_2(:));
        
                % add in match_select_frags
                group_length_1 = 0;
                group_length_2 = 0;
                
                for f1 = 1:size(group_gt_fragments, 2)
                    group_length_1 = group_length_1 + contour_length_mex(group_gt_fragments{1,f1}');
                end
                for f2 = 1:size(group_cp_fragments, 2)
                    group_length_2 = group_length_2 + contour_length_mex(group_cp_fragments{1,f2}');
                end                
                
                if(size(group_gt_fragments, 2) < size(group_cp_fragments, 2))
                    match_select_frags = [match_select_frags group_gt_fragments];
                elseif(size(group_gt_fragments, 2) > size(group_cp_fragments, 2))
                    match_select_frags = [match_select_frags group_cp_fragments];
                else
                    if(group_length_1 > group_length_2)
                        match_select_frags = [match_select_frags group_gt_fragments];
                    else
                        match_select_frags = [match_select_frags group_cp_fragments];
                    end
                end
                
                if(display_image>0)
                    draw_contour_pairs(group_gt_fragments, group_cp_fragments);
                end
            end
        end
    end
    

end

% start checking row for the remaining
for i = 1:row
    if(~isempty(find(row_checked==i, 1)))
        continue;
    end
    group_gt_id = zeros(1, row);
    group_cp_id = zeros(1, col);
    row2check = [];
    group_gt_id(1,i) = 1;
    row_checked = [row_checked, i];
    row2check = [row2check, i];
    
    row_v = cost_cp2gt(i, :);
    
    if (min(row_v)==1000)
        row_checked = row_checked(find(row_checked~=i));
        group_gt_id(1,i) = 0;
        continue;
    end
    
    [new_row_checked, new_column_checked, new_group_gt_id, new_group_cp_id] = check_row(cost_gt2cp, cost_cp2gt, row2check, row_checked, column_checked, group_gt_id, group_cp_id, cost_thresh);
    row_checked = new_row_checked;
    column_checked = new_column_checked;
    group_gt_id = new_group_gt_id;
    group_cp_id = new_group_cp_id;
    
    if(max(group_gt_id)==0 || max(group_cp_id)==0)
        continue;
    else
        %   display the grouped fragments
        [ones, gt_index] = find(new_group_gt_id==1);
        [ones, cp_index] = find(new_group_cp_id==1);
        if(bgroup_only>0)
%             group_gt_fragments = New_CFG_1(gt_index);
%             group_cp_fragments = New_CFG_2(cp_index);
%             gt_matched_ind = [gt_matched_ind gt_index];
%             cp_matched_ind = [cp_matched_ind cp_index];
%             match_GT_frags = [match_GT_frags group_gt_fragments];
%             match_CP_frags = [match_CP_frags group_cp_fragments];
%             for l = 1:size(group_gt_fragments, 2)
%                 gt_match_length = gt_match_length + contour_length_mex(group_gt_fragments{1,l}');
%             end            
%             for l = 1:size(group_cp_fragments, 2)
%                 cp_match_length = cp_match_length + contour_length_mex(group_cp_fragments{1,l}');
%             end            
%             if(display_image>0)
%                 draw_contour_pairs(group_gt_fragments, group_cp_fragments); 
%             end
        else
            [gt_groups_ind, cp_groups_ind] = edit_distance_process(gt_index, cp_index, New_CFG_1, New_CFG_2, local_dist, edit_thresh);

            
            % this part count for the fragments pruned during edit distance
            % process
            A = unique(gt_groups_ind);
            B = unique(cp_groups_ind);
            if(~isempty(find(A==0)))
                A = A(2:end);
            end
            if(~isempty(find(B==0)))
                B = B(2:end);
            end            
            
            gt_matched_ind = [gt_matched_ind; A(:)];
            cp_matched_ind = [cp_matched_ind; B(:)];
            
%             miss_ind = setdiff(A, gt_index);
%             if(~isempty(miss_ind))
%                 miss_GT_frags = [miss_GT_frags New_CFG_1(miss_ind)];                
%             end
%             
%             extra_ind = setdiff(B, cp_index);
%             if(~isempty(extra_ind))
%                 extra_CP_frags = [extra_CP_frags New_CFG_2(extra_ind)];                
%             end
%             
            
            for k = 1: size(gt_groups_ind, 1)
                gt_ind = gt_groups_ind(k,:);
                gt_ind = gt_ind(gt_ind~=0);
                cp_ind = cp_groups_ind(k,:);
                cp_ind = cp_ind(cp_ind~=0);
                group_gt_fragments = New_CFG_1(gt_ind);
%                 for l = 1:size(group_gt_fragments, 2)
%                     gt_match_length = gt_match_length + contour_length_mex(group_gt_fragments{1,l}');
%                 end
                group_cp_fragments = New_CFG_2(cp_ind);
%                 for l = 1:size(group_cp_fragments, 2)
%                     cp_match_length = cp_match_length + contour_length_mex(group_cp_fragments{1,l}');
%                 end
                match_GT_frags = [match_GT_frags group_gt_fragments];
                match_CP_frags = [match_CP_frags group_cp_fragments];
                
                % compute new version gt_match_length and cp_match_length,
                % only consider the matched poition rather than the whold
                % lengt of matched curve fragments
                CEM_1_1 = CEM_1;
                CEM_1_1{2} = group_gt_fragments;
                
                CEM_2_1 = CEM_2;
                CEM_2_1{2} = group_cp_fragments;
                
                [binary_1, orient_1] = convert_contours(CEM_1_1,3);
                [binary_2, orient_2] = convert_contours(CEM_2_1,3);
 
                        % edge evaluation with orientation
                [match1,match2] = correspondPixels_theta(double(binary_2),double(binary_1),orient_2,orient_1,theta_diff,maxDist);
                
%                 keyboard;
                accP = zeros(size(binary_1));
                accR = zeros(size(binary_2));

                % accumulate machine matches
                accP = accP | match1;
                % accumulate ground truth matches
                accR = accR | match2;

                % save the counts of matched pixels
                cntR = sum(accR(:));
                cntP = sum(accP(:)); 
                
                gt_match_length = gt_match_length + cntR;
                cp_match_length = cp_match_length + cntP;
                gt_total_length = gt_total_length + sum(binary_1(:));
                cp_total_length = cp_total_length + sum(binary_2(:));
                
                % add in match_select_frags
                group_length_1 = 0;
                group_length_2 = 0;
                
                for f1 = 1:size(group_gt_fragments, 2)
                    group_length_1 = group_length_1 + contour_length_mex(group_gt_fragments{1,f1}');
                end
                for f2 = 1:size(group_cp_fragments, 2)
                    group_length_2 = group_length_2 + contour_length_mex(group_cp_fragments{1,f2}');
                end                
                
                if(size(group_gt_fragments, 2) < size(group_cp_fragments, 2))
                    match_select_frags = [match_select_frags group_gt_fragments];
                elseif(size(group_gt_fragments, 2) > size(group_cp_fragments, 2))
                    match_select_frags = [match_select_frags group_cp_fragments];
                else
                    if(group_length_1 > group_length_2)
                        match_select_frags = [match_select_frags group_gt_fragments];
                    else
                        match_select_frags = [match_select_frags group_cp_fragments];
                    end
                end
                
                % display groups
                if(display_image>0)
                    draw_contour_pairs(group_gt_fragments, group_cp_fragments);
                end
            end
        end
    end
        
    
end

% Recall = gt_match_length/gt_total_length
% Precision = cp_match_length/cp_total_length
% F = 2*Recall*Precision/(Recall+Precision)

% push the miss and extra in lists
tp_rm_ind2 = setdiff(1:row, gt_matched_ind);
cp_rm_ind2 = setdiff(1:col, cp_matched_ind);

miss_GT_frags = [miss_GT_frags New_CFG_1(1, tp_rm_ind2)];
extra_CP_frags = [extra_CP_frags New_CFG_2(1,cp_rm_ind2)];

if (display_image>0)
    imgBW = zeros(CEM_1{1,1}(1,2), CEM_1{1,1}(1,1),3);
    h = figure(3);
    imshow(imgBW, 'border', 'tight');hold on;
    draw_red_contours(miss_GT_frags);
    draw_green_contours(extra_CP_frags);
    hold off;
end

% Recall = gt_match_length/gt_total_length
% Precision = cp_match_length/cp_total_length
% F = 2*Recall*Precision/(Recall+Precision)



TP_GT_L = gt_match_length;
TP_CP_L = cp_match_length;
GT_L = gt_total_length; 
CP_L = cp_total_length;

recall = TP_GT_L/GT_L;
precision = TP_CP_L/CP_L;
F = 2*recall*precision/(recall+precision);
fprintf(['\n recall: ' num2str(recall) '\n precision: ' num2str(precision) '\n F-measure: ' num2str(F) '\n Evaluation Done... \n\n']);

end
