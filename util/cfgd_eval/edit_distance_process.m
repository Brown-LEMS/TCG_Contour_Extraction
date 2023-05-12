function [gt_groups, cp_groups] = edit_distance_process(gt_index, cp_index, CFG_1, CFG_2, local_dist, edit_thresh)

gt_groups = [];
cp_groups = [];

gt_fragments = CFG_1(gt_index);
cp_fragments = CFG_2(cp_index);

gt_cp_groups = cell(2,1);

% start putting all frags as one combo
gt_all_combinations = zeros(length(gt_index), length(gt_index));
cp_all_combinations = zeros(length(cp_index), length(cp_index));

gt_all_combinations(:, 1) = gt_index';
cp_all_combinations(:, 1) = cp_index';


if (length(gt_index)>1)
    % build up the basic one fragment combo list
    gt_1 = struct('ind', gt_index', 'start_point', zeros(length(gt_index), 3), 'end_point', zeros(length(gt_index),3));

    for i = 1:length(gt_index)
        gt_1.start_point(i,:) = gt_fragments{1,i}(1, 1:3);
        gt_1.end_point(i,:) = gt_fragments{1, i}(end, 1:3);
    end

    % iteratively generate feasible combinations
    prev_c = gt_1;
    next_c = generate_combination(prev_c, gt_1, CFG_1, local_dist);
    while (~isempty(next_c.ind))
        gt_all_combinations = [gt_all_combinations; next_c.ind];
        if(size(gt_all_combinations, 1) > 500)
            gt_all_combinations = [gt_all_combinations; gt_index];       
            break;
        end
        prev_c = next_c;
        next_c = generate_combination(prev_c, gt_1, CFG_1, local_dist);
    end
        

end


if (length(cp_index)>1)
    % build up the basic one fragment combo list
    cp_1 = struct('ind', cp_index', 'start_point', zeros(length(cp_index), 3), 'end_point', zeros(length(cp_index),3));

    for i = 1:length(cp_index)
        cp_1.start_point(i,:) = cp_fragments{1,i}(1, 1:3);
        cp_1.end_point(i,:) = cp_fragments{1, i}(end, 1:3);
    end

    % iteratively generate feasible combinations
    prev_c = cp_1;
    next_c = generate_combination(prev_c, cp_1, CFG_2, local_dist);

    while (~isempty(next_c.ind))
        cp_all_combinations = [cp_all_combinations; next_c.ind];
        if(size(cp_all_combinations, 1) > 500)
            cp_all_combinations = [cp_all_combinations; cp_index];
            break;
        end 
        prev_c = next_c;
        next_c = generate_combination(prev_c, cp_1, CFG_2, local_dist);
    end



end

% ******************* only push back full combos in 1 vs multipy || multiple vs 1 case *****************************************

    if(~isempty(find(gt_all_combinations(end,:)==0)) && length(cp_index)==1)
        gt_all_combinations = [gt_all_combinations; gt_index];       
    end
    if(~isempty(find(cp_all_combinations(end,:)==0)) && length(gt_index)==1)
        cp_all_combinations = [cp_all_combinations; cp_index];
    end

gt_size = size(gt_all_combinations, 1);
cp_size = size(cp_all_combinations, 1);

% iteratively compute the cost matrix,, extract the best match for each
% iteration, remover the sub combos in all_gt_combinations and
% all_cp_combinations
while (gt_size > 0 && cp_size > 0)
    transform_cost_mat = zeros(gt_size, cp_size);
    for i = 1: gt_size
        ind_gt = gt_all_combinations(i, :);
        combo_gt = struct('ind', ind_gt(ind_gt~=0), 'start_point', [], 'end_point', []);
        for j = 1: cp_size
            ind_cp = cp_all_combinations(j, :);
            combo_cp = struct('ind', ind_cp(ind_cp~=0), 'start_point', [], 'end_point', []);
            transform_cost_mat(i,j) = compute_transform_cost(combo_gt, combo_cp, CFG_1, CFG_2, local_dist);
        end
    end
    min_cost = min(min(transform_cost_mat));
    if(min_cost > edit_thresh)
        break;
    end
        
    [min_gt_ind,min_cp_ind] = find(transform_cost_mat==min_cost);
    
    % push the min cost combos into output list
    min_gt_frags_ind = gt_all_combinations(min_gt_ind, :);
    min_cp_frags_ind = cp_all_combinations(min_cp_ind, :);
    gt_groups = [gt_groups; min_gt_frags_ind];
    cp_groups = [cp_groups; min_cp_frags_ind];
    
    % delete the sub combos of the min cost combos
    gt_2_delete = [];
    for i = 1: gt_size
        if (is_sub_combo(gt_all_combinations(i, :), min_gt_frags_ind)==1)
            gt_2_delete = [gt_2_delete, i];
        end       
    end
    
    count = 0;
    for k = 1: length(gt_2_delete)
        gt_all_combinations(gt_2_delete(k)-count, :) = [];
        count = count + 1;
    end
    
    cp_2_delete = [];
    for j = 1: cp_size
        if (is_sub_combo(cp_all_combinations(j, :), min_cp_frags_ind)==1)
            cp_2_delete = [cp_2_delete, j];
        end       
    end
    
    count = 0;
    for k = 1: length(cp_2_delete)
        cp_all_combinations(cp_2_delete(k)-count, :) = [];
        count = count + 1;
    end
    %update gt_size and cp_size
    gt_size = size(gt_all_combinations, 1);
    cp_size = size(cp_all_combinations, 1);
    
end

end