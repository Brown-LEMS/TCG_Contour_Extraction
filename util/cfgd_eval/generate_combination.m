% combine the combinations in prev_list, and one_list to generate a
% new_list of geomatrically constrained feasible combinations
function new_combo_list = generate_combination(prev_list, one_list, CFG, local_dist)
is_feasible = 1;
if(size(one_list.ind, 1)>30 || size(prev_list.ind, 1)>80)
    is_feasible = 0;
end

new_combo_list = struct('ind', [], 'start_point', [], 'end_point', []);

ind_array_length = size(one_list.ind, 1);


for i = 1: size(prev_list.ind, 1)
    index_0 = prev_list.ind(i,:);
    index_0 = index_0(index_0~=0);
    fragments_0 = CFG(index_0);
    combo_0 = struct('ind', index_0, 'start_point', prev_list.start_point(i,:), 'end_point', prev_list.end_point(i,:));
    
    for j = 1:size(one_list.ind, 1)
        index_1 = one_list.ind(j,:);
        index_1 = index_1(index_1~=0);
        fragments_1 = CFG(index_1);
        combo_1 = struct('ind', index_1, 'start_point', one_list.start_point(j,:), 'end_point', one_list.end_point(j,:));
        
        % proceed for combo_1is not include in combo_0 and not overlapping
        % part of combo_o either
        if (isempty(intersect(index_0, index_1)) && ~overlap_combos(combo_0, combo_1, CFG, local_dist))
            new_combo = struct('ind', [], 'start_point', [], 'end_point', []);
            min_dist = 1000;
            
            merge_point_0 = [];
            merge_point_1 = [];
            
            % find the closest two end points of the two input combos, and
            % assign th other two end points as new combo's end points
            dist_0 = compute_edgel_dist(combo_0.start_point, combo_1.start_point);
            if(dist_0<min_dist)
                new_combo.start_point = combo_0.end_point;
                new_combo.end_point = combo_1.end_point;
                min_dist = dist_0;
                merge_point_0 = combo_0.start_point;
                merge_point_1 = combo_1.start_point;
            end
            dist_1 = compute_edgel_dist(combo_0.start_point, combo_1.end_point);
            if(dist_1<min_dist)
                new_combo.start_point = combo_0.end_point;
                new_combo.end_point = combo_1.start_point;
                min_dist = dist_1;
                merge_point_0 = combo_0.start_point;
                merge_point_1 = combo_1.end_point;
            end
            dist_2 = compute_edgel_dist(combo_0.end_point, combo_1.start_point);
            if(dist_2<min_dist)
                new_combo.start_point = combo_0.start_point;
                new_combo.end_point = combo_1.end_point;
                min_dist = dist_2;
                merge_point_0 = combo_0.end_point;
                merge_point_1 = combo_1.start_point;
            end
            dist_3 = compute_edgel_dist(combo_0.end_point, combo_1.end_point);
            if(dist_3<min_dist)
                new_combo.start_point = combo_0.start_point;
                new_combo.end_point = combo_1.start_point;
                min_dist = dist_3;
                merge_point_0 = combo_0.end_point;
                merge_point_1 = combo_1.end_point;
            end
            
            new_combo.ind = [new_combo.ind, combo_0.ind];
            new_combo.ind = [new_combo.ind, combo_1.ind];
            
            % check if the new formed combo already exist in the combo
            is_exist = 0;
            if(~isempty(new_combo_list.ind))
                for c = 1:size(new_combo_list.ind, 1)
                    temp_ind = new_combo_list.ind(c,:);
                    temp_ind = temp_ind(temp_ind ~= 0);
                    if (sort(new_combo.ind) == sort(temp_ind))
                        is_exist = 1;
                        break;
                    end
                end
            end
            
            if(is_exist)
                continue;
            end
            
            % check if new combo can form a closee contour, or end points
            % are close enough, just push new_combo into output list
            if((dist_1<local_dist && dist_2 < local_dist) || (dist_0<local_dist && dist_3 < local_dist) || min_dist < 2*local_dist)
                new_ind = zeros(1, ind_array_length);
                new_ind(1, 1:length(new_combo.ind)) = new_combo.ind;
                new_combo_list.ind = [new_combo_list.ind; new_ind];
                new_combo_list.start_point = [new_combo_list.start_point; new_combo.start_point];
                new_combo_list.end_point = [new_combo_list.end_point; new_combo.end_point];
            elseif(min_dist >=2*local_dist && min_dist < 10 && is_feasible == 1)
                d_x_0=0; d_y_0=0; d_x_1=0; d_y_1=0;
                d_x_merge = merge_point_1(1,1) - merge_point_0(1,1);
                d_y_merge = merge_point_1(1,2) - merge_point_0(1,2);
                
                for k = 1:size(fragments_0, 2)
                    fragment_0 = fragments_0{1,k};
                    if (merge_point_0 == fragment_0(1, 1:3))
                        if( size(fragment_0, 1) >=5)
                            d_x_0 = merge_point_0(1,1) - fragment_0(5,1);
                            d_y_0 = merge_point_0(1,2) - fragment_0(5,2);
                        else
                            d_x_0 = merge_point_0(1,1) - fragment_0(end,1);
                            d_y_0 = merge_point_0(1,2) - fragment_0(end,2);
                        end
                    elseif (merge_point_0 == fragment_0(end, 1:3))
                        if( size(fragment_0, 1) >=5)
                            d_x_0 = merge_point_0(1,1) - fragment_0(end-4,1);
                            d_y_0 = merge_point_0(1,2) - fragment_0(end-4,2);
                        else
                            d_x_0 = merge_point_0(1,1) - fragment_0(1,1);
                            d_y_0 = merge_point_0(1,2) - fragment_0(1,2);
                        end                        
                    end
                end
                
                for k = 1:size(fragments_1, 2)
                    fragment_1 = fragments_1{1,k};
                    if (merge_point_1 == fragment_1(1, 1:3))
                        if( size(fragment_1, 1) >=5)
                            d_x_1 = merge_point_1(1,1) - fragment_1(5,1);
                            d_y_1 = merge_point_1(1,2) - fragment_1(5,2);
                        else
                            d_x_1 = merge_point_1(1,1) - fragment_1(end,1);
                            d_y_1 = merge_point_1(1,2) - fragment_1(end,2);
                        end
                    elseif (merge_point_1 == fragment_1(end, 1:3))
                        if( size(fragment_1, 1) >=5)
                            d_x_1 = merge_point_1(1,1) - fragment_1(end-4,1);
                            d_y_1 = merge_point_1(1,2) - fragment_1(end-4,2);
                        else
                            d_x_1 = merge_point_1(1,1) - fragment_1(1,1);
                            d_y_1 = merge_point_1(1,2) - fragment_1(1,2);
                        end                        
                    end
                end                
                
                % check if there is other contours lie between them, if not
                %, approve the combination
                is_approved = 1;
                if((d_x_0*d_x_merge + d_y_0*d_y_merge)/sqrt(d_x_0^2 + d_y_0^2)/sqrt(d_x_merge^2+d_y_merge^2) < 0.866 &&...
                        (d_x_merge*d_x_1 + d_y_merge*d_y_1)/sqrt(d_x_merge^2 + d_y_merge^2)/sqrt(d_x_1^2 + d_y_1^2) < 0.866)
                    is_approved = 0;
                end
                if((d_x_0*d_x_1 + d_y_0*d_y_1)/sqrt(d_x_0^2+d_y_0^2)/sqrt(d_x_1^2+d_y_1^2)>0.866)
                    is_approved = 1;
                end
                
                for k = 1:size(one_list.ind, 1)
                    index_k = one_list.ind(k,1);                    
                    if(is_approved == 0)
                        break;
                    end
                    % exclued the situation combo_k is combo_1, or it is
                    % included in combo_0
                    if(k==j || ~isempty(intersect(index_0, index_k)))
                        continue;
                    end
                    
                    % if there are other contour lying between the merging
                    % points
                    Edges = CFG(index_k);
                    if(is_between(Edges', merge_point_0', merge_point_1', min_dist)>0)
                        is_approved = 0;
                        break;
                    end
                    
                end
                
                if(is_approved ==1)
                    new_ind = zeros(1, ind_array_length);
                    new_ind(1, 1:length(new_combo.ind)) = new_combo.ind;
                    new_combo_list.ind = [new_combo_list.ind; new_ind];
                    new_combo_list.start_point = [new_combo_list.start_point; new_combo.start_point];
                    new_combo_list.end_point = [new_combo_list.end_point; new_combo.end_point];                    
                end
                     
            end
            
                
        end
        
    end
end
end