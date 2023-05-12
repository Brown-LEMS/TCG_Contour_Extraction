function [new_row_checked, new_column_checked, new_group_gt_id, new_group_cp_id] = check_column(cost_gt2cp, cost_cp2gt, column2check, row_checked, column_checked, group_gt_id, group_cp_id, cost_thresh)

% new_column2check = column2check;
new_row_checked = row_checked;
new_column_checked = column_checked;
new_group_gt_id = group_gt_id;
new_group_cp_id = group_cp_id;

    for v = 1: length(column2check)
        col_v = cost_gt2cp(:, column2check(v));
        if(min(col_v == 1000))
            continue;
        end
        row2check = [];
        
        temp_min = 1000;
        for i = 1: length(col_v)
            if(~isempty(find(new_row_checked==i)))
                continue;
            end
            
            if(col_v(i)~=1000)
                row_i = cost_gt2cp(i,:);
                
                if(min(row_i)~= col_v(i) && col_v(i) < cost_thresh) % if this short frag matching other long frag better than cc(j), leave it for future
                    continue;
                end
                
                new_group_gt_id(1, i) = 1;
                if(col_v(i)<temp_min)
                    temp_min = col_v(i);
                    row2check = [i,row2check];
                else
                    row2check = [row2check, i];
                end
                new_row_checked = [new_row_checked, i];
            end
        end
        
        if(~isempty(row2check))
            [new_row_checked2, new_column_checked2, new_group_gt_id2, new_group_cp_id2] = check_row(cost_gt2cp, cost_cp2gt, row2check, new_row_checked, new_column_checked, new_group_gt_id, new_group_cp_id, cost_thresh);
            new_row_checked = new_row_checked2;
            new_column_checked = new_column_checked2;
            new_group_gt_id = new_group_gt_id2;
            new_group_cp_id = new_group_cp_id2;

        end
    end
%     new_column2check = [];
end