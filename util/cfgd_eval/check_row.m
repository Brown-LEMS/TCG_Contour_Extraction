function [new_row_checked, new_column_checked, new_group_gt_id, new_group_cp_id] = check_row(cost_gt2cp, cost_cp2gt, row2check, row_checked, column_checked, group_gt_id, group_cp_id, cost_thresh)
new_row_checked = row_checked;
new_column_checked = column_checked;
new_group_gt_id = group_gt_id;
new_group_cp_id = group_cp_id;

for v = 1:length(row2check)
    row_v = cost_cp2gt(row2check(v), :);
    if(min(row_v)==1000)
        continue;
    end
    column2check = [];
    temp_min = 1000;
    for j = 1: length(row_v)
        if(~isempty(find(new_column_checked==j)))
            continue;
        end
        if(row_v(j)~=1000)
            col_j = cost_cp2gt(:,j);
            if(min(col_j)~=row_v(j) && row_v(j)<cost_thresh)
                continue;
            end
            new_group_cp_id(1,j)=1;
            if (row_v(j)<temp_min)
                temp_min = row_v(j);
                column2check = [j, column2check];
            else
                column2check = [column2check, j];
            end
            new_column_checked = [new_column_checked, j];
        end
    end
    if(~isempty(column2check))
        [new_row_checked2, new_column_checked2, new_group_gt_id2, new_group_cp_id2] = check_column(cost_gt2cp, cost_cp2gt, column2check, new_row_checked, new_column_checked, new_group_gt_id, new_group_cp_id, cost_thresh);        
        new_row_checked = new_row_checked2;
        new_column_checked = new_column_checked2;
        new_group_gt_id = new_group_gt_id2;
        new_group_cp_id = new_group_cp_id2;
    end
end
end