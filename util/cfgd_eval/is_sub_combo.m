% this function check if combo_1 and combo_2 share any curve fragment
% but combo_1 and combo_2 must only contain fragments index from the same CFG
function b = is_sub_combo(combo_1_ind, combo_2_ind)
    new_combo_1_ind = combo_1_ind(combo_1_ind~=0);
    new_combo_2_ind = combo_2_ind(combo_2_ind~=0);
    common = intersect(new_combo_1_ind, new_combo_2_ind);
    if(isempty(common))
        b = 0;
    else
        b = 1;
    end
end