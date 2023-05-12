% compute the total length of contours in combo_1
function l = combo_length(combo_1, CFG_1)

l=0;
ind_1 = combo_1.ind;
frags_1 = CFG_1(ind_1);

for i = 1:size(frags_1, 2)
    l = l + contour_length_mex(frags_1{1,i}');
end

end