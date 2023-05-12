function d_theta = compute_ori_diff(e1,e2)
% compute the orientation differece between two edges
dir_1 = e1(1,3);
dir_2 = e2(1,3);

if(dir_1 >= pi)
    dir_1 = dir_1 - pi;
end
if(dir_2 >= pi)
    dir_2 = dir_2 - pi;
end

d_theta = abs(dir_2-dir_1);
if(d_theta >= pi*0.5)
    d_theta = pi - d_theta;
end

end