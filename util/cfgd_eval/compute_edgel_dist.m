function d = compute_edgel_dist(e1, e2)
% compute the localization distance between two edges
p1 = e1(1, 1:2);
p2 = e2(1, 1:2);
d = sqrt(sum((p1-p2).^2));


end