function bry_map_trans = transform_bry_map_given_GT_disparity( bry_map, GT_disparity_map, is_left)

%% translate edges x coordinates
[h,w] = size(GT_disparity_map);
% x_coords = round(edges(:,1))+1;
% y_coords = round(edges(:,2))+1;
[Y, X] = find(bry_map>0);
edges = [X, Y]; % in matlab coordinates

d_vec = [];
for dx = [-1 0 1]
    x_coords = round(edges(:,1)) + dx;
    x_coords =  max(x_coords, 1);
    x_coords = min(x_coords, w);
    for dy = [-1 0 1]
        y_coords = round(edges(:,2)) + dy;
        y_coords =  max(y_coords, 1);
        y_coords = min(y_coords, h);

        d_vec = [d_vec GT_disparity_map(sub2ind([h,w], y_coords, x_coords))];
    end
end

% ind = sub2ind([h,w], y_coords, x_coords);
% d = GT_disparity_map(ind);
d = max(d_vec, [], 2);

edges_trans =  edges;
if(is_left)
    edges_trans(:,1) = edges_trans(:,1) - d;
else
    edges_trans(:,1) = edges_trans(:,1) + d;
end

edges_trans(:,1) =  max(1, round(edges_trans(:,1)));
edges_trans(:,1) =  min(w, round(edges_trans(:,1)));


bry_map_trans = zeros(size(bry_map));
bry_map_trans(sub2ind([h,w], Y, edges_trans(:,1))) = bry_map(sub2ind([h, w], Y, X));







end