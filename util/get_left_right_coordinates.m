function [left_x, left_y, right_x, right_y] = get_left_right_coordinates(V, N, local_dist, w, h)

% x = round(V(:,1));
% y = round(V(:,2));

% N's direction matters
left_x = round(V(:,1) - local_dist*N(:,1));
left_y = round(V(:,2) - local_dist*N(:,2));
right_x = round(V(:,1) + local_dist*N(:,1));
right_y = round(V(:,2) + local_dist*N(:,2));

% deal with the local exceeding the border
left_x = max(left_x, ones(size(left_x)));
left_y = max(left_y, ones(size(left_y)));
left_x = min(left_x, w*ones(size(left_x)));
left_y = min(left_y, h*ones(size(left_y)));
right_x = max(right_x, ones(size(right_x)));
right_y = max(right_y, ones(size(right_y)));
right_x = min(right_x, w*ones(size(right_x)));
right_y = min(right_y, h*ones(size(right_y)));     

end