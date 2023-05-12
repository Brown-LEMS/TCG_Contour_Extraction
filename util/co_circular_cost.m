function diff_alpha = co_circular_cost(c1, c2)

%%%%% compute the probability that two curve fragment belongs to the same
%%%%% circular, order of edges: c1 ------> <---------c2

c2 = fliplr(c2');
c2 = c2'; 
%%%%% order of edges: c1 ------> -------> c2

%%%%% extremely short curve could cause error, need to be ruled out advance
%%%%% when one of the contours is super short, just return the ori diff
len = min([size(c1, 1), size(c2, 1), 5]);
dx1 = c1(end,1) - c1(end-len+1,1);
dy1 = c1(end,2) - c1(end-len+1,2);
dx2 = c2(len,1) - c2(1,1);
dy2 = c2(len,2) - c2(1,2);

local_diff_alpha = acos((dx1*dx2 + dy1*dy2) / sqrt(dx1*dx1 + dy1*dy1) / sqrt(dx2*dx2 + dy2*dy2));
if(size(c1, 1)<=5 || size(c2, 1) <=5 )
    diff_alpha = local_diff_alpha;
    return;
end

%%%%% interpolation seems necessuary as well
[x_coords1, y_coords1] = interpolate_cfrag(c1);
[x_coords2, y_coords2] = interpolate_cfrag(c2);
c1= [x_coords1(:) y_coords1(:)];
c2= [x_coords2(:) y_coords2(:)];

% another local measure on interpolated coordinates
len = min([size(c1, 1), size(c2, 1), 15]);
dx1 = c1(end,1) - c1(end-len+1,1);
dy1 = c1(end,2) - c1(end-len+1,2);
dx2 = c2(len,1) - c2(1,1);
dy2 = c2(len,2) - c2(1,2);
local_diff_alpha = acos((dx1*dx2 + dy1*dy2) / sqrt(dx1*dx1 + dy1*dy1) / sqrt(dx2*dx2 + dy2*dy2));

% actual computation
x1 = c1(1:end-2,1);
y1 = c1(1:end-2,2);
dx1 = c1(3:end, 1) - x1;
dy1 = c1(3:end, 2) - y1;
len1 = length(x1);

x2 = c2(3:end, 1);
y2 = c2(3:end, 2);
dx2 = x2-c2(1:end-2, 1);
dy2 = y2-c2(1:end-2, 2);
len2 = length(x2);    

%%%%%  Apporoach 1: full cross
%%%%% matrix computation for efficiency
x1 = repmat(x1, [1, len2]);
y1 = repmat(y1, [1, len2]);
dx1 = repmat(dx1, [1, len2]);
dy1 = repmat(dy1, [1, len2]);

x2 = repmat(x2', [len1, 1]);
y2 = repmat(y2', [len1, 1]);
dx2 = repmat(dx2', [len1, 1]);
dy2 = repmat(dy2', [len1, 1]);

dx = x2-x1;
dy = y2-y1;

% %%%%% compute the angle between circular tangent
% alpha1 = acos((dx1.*dx + dy1.*dy) ./ sqrt(dx1.*dx1 + dy1.*dy1) ./ sqrt(dx.*dx + dy.*dy));
% alpha2 = acos((dx2.*dx + dy2.*dy) ./ sqrt(dx2.*dx2 + dy2.*dy2) ./ sqrt(dx.*dx + dy.*dy));
% 
% diff_alpha = abs(alpha1 - alpha2);
% diff_alpha = mean(diff_alpha(:));

diff_alpha = abs(2*atan2(dy, dx) - atan2(dy1, dx1) - atan2(dy2, dx2));
diff_alpha = mod(diff_alpha, 2*pi);
diff_alpha(diff_alpha>pi) =  2*pi - diff_alpha(diff_alpha>pi);
% diff_alpha(diff_alpha<-pi) =  -2*pi - diff_alpha(diff_alpha<-pi);
diff_alpha = mean(diff_alpha(:));
diff_alpha = min(local_diff_alpha, diff_alpha);

%%%%%  Apporoach 2: shift one round
% len = min(len1, len2);
% x1 = x1(end-len+1:end);
% y1 = y1(end-len+1:end);
% dx1 = dx1(end-len+1:end);
% dy1 = dy1(end-len+1:end);
% 
% x2 = x2(1:len);
% y2 = y2(1:len);
% dx2 = dx2(1:len);
% dy2 = dy2(1:len);
% 
% dx = x2-x1;
% dy = y2-y1;
% 
% % %%%%% compute the angle between circular tangent
% % alpha1 = acos((dx1.*dx + dy1.*dy) ./ sqrt(dx1.*dx1 + dy1.*dy1) ./ sqrt(dx.*dx + dy.*dy));
% % alpha2 = acos((dx2.*dx + dy2.*dy) ./ sqrt(dx2.*dx2 + dy2.*dy2) ./ sqrt(dx.*dx + dy.*dy));
% % 
% % diff_alpha = abs(alpha1 - alpha2);
% % diff_alpha = mean(diff_alpha(:));
% 
% diff_alpha = abs(2*atan2(dy, dx) - atan2(dy1, dx1) - atan2(dy2, dx2));
% diff_alpha = mod(diff_alpha, 2*pi);
% diff_alpha(diff_alpha>pi) =  2*pi - diff_alpha(diff_alpha>pi);
% % diff_alpha(diff_alpha<-pi) =  -2*pi - diff_alpha(diff_alpha<-pi);
% diff_alpha = mean(diff_alpha(:));