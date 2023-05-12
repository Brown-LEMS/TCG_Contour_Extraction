function [ori_diff_vec, Dtheta] = filter_co_circular_along_cfrag(cfrag, nbr_range_th)
%%%%%%%% ori_diff_vec is in the range [0, pi];

% nbr_range_th = min(nbr_range_th, 5); % does not make sense for too long range

len = size(cfrag,1);
ori_diff_vec = zeros(len,1);

V = cfrag(:, 1:2);
V = V +1; % change CXX to the matlab coordinates

X = V(:,1);
% X = smooth(X);
Y = V(:,2);
% Y = smooth(Y);
DX = X(3:end) - X(1:end-2);
DY = Y(3:end) - Y(1:end-2);
DX = [DX(1); DX; DX(end)];
DY = [DY(1); DY; DY(end)];

theta_vec = atan2(DY, DX);
Dtheta = abs(diff(theta_vec));
Dtheta = [0; Dtheta];
Dtheta(Dtheta>pi) = 2*pi - Dtheta(Dtheta>pi);
%% compute prob along each curve
%%%%%%%%%%%%%%%%%%  TODO: can be converted into matrix
%%%%%%%%%%%%%%%%%%  computation to get rid of loop
for i = nbr_range_th +1: len-nbr_range_th
    
    %%%%%%%%%%%% compute ori diff at the connecting points

    dx1 = DX(i-nbr_range_th:i-1);
    dy1 = DY(i-nbr_range_th:i-1);
    
    dx2 = DX(i+1:i+nbr_range_th);
    dy2 = DY(i+1:i+nbr_range_th);
    
    dx = X(i+1:i+nbr_range_th) - X(i-nbr_range_th:i-1);
    dy = Y(i+1:i+nbr_range_th) - Y(i-nbr_range_th:i-1);
    
%     alpha1 = acos((dx1.*dx + dy1.*dy) ./ sqrt(dx1.*dx1 + dy1.*dy1) ./ sqrt(dx.*dx + dy.*dy));
%     alpha2 = acos((dx2.*dx + dy2.*dy) ./ sqrt(dx2.*dx2 + dy2.*dy2) ./ sqrt(dx.*dx + dy.*dy));
% 
%     diff_alpha = abs(alpha1 - alpha2);
%     diff_alpha = mean(diff_alpha(:));
    
    diff_alpha = abs(2*atan2(dy, dx) - atan2(dy1, dx1) - atan2(dy2, dx2));
    diff_alpha = mod(diff_alpha, 2*pi);
    diff_alpha(diff_alpha>pi) =  2*pi - diff_alpha(diff_alpha>pi);
%     diff_alpha(diff_alpha<-pi) =  -2*pi - diff_alpha(diff_alpha<-pi);
    diff_alpha = mean(diff_alpha(:));
    
    ori_diff_vec(i) = diff_alpha;
    
end

end