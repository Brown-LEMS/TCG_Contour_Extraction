function [bg_grad, abs_k, edge_sparsity, wigg, c_len, high_brg, low_brg, bg_grad_std, high_brg_std, low_brg_std] = curve_fragment_cues_gray(cfrag, gray_img, edge_map)
% cfrags is a list of edges which form a single curve fragment, edges in
% CXX coordinates,
% hsv_map is the hsv space matrix from input image
% edge_map is the binary map of detected edges before edge linking

local_dist = 3; % for general image
nbr_width = 3; % region for edge sparsity
gray_img = double(gray_img);

[h w ~] = size(gray_img);
len = size(cfrag,1);

V = cfrag(:, 1:2);
V = V +1; % change CXX to the matlab coordinates
K = LineCurvature2D(V);


if( sum(isnan(K))>0)
%     keyboard;
    K(find(isnan(K))) =0;
end

% wiggliness is the time curvature change sign

sign_vec = K(1:end-1).*K(2:end);
wigg = sum(find(sign_vec<0))/len;

N = LineNormals2D(V); % N's direction seems stable and values are normalized to unit already

x = round(V(:,1)); x = min(w, x); x = max(1,x);
y = round(V(:,2)); y = min(h, y); y = max(1,y);

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

% compute the attribute along left side and right side of the curve 

left_brg = diag(gray_img(left_y, left_x));
right_brg = diag(gray_img(right_y, right_x));
K = abs(K); % only use the abs, ignore the direction

% integral of "attribute" gradient along the curve sum of abs can 
% remove the influence of changing direction in the Norm computation
bg_grad = mean(abs(left_brg-right_brg));
bg_grad_std = std(abs(left_brg-right_brg));
left_brg_mean = mean(left_brg);
right_brg_mean = mean(right_brg);


if(left_brg_mean > right_brg_mean)
    high_brg = left_brg_mean;
    low_brg = right_brg_mean;
    high_brg_std = std(left_brg);
    low_brg_std = std(right_brg);
else
    high_brg = right_brg_mean;
    low_brg = left_brg_mean;
    high_brg_std = std(right_brg);
    low_brg_std = std(left_brg);
end


% absolute curvature
abs_k = sum(K)/len;

% % compute brightness consistency along the curve
% % this computation might be problomatic when the Norm vector direction instable
% left_brg_diff = abs(diff(left_brg));
% right_brg_diff = abs(diff(right_brg));
% 
% cost = sqrt(left_brg_diff.^2+right_brg_diff.^2);
% cost = cost./bg_consis_th;
% 
% bg_consist = sum(cost.^2./(ones(size(cost))-cost.^2))/length(cost);

% lateral edge sparsity
mask = zeros(h,w); % build a mask for neighborehoodregions
for k = 1:size(V,1)
    mask(max(y(k)-nbr_width, 1) : min(y(k)+nbr_width, h), max(x(k)-nbr_width,1) : min(x(k)+nbr_width, w)) = 1;    
end



mask(y,x)=0;

nbr_edge_map = edge_map.*mask;
total_edges = sum(sum((nbr_edge_map>0)));
c_len = contour_length_mex(cfrag');
edge_sparsity = total_edges/c_len;
% keyboard;

% normalize c_len by diag of the input image
c_len = c_len/sqrt(h^2+w^2);


end