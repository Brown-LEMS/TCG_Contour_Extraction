function [bg_grad, sat_grad, hue_grad, abs_k, edge_sparsity, wigg, c_len, mean_conf, area] = curve_fragment_cues(cfrag, hsv_img, edge_map)
% cfrags is a list of edges which form a single curve fragment, edges in
% CXX coordinates,
% hsv_map is the hsv space matrix from input image
% edge_map is the binary map of detected edges before edge linking

local_dist = 2; % for general image, left and right point to compute features
nbr_width = 3; % region for edge sparsity


[h w ~] = size(hsv_img);
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
left_hue = diag(hsv_img(left_y, left_x, 1));
right_hue = diag(hsv_img(right_y, right_x, 1));
left_sat = diag(hsv_img(left_y, left_x, 2));
right_sat = diag(hsv_img(right_y, right_x, 2));
left_brg = diag(hsv_img(left_y, left_x, 3));
right_brg = diag(hsv_img(right_y, right_x, 3));
K = abs(K); % only use the abs, ignore the direction

% integral of "attribute" gradient along the curve sum of abs can 
% remove the influence of changing direction in the Norm computation
%%%%%%%%%%%%%%%% TODO: check, the abs should not be used
% bg_grad = sum(abs(left_brg-right_brg))/len;
% sat_grad = sum(abs(left_sat-right_sat))/len;
% hue_grad = sum(abs(left_hue-right_hue))/len;
bg_grad = sum((left_brg-right_brg))/len;
sat_grad = sum((left_sat-right_sat))/len;
hue_grad = sum((left_hue-right_hue))/len;

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
if(isinf(edge_sparsity))
    edge_sparsity = total_edges;
end
% keyboard;

% % normalize c_len by diag of the input image
% c_len = c_len/sqrt(h^2+w^2);

mean_conf = mean(cfrag(:,4));
area = (max(cfrag(:,1))-min(cfrag(:,1))) * (max(cfrag(:,2))-min(cfrag(:,2)));
end