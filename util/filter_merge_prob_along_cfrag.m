function merge_prob_vec = filter_merge_prob_along_cfrag(cfrag, hsv_img, edge_map, tmap, beta, fmean, nbr_range_th)

%%%%%%%%%%%%%%%%   IMPORTANT (diff from later process):
%%%%%%%%%%%%%%%%   Down sampled to pixel-level density, compute
%%%%%%%%%%%%%%%%   all the features and probability, then project back to
%%%%%%%%%%%%%%%%   the original sample density: a lot of region-based
%%%%%%%%%%%%%%%%   features do not need sub-pixel localization


% merge_prob_vec = ones(size(cfrag,1),1);
    
local_dist = 1; % for general image
% bg_consis_th = 2;% sensitive value
nbr_width = 3; % region for edge sparsity

[h w ~] = size(hsv_img);

V = cfrag(:, 1:2);
V = V +1; % change CXX to the matlab coordinates

%% if edges are subpixel level dense,  down sample to pixel level density (important)
V = round(V);
[V, ia, ib] = unique(V, 'rows', 'stable');
len = size(V,1);
merge_prob_vec_down_sampled = ones(len,1);

%%%%%%%%%%%%%%%% compute curveture vector
K = LineCurvature2D(V);

if( sum(isnan(K))>0)
%     keyboard;
    K(find(isnan(K))) =0;
end

%%%%%%%%%%%%%%%%%% wiggliness is the time curvature change sign

sign_vec = K(1:end-1).*K(2:end);
sign_vec = [sign_vec;0];

%%%%%%%%%%%%  compute feature given coordinate of left and right points

N = LineNormals2D(V); % N's direction are normalized to unit already

x = V(:,1); x = min(w, x); x = max(1,x);
y = V(:,2); y = min(h, y); y = max(1,y);

% N's direction matters
left_x = round(x - local_dist*N(:,1));
left_y = round(y - local_dist*N(:,2));
right_x = round(x + local_dist*N(:,1));
right_y = round(y + local_dist*N(:,2));

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

%% Compute interest features using integral vector
% left_hue_cum = cumsum(left_hue);
% right_hue_cum = cumsum(right_hue);
% left_sat_cum = cumsum(left_sat);
% right_sat_cum = cumsum(right_sat);
% left_brg_cum = cumsum(left_brg);
% right_brg_cum = cumsum(right_brg);

 %************** INDEED, abs should not be used, the grads with different signs differs
bg_grad_cum = cumsum((left_brg - right_brg));
sat_grad_cum = cumsum((left_sat - right_sat));
hue_grad_cum = cumsum((left_hue - right_hue));
wigg_cum = cumsum(sign_vec<0);
K_cum = cumsum(K);
[~, edge_sparsity_cum] = compute_edge_sparsity_integral(x,y,N, nbr_width, edge_map);
[~, ~, texton_hist_left_cum, texton_hist_right_cum] = compute_texture_hist_integral(x,y,N, nbr_width, tmap);

%% compute prob along each curve
%%%%%%%%%%%%%%%%%%  TODO: can be converted into matrix
%%%%%%%%%%%%%%%%%%  computation to get rid of loop
for i = nbr_range_th +1: len-nbr_range_th
    
    %%%%%%%%%%%% compute c1 cues
    bg_grad_c1 = (bg_grad_cum(i) - bg_grad_cum(i-nbr_range_th))/nbr_range_th;
    sat_grad_c1 = (sat_grad_cum(i) - sat_grad_cum(i-nbr_range_th))/nbr_range_th;
    hue_grad_c1 = (hue_grad_cum(i) - hue_grad_cum(i-nbr_range_th))/nbr_range_th;
    abs_k_c1 = (K_cum(i) - K_cum(i-nbr_range_th))/nbr_range_th;
    edge_sparsity_c1 = (edge_sparsity_cum (i) - edge_sparsity_cum (i-nbr_range_th)) / nbr_range_th;
    wigg_c1 = (wigg_cum(i) - wigg_cum(i-nbr_range_th))/nbr_range_th;
    
    c1_cues = [bg_grad_c1; sat_grad_c1; hue_grad_c1; abs_k_c1; edge_sparsity_c1; wigg_c1];

    %%%%%%%%%%%%% compute c2 cues
    bg_grad_c2 = (bg_grad_cum(i+nbr_range_th) - bg_grad_cum(i))/nbr_range_th;
    sat_grad_c2 = (sat_grad_cum(i+nbr_range_th) - sat_grad_cum(i))/nbr_range_th;
    hue_grad_c2 = (hue_grad_cum(i+nbr_range_th) - hue_grad_cum(i))/nbr_range_th;
    abs_k_c2 = (K_cum(i+nbr_range_th) - K_cum(i))/nbr_range_th;
    edge_sparsity_c2 = (edge_sparsity_cum(i+nbr_range_th) -  edge_sparsity_cum(i))/nbr_range_th;
    wigg_c2 = (wigg_cum(i+nbr_range_th) - wigg_cum(i))/nbr_range_th;
    
    c2_cues = [bg_grad_c2; sat_grad_c2; hue_grad_c2; abs_k_c2; edge_sparsity_c2; wigg_c2];
    
    
    %%%%%%%%%%%% compute cues difference
    cues_diff = abs(c1_cues-c2_cues);
    
    
    %%%%%%%%%%%% compute texture difference
    hist_1_left = (texton_hist_left_cum(:,i) - texton_hist_left_cum(:,i-nbr_range_th))/nbr_range_th;
    hist_2_left = (texton_hist_left_cum(:,i+nbr_range_th) - texton_hist_left_cum(:, i))/nbr_range_th;
    left_dist = pdist2(hist_1_left',hist_2_left','chisq');
    
    hist_1_right = (texton_hist_right_cum(:,i) - texton_hist_right_cum(:,i-nbr_range_th))/nbr_range_th;
    hist_2_right = (texton_hist_right_cum(:,i+nbr_range_th) - texton_hist_right_cum(:, i))/nbr_range_th;
    right_dist = pdist2(hist_1_right',hist_2_right','chisq');    
    
    texture_diff = left_dist + right_dist;
    
    %%%%%%%%%%%% compute ori diff at the connecting points
%     [geom_diff, texture_diff] = degree_2_node_cues(cfrag(idx_c1, :),cfrag(idx_c2, :), tmap);
    % the range of ori diff in semantic should also be adaptive
    ori_1 = V(i, 1:2) - V(i-nbr_range_th, 1:2);
    ori_2 = V(i+nbr_range_th, 1:2) - V(i, 1:2);

    geom_diff = ori_1*ori_2'/sqrt(ori_1*ori_1')/sqrt(ori_2*ori_2');
    
    
    %%%%%%%%%%% compute the probability given difference between cues
    cues_diff = [cues_diff; geom_diff; texture_diff];

     p = 1 / (1 + exp(-([1 cues_diff']-fmean)*beta'));
     merge_prob_vec_down_sampled(i) = p;
end

%%%%%%%%%%%%% covert back to the original density

merge_prob_vec = merge_prob_vec_down_sampled(ib);


end

% function edge_sparsity = compute_edge_sparsity(x,y, nbr_width, edge_map, cfrag)
%     % lateral edge sparsity
%     [h,w, ~] = size(edge_map);
%     V = cfrag(:, 1:2);
%     V = V +1; % change CXX to the matlab coordinates
%     mask = zeros(size(edge_map)); % build a mask for neighborehood regions
%     for k = 1:size(V,1)
%         mask(max(y(k)-nbr_width, 1) : min(y(k)+nbr_width, h), max(x(k)-nbr_width,1) : min(x(k)+nbr_width, w)) = 1;    
%     end
% 
%     mask(y,x)=0;
% 
%     nbr_edge_map = edge_map.*mask;
%     total_edges = sum(sum((nbr_edge_map>0)));
%     c_len = contour_length_mex(cfrag');
%     edge_sparsity = total_edges/c_len;    
%     
% end
