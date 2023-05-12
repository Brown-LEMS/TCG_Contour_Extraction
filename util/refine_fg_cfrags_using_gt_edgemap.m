function [fg_contours, fg_edgemap, bg_contours, fg_probs, bg_probs] = refine_fg_cfrags_using_gt_edgemap(CEM, gt_edgemap, patio_th, len_th)
% only care about contours with len_th or more edges, for both fg and bg
% patio_th = 0.3;
fg_contours = cell(1, 0);
bg_contours = cell(1, 0);
fg_probs = [];
bg_probs = [];

w = CEM{1,1}(1);
h = CEM{1,1}(2);


gt_edgemap = (gt_edgemap>0);
fg_edgemap = zeros(size(gt_edgemap));

for c = 1:length(CEM{2,1})
    cfrag = CEM{2,1}{c};
    if(size(cfrag, 1)< len_th)
        continue;
    end
%     c_len = contour_length_mex(cfrag');
%     if(c_len>80)
%         continue;
%     end
    V = cfrag(:, 1:2);
    V = V +1; % change CXX to the matlab coordinates
    x = round(V(:,1)); x = min(w, x); x = max(1,x);
    y = round(V(:,2)); y = min(h, y); y = max(1,y);
    mask = zeros(h,w);
    linearInd = sub2ind(size(mask),y,x);
    mask(linearInd) = 1;
    pixel_len = sum(sum(mask));
    tp = sum(sum(mask.*gt_edgemap));
%     imshow(mask);
    if(tp/pixel_len >patio_th)
        fg_contours = [fg_contours, cfrag];
        fg_probs = [fg_probs tp/pixel_len];
        fg_edgemap = fg_edgemap + mask;
    else
        bg_contours = [bg_contours, cfrag];
        bg_probs = [bg_probs tp/pixel_len];
    end
    
end

fg_edgemap = double(fg_edgemap>0);

% imshow(fg_edgemap);



