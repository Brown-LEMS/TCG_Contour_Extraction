function [sum_tp_gt, sum_tp_cp, sum_gt, sum_cp] = cfrags_vs_gt_edgemap(CEM, edgemap, gt_edgemap, maxDist, only_BB)
% match a set of curve frags to a binary ground-truth edgemap
if (nargin < 5)
    only_BB = 0;
end

[h,w] = size(gt_edgemap);
union_gt_edgemap = zeros(size(gt_edgemap));

[Y_gt,X_gt] = find(gt_edgemap>0);
min_y = min(Y_gt)-5; min_y = min(h, min_y); min_y = max(1,min_y); 
max_y = max(Y_gt)+5; max_y = min(h, max_y); max_y = max(1,max_y);
min_x = min(X_gt)-5; min_x = min(w, min_x); min_x = max(1,min_x);
max_x = max(X_gt)+5; max_x = min(w, max_x); max_x = max(1,max_x);
object_BB_mask = zeros(h,w);
object_BB_mask(min_y:max_y, min_x:max_x) = 1;
% keyboard;

sum_gt = sum(sum(gt_edgemap));
sum_tp_cp = 0;
sum_cp = 0;

for c = 1:length(CEM)
    cfrag = CEM{c};

    V = cfrag(:, 1:2);
    V = V +1; % change CXX to the matlab coordinates
    x = round(V(:,1)); x = min(w, x); x = max(1,x);
    y = round(V(:,2)); y = min(h, y); y = max(1,y);
    
    mask = zeros(h,w);
    linearInd = sub2ind(size(mask),y,x);
    mask(linearInd) = 1;
    
    if(only_BB)
       mask = mask.*object_BB_mask;
       if(sum(mask(:))==0)
           continue;
       end
    end
    
    [match1,match2] = correspondPixels(255*gt_edgemap, 255*mask, maxDist);
    union_gt_edgemap = union_gt_edgemap | match1;

%     % only consider object contour precision
%     if(sum(sum(match2>0))>0)
        sum_tp_cp = sum_tp_cp + sum(sum(match2>0));
        sum_cp = sum_cp + sum(sum(mask));
%     end
    
end

[match1,match2] = correspondPixels(255*gt_edgemap, 255*union_gt_edgemap, maxDist);
sum_tp_gt = sum(sum(match1>0));


