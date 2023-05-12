function Y = b_merge_node(c1, c2, gt_edgemap_objects, maxDist)
[h,w] = size(gt_edgemap_objects);

V1 = c1(:, 1:2);
V1 = V1 +1; % change CXX to the matlab coordinates
x1 = round(V1(:,1)); x1 = min(w, x1); x1 = max(1,x1);
y1 = round(V1(:,2)); y1 = min(h, y1); y1 = max(1,y1);
mask1 = zeros(h,w);
linearInd = sub2ind(size(mask1),y1,x1);
mask1(linearInd) = 1;

[match1,match2] = correspondPixels(255*(gt_edgemap_objects>0), 255*mask1, maxDist);

[gt_y1, gt_x1] = find(match1>0);
linearInd = sub2ind(size(match1),gt_y1,gt_x1);
object_ids_1 = gt_edgemap_objects(linearInd);


V2 = c2(:, 1:2);
V2 = V2 +1; % change CXX to the matlab coordinates
x2 = round(V2(:,1)); x2 = min(w, x2); x2 = max(1,x2);
y2 = round(V2(:,2)); y2 = min(h, y2); y2 = max(1,y2);
mask2 = zeros(h,w);
linearInd = sub2ind(size(mask2),y2,x2);
mask2(linearInd) = 1;

[match1,match2] = correspondPixels(255*(gt_edgemap_objects>0), 255*mask2, maxDist);

[gt_y2, gt_x2] = find(match1>0);
linearInd = sub2ind(size(match1),gt_y2,gt_x2);
object_ids_2 = gt_edgemap_objects(linearInd);

% keyboard;

if(isempty(object_ids_1) && isempty(object_ids_2))
    Y = -1;
elseif(isempty(object_ids_1) && ~isempty(object_ids_2))
    Y = -1;
elseif(~isempty(object_ids_1) && isempty(object_ids_2))
    Y = -1;
elseif(~isempty(object_ids_1) && ~isempty(object_ids_2))
    if(length(unique(object_ids_1)) ~= length(unique(object_ids_2)))
        Y = 0;
    elseif(size(unique(object_ids_1)) == size(unique(object_ids_2))  )
        if(unique(object_ids_1) == unique(object_ids_2))
            Y = 1;
        else
            Y = 0;
        end
    end
else
    Y = -1;
end


end