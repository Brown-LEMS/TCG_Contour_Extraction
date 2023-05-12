function [edges_trans, cfrags_trans, cfrags_idx_trans, cfrags_prob_trans] = transform_Contours_given_GT_disparity( cfrags_idx, edges, GT_disparity_map, is_left, cfrags_prob)

%% translate edges x coordinates
[h,w] = size(GT_disparity_map);
% x_coords = round(edges(:,1))+1;
% y_coords = round(edges(:,2))+1;


d_vec = [];
for dx = -2:2
    x_coords = round(edges(:,1))+1 + dx;
    x_coords =  max(x_coords, 1);
    x_coords = min(x_coords, w);
    for dy = -2:2
        y_coords = round(edges(:,2))+1 + dy;
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

edges_trans(:,1) =  max(0, edges_trans(:,1));
edges_trans(:,1) =  min(w-1, edges_trans(:,1));

%% Get Contours
cfrags_trans = cell(1,length(cfrags_idx));
cfrags_idx_trans =  cfrags_idx;
cfrags_prob_trans = cfrags_prob;
idx2delete = [];
for i = 1:length(cfrags_idx)
   cur_edges_org = edges(cfrags_idx{i}, :);
   cur_edges = edges_trans(cfrags_idx{i}, :);
   cfrags_trans{i} = cur_edges;
   
   %%%%%%% check artifact junctions, and remove them
   if(length(cfrags_idx{i})<=1)
       cfrags_trans{i} = [];
       cfrags_idx_trans{i} = [];
       idx2delete = [idx2delete i];
       continue;
   end
   
   if( sqrt(sum((cur_edges(1, 1:2) - cur_edges(2, 1:2)).^2)) > 5 )
%    if( abs(d_vec(cfrags_idx_trans{i}(1)) - d_vec(cfrags_idx_trans{i}(2)))  > 3 )
       cur_edges(1,:) = [];
       cfrags_idx{i}(1) = [];
   end
   if(length(cfrags_idx{i})<=1)
       cfrags_trans{i} = [];
       cfrags_idx_trans{i} = [];
       idx2delete = [idx2delete i];
       continue;
   end

   if( sqrt(sum((cur_edges(end, 1:2) - cur_edges(end-1, 1:2)).^2)) > 5 )
%    if( abs(d_vec(cfrags_idx_trans{i}(end)) - d_vec(cfrags_idx_trans{i}(end-1)))  > 3 )
       cur_edges(end,:) = [];
       cfrags_idx{i}(end) = [];
   end
   if(length(cfrags_idx{i})<=1)
       cfrags_trans{i} = [];
       cfrags_idx_trans{i} = [];
       idx2delete = [idx2delete i];
       continue;
   end
   
end

cfrags_trans = cfrags_trans(~cellfun('isempty',cfrags_trans));
cfrags_idx_trans = cfrags_idx_trans(~cellfun('isempty',cfrags_idx_trans)); 
cfrags_prob_trans(idx2delete) = [];
end