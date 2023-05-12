function [edges_trans, cfrags_trans, cfrags_idx_trans, cfrags_prob_trans] = transform_Contours_given_GT_disparity_v2( cfrags_idx, edges, GT_disparity_map, is_left, cfrags_prob)
%%%%%%% This version transform contour by contour
%%%%%%% edges belonging to each contour is transformed considering spatial
%%%%%%% continuity in disparity

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
% add disparity to it;
edges = [edges d];

edges_trans =  edges;
% if(is_left)
%     edges_trans(:,1) = edges_trans(:,1) - d;
% else
%     edges_trans(:,1) = edges_trans(:,1) + d;
% end



%% Get Contours
cfrags_trans = cell(1,length(cfrags_idx));
cfrags_idx_trans =  cfrags_idx;
cfrags_prob_trans = cfrags_prob;
cid2delete = [];
moved_idx = zeros(size(edges,1), 1);
for i = 1:length(cfrags_idx)
   cur_edges = edges(cfrags_idx{i}, :);
   cur_d = cur_edges(:, end); 
   
   
%    occ_idx = cur_d==0;
%    cur_d(occ_idx) = [];
%    cur_edges(occ_idx, :) = [];
%    cfrags_idx{i}(occ_idx) = [];
   
%    smooth_d = cur_d;
   smooth_d = hampel(cur_d, 10); % remove outlier
%    smooth_d = movmax(cur_d, 5);
%    smooth_d = medfilt1(cur_d, 5);
   
   if(is_left)
       cur_edges(:,1) = cur_edges(:,1) - smooth_d;
       edges_trans(cfrags_idx{i},1) = edges(cfrags_idx{i},1) - smooth_d;     
   else
       cur_edges(:,1) = cur_edges(:,1) + smooth_d;
       edges_trans(cfrags_idx{i},1) = edges(cfrags_idx{i},1) + smooth_d;
   end
   
   % if the two times of movements of an end detach from each other, add a
   % seperate transfered edge, and update cfrags_idx_trans
   if( moved_idx(cfrags_idx{i}(1)) && abs(smooth_d(1) - cur_d(1))>10)
       edges_trans = [edges_trans; cur_edges(1,:)];
       cfrags_idx_trans{i}(1) = size(edges_trans,1);
   end
   if( moved_idx(cfrags_idx{i}(end)) && abs(smooth_d(end) - cur_d(end))>10)
       edges_trans = [edges_trans; cur_edges(end,:)];
       cfrags_idx_trans{i}(end) = size(edges_trans,1);
   end  
   
   
   cfrags_trans{i} = cur_edges;
   moved_idx(cfrags_idx{i}) = 1;
   
   %%%%%%% check artifact junctions, and remove them
   if(length(cfrags_idx{i})<=1)
       cfrags_trans{i} = [];
       cfrags_idx_trans{i} = [];
       cid2delete = [cid2delete i];
       continue;
   end
   
   if( sqrt(sum((cur_edges(1, 1:2) - cur_edges(2, 1:2)).^2)) > 10 )
%    if( abs(d_vec(cfrags_idx_trans{i}(1)) - d_vec(cfrags_idx_trans{i}(2)))  > 3 )
       cur_edges(1,:) = [];
       cfrags_idx{i}(1) = [];
   end
   if(length(cfrags_idx{i})<=1)
       cfrags_trans{i} = [];
       cfrags_idx_trans{i} = [];
       cid2delete = [cid2delete i];
       continue;
   end

   if( sqrt(sum((cur_edges(end, 1:2) - cur_edges(end-1, 1:2)).^2)) > 10 )
%    if( abs(d_vec(cfrags_idx_trans{i}(end)) - d_vec(cfrags_idx_trans{i}(end-1)))  > 3 )
       cur_edges(end,:) = [];
       cfrags_idx{i}(end) = [];
   end
   if(length(cfrags_idx{i})<=1)
       cfrags_trans{i} = [];
       cfrags_idx_trans{i} = [];
       cid2delete = [cid2delete i];
       continue;
   end
   
end

edges_trans(:,1) =  max(0, edges_trans(:,1));
edges_trans(:,1) =  min(w-1, edges_trans(:,1));
cfrags_trans = cfrags_trans(~cellfun('isempty',cfrags_trans));
cfrags_idx_trans = cfrags_idx_trans(~cellfun('isempty',cfrags_idx_trans)); 
cfrags_prob_trans(cid2delete) = [];

%%%%%%%%%%% iterate through all curves, break at discontinuous spots
cid2delete = [];
for cid = 1:length(cfrags_idx_trans)

    cur_c = cfrags_trans{cid};
    cur_c_idx = cfrags_idx_trans{cid};
    c_size = length(cur_c_idx);
    cur_prob = cfrags_prob_trans(cid);
    
    % look for discontinuous spots
    pairwise_dists= sqrt(diff(cur_c(:,1)).^2 + diff(cur_c(:,2)).^2);
    cur_break_e_ids = find(pairwise_dists>10);
    
    %%%%%%%%%%%%%  break the cfrag at junctions (one or more)
    %%%%%%%%%%%%% iteratively update at cfrags, finally replace new
    if(~isempty(cur_break_e_ids))
%         cur_break_e_ids = sort(cur_break_e_ids, 'ascend');

        if(c_size - cur_break_e_ids(end)>2)
            cfrags_trans{cid} = cur_c(cur_break_e_ids(end)+1:end, :);
            cfrags_idx_trans{cid} = cur_c_idx(cur_break_e_ids(end)+1:end);
        else
            cid2delete = [cid2delete cid];
        end

        if(cur_break_e_ids(1)>2)
            cfrags_trans = [cfrags_trans cur_c(1:cur_break_e_ids(1), :)];
            cfrags_idx_trans = [cfrags_idx_trans cur_c_idx(1:cur_break_e_ids(1))];
            cfrags_prob_trans = [cfrags_prob_trans cur_prob];
        end
        
        %%%%%%%% break at multiple pts
        prev_e_id = cur_break_e_ids(1);
        if(length(cur_break_e_ids)>=2)
            for j = 1:length(cur_break_e_ids)-1  
                idx_to_add = prev_e_id+1:cur_break_e_ids(j+1);
                
                if(length(idx_to_add)>2)
                    cfrags_trans = [cfrags_trans cur_c(idx_to_add, :)];
                    cfrags_idx_trans = [cfrags_idx_trans cur_c_idx(idx_to_add)];
                    cfrags_prob_trans = [cfrags_prob_trans cur_prob];
                end

                prev_e_id = cur_break_e_ids(j+1);
            end
        end
    end
end
cfrags_trans(cid2delete) = [];
cfrags_idx_trans(cid2delete) = [];
cfrags_prob_trans(cid2delete) = [];

cfrags_trans = cfrags_trans(~cellfun('isempty',cfrags_trans));
cfrags_idx_trans = cfrags_idx_trans(~cellfun('isempty',cfrags_idx_trans)); 

end