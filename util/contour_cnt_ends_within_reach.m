function [new_cfrags, new_cfrags_idx, T_junctions, cnt_ends_within_reach, free_end_map] = contour_cnt_ends_within_reach(cfrags, cfrags_idx, h, w, edges)
% locate the conners and break the contour
search_nbr_range = 20; % search range from an end of curve fragment 
jct_nbr_range = 5; % supression range of an existing junction
search_ori_range = pi/18; % search angle at an end point of cfrag
approach_clen_th = 2; % need to be small or the gap will not be filled for very broken short cfrags
c_len_th = 20;
nbr_xx_coordinates = repmat(-search_nbr_range:search_nbr_range, [2*search_nbr_range+1, 1]);
nbr_yy_coordinates = nbr_xx_coordinates';
nbr_rr = sqrt(nbr_xx_coordinates.^2 + nbr_yy_coordinates.^2);
nbr_theta = atan2(nbr_yy_coordinates, nbr_xx_coordinates);

% max_iter = params.max_iter;

disp('introduce junctions at gaps');

G_test = construct_fac_graph_from_curve_fragments (cfrags_idx, cfrags);
% introduced_num_junction_points = 0;
T_junctions = [];
num_org_cfrags = length(cfrags);
%% look for junctions, break the cfrag, extent the approaching cfrag

%%%%%%%%%%%%%%%%%%%%% contruct tables refering end pts %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% extract start/end pts:  format: x, y, theta
start_pts = zeros(num_org_cfrags, 3);
end_pts = zeros(num_org_cfrags, 3);
c_lens = zeros(num_org_cfrags, 1);
ref_table_start_pts = zeros([h, w]); % save the contour id of each start point
ref_table_end_pts = zeros([h, w]); % save the contour id of each end point
free_end_search_coods = cell(1,length(G_test.var));
free_end_map = zeros([h, w]);
is_free_end = zeros(1, length(G_test.var));
is_starting_point = zeros(1, length(G_test.var));

cnt_ends_within_reach = [];

%%%%%%%%%% construct free end binary map
for vid = 1:length(G_test.var)
   if(length(G_test.var(vid).nbrs_fac)==1)
        actual_edge = edges(G_test.var(vid).actual_edge_id,:);
        free_end_map(round(actual_edge(1,2)+1),round(actual_edge(1,1)+1)) = vid;
   end
end




for cid = 1: num_org_cfrags
    cur_c = cfrags{cid};
    
    c_len = contour_length_mex(cur_c');
    c_lens(cid) = c_len;
    start_pts(cid,:) = cur_c(1,1:3);
    end_pts(cid,:) = cur_c(end,1:3);

    %%%%%% skip circles
    if(cur_c(1,:) == cur_c(end,:))
        continue;
    end
    
    %%%%%% only consider length>5 contour as approaching cfrags
    if(c_len > approach_clen_th || size(cur_c, 1) >approach_clen_th)
        e_size = size(cur_c,1);

        %%%%%% constrain the angle of area from a end
        angle = atan2(cur_c(1,2) - cur_c(min(e_size,3),2), cur_c(1,1) - cur_c(min(e_size,3),1));
        angle_max = angle+search_ori_range;
        angle_min = angle-search_ori_range;
        
        if(angle_max>pi)
            qualify_idx = find((nbr_theta(:) >= angle_min & nbr_theta(:) <= pi) | ...
                            (nbr_theta(:) > -pi & nbr_theta(:) <= angle_max-2*pi) & nbr_rr(:)<=search_nbr_range);
        elseif(angle_min<-pi)
            qualify_idx = find((nbr_theta(:) >= angle_min +2*pi & nbr_theta(:) <= pi) | ...
                            (nbr_theta(:) > -pi & nbr_theta(:) <= angle_max) & nbr_rr(:)<=search_nbr_range);
        else
            qualify_idx = find(nbr_theta(:) >= angle_min & nbr_theta(:) <= angle_max  & nbr_rr(:)<=search_nbr_range);
        end
        x_coords_start = nbr_xx_coordinates(qualify_idx)+round(start_pts(cid,1)+1);
        y_coords_start = nbr_yy_coordinates(qualify_idx)+round(start_pts(cid,2)+1);
        x_coords_start = max(x_coords_start, ones(size(x_coords_start)));
        y_coords_start = max(y_coords_start, ones(size(y_coords_start)));
        x_coords_start = min(x_coords_start, w*ones(size(x_coords_start)));
        y_coords_start = min(y_coords_start, h*ones(size(y_coords_start)));
        % only consider free end
        start_vid = G_test.fac(cid).nbrs_var(1);
        if(length(G_test.var(start_vid).nbrs_fac)==1)% && sum(sum(jct_map(y_coords_start, x_coords_start)))>0)
            
            % look for covered free_end
            ends_within_reach = free_end_map(sub2ind([h,w], y_coords_start, x_coords_start));
            cnt_ends_within_reach = [cnt_ends_within_reach; sum(ends_within_reach>0)];
            
            ref_table_start_pts(sub2ind([h,w], y_coords_start, x_coords_start)) = cid;
            is_free_end(start_vid)=1;
            is_starting_point(start_vid)=1;
            free_end_search_coods{start_vid} = [x_coords_start y_coords_start];
        end
        
        %%%%%% constrain the angle of area from a end
        angle = atan2(cur_c(end,2) - cur_c(max(e_size-2, 1),2), cur_c(end,1) - cur_c(max(e_size-2, 1),1));
        angle_max = angle+search_ori_range;
        angle_min = angle-search_ori_range;
        
        if(angle_max>pi)
            qualify_idx = find((nbr_theta(:) >= angle_min & nbr_theta(:) <= pi) | ...
                            (nbr_theta(:) > -pi & nbr_theta(:) <= angle_max-2*pi) & nbr_rr(:)<=search_nbr_range);
        elseif(angle_min<-pi)
            qualify_idx = find((nbr_theta(:) >= angle_min +2*pi & nbr_theta(:) <= pi) | ...
                            (nbr_theta(:) > -pi & nbr_theta(:) <= angle_max) & nbr_rr(:)<=search_nbr_range);
        else
            qualify_idx = find(nbr_theta(:) >= angle_min & nbr_theta(:) <= angle_max & nbr_rr(:)<=search_nbr_range);
        end
        x_coords_end = nbr_xx_coordinates(qualify_idx)+round(end_pts(cid,1)+1);
        y_coords_end = nbr_yy_coordinates(qualify_idx)+round(end_pts(cid,2)+1);
        x_coords_end = max(x_coords_end, ones(size(x_coords_end)));
        y_coords_end = max(y_coords_end, ones(size(y_coords_end)));
        x_coords_end = min(x_coords_end, w*ones(size(x_coords_end)));
        y_coords_end = min(y_coords_end, h*ones(size(y_coords_end)));
        % only consider free end
        end_vid = G_test.fac(cid).nbrs_var(2);
        if(length(G_test.var(end_vid).nbrs_fac)==1)% && sum(sum(jct_map(y_coords_end, x_coords_end)))>0)
            
            % look for covered free_end
            ends_within_reach = free_end_map(sub2ind([h,w], y_coords_end, x_coords_end));
            cnt_ends_within_reach = [cnt_ends_within_reach; sum(ends_within_reach>0)];
            
            ref_table_end_pts(sub2ind([h,w], y_coords_end, x_coords_end)) = cid;  
            is_free_end(end_vid)=1;
            is_starting_point(end_vid)=0;
            free_end_search_coods{end_vid} = [x_coords_end y_coords_end];
        end
    end
end


new_cfrags = cfrags;
new_cfrags_idx = cfrags_idx;
% only keep the non-empty contours
new_cfrags = new_cfrags(~cellfun('isempty',new_cfrags));
new_cfrags_idx = new_cfrags_idx(~cellfun('isempty',new_cfrags_idx)); 


end