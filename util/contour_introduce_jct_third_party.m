function [new_cfrags, new_cfrags_idx, junction_pts] = contour_introduce_jct_third_party(cfrags, cfrags_idx, jct_map, edges, params)
% locate the conners and break the contour
ref_nbr_range = 5; % search range from an end of curve fragment 
jct_nbr_range = 5; % supression range of an existing junction
ori_range = pi/3;
approach_clen_th = 5;
nbr_xx_coordinates = repmat(-ref_nbr_range:ref_nbr_range, [2*ref_nbr_range+1, 1]);
nbr_yy_coordinates = nbr_xx_coordinates';
nbr_rr = sqrt(nbr_xx_coordinates.^2 + nbr_yy_coordinates.^2);
nbr_theta = atan2(nbr_yy_coordinates, nbr_xx_coordinates);

diag_ratio= params.diag_ratio;
merge_th_geom = params.merge_th_geom; 
nbr_num_edges = params.nbr_num_edges;
beta_1 = params.beta_1;
fmean_1 = params.fmean_1;
% max_iter = params.max_iter;

disp('introduce junctions based on complimentary junction deteciton');
[h,w] = size(jct_map);

G_test = construct_fac_graph_from_curve_fragments (cfrags_idx, cfrags);
% introduced_num_junction_points = 0;
junction_pts = [];
num_org_cfrags = length(cfrags);
%% look for junctions, break the cfrag, extent the approaching cfrag

%%%%%%%%%%%%%%%%%%%%% contruct tables refering end pts %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% extract start/end pts:  format: x, y, theta
start_pts = zeros(num_org_cfrags, 3);
end_pts = zeros(num_org_cfrags, 3);
c_lens = zeros(num_org_cfrags, 1);
ref_table_start_pts = zeros(size(jct_map)); % save the contour id of each start point
ref_table_end_pts = zeros(size(jct_map)); % save the contour id of each end point
free_end_search_coods = cell(1,length(G_test.var));
free_end_map = zeros(size(jct_map));
is_free_end = zeros(1, length(G_test.var));
is_starting_point = zeros(1, length(G_test.var));

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
    if(c_len > approach_clen_th*diag_ratio || size(cur_c, 1) >approach_clen_th)
        
        %%%%%% constrain the angle of area from a end
        angle = atan2(cur_c(1,2) - cur_c(2,2), cur_c(1,1) - cur_c(2,1));
        angle_max = angle+ori_range;
        angle_min = angle-ori_range;
        
        if(angle_max>pi)
            qualify_idx = find((nbr_theta(:) >= angle_min & nbr_theta(:) <= pi) | ...
                            (nbr_theta(:) > -pi & nbr_theta(:) <= angle_max-2*pi) & nbr_rr(:)<=ref_nbr_range);
        elseif(angle_min<-pi)
            qualify_idx = find((nbr_theta(:) >= angle_min +2*pi & nbr_theta(:) <= pi) | ...
                            (nbr_theta(:) > -pi & nbr_theta(:) <= angle_max) & nbr_rr(:)<=ref_nbr_range);
        else
            qualify_idx = find(nbr_theta(:) >= angle_min & nbr_theta(:) <= angle_max  & nbr_rr(:)<=ref_nbr_range);
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
            ref_table_start_pts(sub2ind([h,w], y_coords_start, x_coords_start)) = cid;
            free_end_map(round(start_pts(cid,2)+1),round(start_pts(cid,1)+1)) = start_vid;
            is_free_end(start_vid)=1;
            is_starting_point(start_vid)=1;
            free_end_search_coods{start_vid} = [x_coords_start y_coords_start];
        end
        
        %%%%%% constrain the angle of area from a end
        angle = atan2(cur_c(end,2) - cur_c(end-1,2), cur_c(end,1) - cur_c(end-1,1));
        angle_max = angle+ori_range;
        angle_min = angle-ori_range;
        
        if(angle_max>pi)
            qualify_idx = find((nbr_theta(:) >= angle_min & nbr_theta(:) <= pi) | ...
                            (nbr_theta(:) > -pi & nbr_theta(:) <= angle_max-2*pi) & nbr_rr(:)<=ref_nbr_range);
        elseif(angle_min<-pi)
            qualify_idx = find((nbr_theta(:) >= angle_min +2*pi & nbr_theta(:) <= pi) | ...
                            (nbr_theta(:) > -pi & nbr_theta(:) <= angle_max) & nbr_rr(:)<=ref_nbr_range);
        else
            qualify_idx = find(nbr_theta(:) >= angle_min & nbr_theta(:) <= angle_max & nbr_rr(:)<=ref_nbr_range);
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
            ref_table_end_pts(sub2ind([h,w], y_coords_end, x_coords_end)) = cid;  
            free_end_map(round(end_pts(cid,2)+1),round(end_pts(cid,1)+1)) = end_vid;
            is_free_end(end_vid)=1;
            is_starting_point(end_vid)=0;
            free_end_search_coods{end_vid} = [x_coords_end y_coords_end];
        end
    end
end

%%%%%%%%%%%%%%%%%%%%% construct the distance transform of existing junctions
ref_table_jct_pts = zeros(size(jct_map)); % save the vid
for vid = 1:length(G_test.var)
    if(length(G_test.var(vid).nbrs_fac)<3)
        continue;
    end
    
    cur_edge = edges(G_test.var(vid).actual_edge_id, :);
    x_min = max(1, round(cur_edge(1))+1-jct_nbr_range);
    x_max = min(w, round(cur_edge(1))+1+jct_nbr_range);
    y_min = max(1, round(cur_edge(2))+1-jct_nbr_range);
    y_max = min(h, round(cur_edge(2))+1+jct_nbr_range);
    ref_table_jct_pts(y_min:y_max, x_min:x_max) = vid;
end

%%%%%%%%%%%%%%%%%%%%% check for the case junction that is formed by connecting 3 end points
for vid = 1:length(G_test.var) 
    if(~is_free_end(vid))
        continue;
    end
    cur_eid = G_test.var(vid).actual_edge_id;
    cur_edge = edges(cur_eid, :);
    xx = round(cur_edge(1))+1;
    yy = round(cur_edge(2))+1;
    %%%% the junction to insert should be detected in jct_map
    x_min = max(1, xx-ref_nbr_range);
    x_max = min(w, xx+ref_nbr_range);
    y_min = max(1, yy-ref_nbr_range);
    y_max = min(h, yy+ref_nbr_range);
%     if(sum(sum(jct_map(y_min:y_max, x_min:x_max)))==0)
%        continue; 
%     end
    %%%% check if there are >=2 end points in the neighbore
    nbr_end_vids = unique(free_end_map(y_min:y_max, x_min:x_max));
    nbr_end_vids(nbr_end_vids==0) = [];
    nbr_end_vids(nbr_end_vids==vid) = [];
    if(length(nbr_end_vids)<2)
       continue; 
    end
    
    %%%% check geometrically qualified approaching ends
    cur_c = cfrags{G_test.var(vid).nbrs_fac};
    if(is_starting_point(vid))
        dx = cur_c(2,1) - cur_c(1,1);
        dy = cur_c(2,2) - cur_c(1,2);
    else
        dx = cur_c(end-1,1) - cur_c(end,1);
        dy = cur_c(end-1,2) - cur_c(end,2);        
    end
    exist_dx = dx;
    exist_dy = dy;
    qualify_end_vids = [];
%     single_cos = 0;
    for vv = nbr_end_vids'
        % skip the used ones
        if(~is_free_end(vv))
            continue;
        end
        
        % current end is not within a nbr end's search range
        if(sum(free_end_search_coods{vv}(:,1)==xx & free_end_search_coods{vv}(:,2)==yy)==0)
            continue;
        end
        % is overlapping with any existing branch 
        approach_edge = edges(G_test.var(vv).actual_edge_id, :);
        dx = approach_edge(1)-cur_edge(1);
        dy = approach_edge(2)-cur_edge(2);
        dx = repmat(dx, size(exist_dx));
        dy = repmat(dy, size(exist_dy));
        cos_vec = (exist_dx.*dx + exist_dy.*dy)./ ...
                   sqrt(exist_dx.^2 + exist_dy.^2)./sqrt(dx.^2 + dy.^2);
        if(sum(cos_vec > cos(pi/4))) % has a branch < 45 degree
            continue;
        end
        exist_dx = [exist_dx; dx];
        exist_dy = [exist_dy; dy];
        qualify_end_vids = [qualify_end_vids vv];
%         single_cos = cos_vec(1);
    end
    
    if(length(qualify_end_vids)<2)
       continue; 
    end
    
%     if(length(qualify_end_vids)==1 && single_cos>-0.85) % the filling gap need to be super smooth
%         continue;
%     end
    
    %%%% form the junction
    for vv = qualify_end_vids
        cid = G_test.var(vv).nbrs_fac;
        if(is_starting_point(vv))
            cfrags{cid} = [cur_edge; cfrags{cid}];
            cfrags_idx{cid} = [cur_eid cfrags_idx{cid}];
            ref_table_start_pts(ref_table_start_pts==cid) = 0;
        else
            cfrags{cid} = [cfrags{cid}; cur_edge];
            cfrags_idx{cid} = [cfrags_idx{cid} cur_eid];
            ref_table_end_pts(ref_table_end_pts==cid) = 0;
        end
    end    
    
    % update markers and maps
    junction_pts = [junction_pts; cur_edge];
    ref_table_jct_pts(y_min:y_max, x_min:x_max) = vid;
    
    
    is_free_end(vid) = 0;
    is_free_end(qualify_end_vids) = 0;
    
end

%%%%%%%%%%%%%%%%%%%%% along each cfrag, look for junctions, break the cfrag, extent the approaching cfrag
for i = 1: num_org_cfrags
    cur_c = cfrags{i};
    c_len = c_lens(i);

    cur_c_idx = cfrags_idx{i};
    % skip closed circle
    if(cur_c(1,:) == cur_c(end,:))
        continue;
    end

    % only introducing breaking points for curves which are long enough
    if(c_len>(10 * diag_ratio) && length(cur_c_idx) > nbr_num_edges+1)
        x_coords = round(cur_c(:,1))+1;
        y_coords = round(cur_c(:,2))+1;
        x_coords = max(x_coords, ones(size(x_coords)));
        y_coords = max(y_coords, ones(size(y_coords)));
        x_coords = min(x_coords, w*ones(size(x_coords)));
        y_coords = min(y_coords, h*ones(size(y_coords)));          

        start_id_vec = ref_table_start_pts(sub2ind([h,w], y_coords, x_coords));
        end_id_vec = ref_table_end_pts(sub2ind([h,w], y_coords, x_coords));

        unique_start_id = unique(start_id_vec);
        unique_end_id = unique(end_id_vec);

        cur_break_e_ids = []; % save the break e_ids into a vector

        prev_e_id = 1;
        for j = 1:length(unique_start_id)
            %%%%% id of approaching cfrag
            c_id = unique_start_id(j);
            
            if(c_id==0 || c_id ==i)
                continue;
            end
               
            %%%%% look for the junction edge, choose the midian
            e_ids = find(start_id_vec==c_id);
            %%%%% find the e_id with min dist to the end pt
            e_locs = cur_c(e_ids, 1:2);
            a_loc = repmat(cfrags{c_id}(1, 1:2), [length(e_ids),1]);
            loc_diff = e_locs - a_loc;
            loc_diff = sqrt(loc_diff(:,1).^2 + loc_diff(:,2).^2);
            [min_dist, min_id] = min(loc_diff);
            e_id = e_ids(min_id);

%             %%%% the junction to insert should be detected in jct_map
%             x_min = max(1, round(cur_c(e_id, 1))+1-jct_nbr_range);
%             x_max = min(w, round(cur_c(e_id, 1))+1+jct_nbr_range);
%             y_min = max(1, round(cur_c(e_id, 2))+1-jct_nbr_range);
%             y_max = min(h, round(cur_c(e_id, 2))+1+jct_nbr_range);
%             if(sum(sum(jct_map(y_min:y_max, x_min:x_max)))==0)
%                continue; 
%             end
            
            %%%% skip locations close to existing junctions
            if(ref_table_jct_pts(round(cur_c(e_id, 2))+1, round(cur_c(e_id, 1))+1))
                continue;
            end
            
            % if there are multiple junctions to insert, they should beyond
            % some distanct
            if(abs(e_id-prev_e_id) < 3)
               continue; 
            end
            
            % when the intersection super close to an end point
            if(e_id<3 )
%                 % simply extend the approaching cfrag to the end
%                 cfrags{c_id} = [cur_c(1, :); cfrags{c_id}];
%                 cfrags_idx{c_id} = [cur_c_idx(1) cfrags_idx{c_id}];
%                 junction_pts = [junction_pts; cur_c(1, :)];
                continue;
            end
            if(e_id > size(cur_c,1)-2)
%                 % simply extend the approaching cfrag to the end
%                 cfrags{c_id} = [cur_c(end, :); cfrags{c_id}];
%                 cfrags_idx{c_id} = [cur_c_idx(end) cfrags_idx{c_id}];
%                 junction_pts = [junction_pts; cur_c(end, :)];
                continue;
            end
            

            % Break existing T-junction
            if(min_dist == 0)
                junction_pts = [junction_pts; cur_c(e_id, :)];
%                 introduced_num_junction_points = introduced_num_junction_points + 1; 
                cur_break_e_ids =  [cur_break_e_ids e_id];
                continue;
            end

            % only consider junction when abs(cos(ori_diff)) < cos(pi/6);
            if(e_id==1)
                c_ori_vec = cur_c(e_id + 1, 1:2) - cur_c(e_id, 1:2);
            elseif (e_id==size(cur_c,1))
                c_ori_vec = cur_c(e_id, 1:2) - cur_c(e_id-1, 1:2);  
            else
                c_ori_vec = cur_c(e_id + 1, 1:2) - cur_c(e_id-1, 1:2);
            end
            a_ori_vec = cfrags{c_id}(min(5, end), 1:2) - cfrags{c_id}(1, 1:2);

            cos_ori_diff = abs(c_ori_vec*a_ori_vec'/norm(c_ori_vec)/norm(a_ori_vec));
            if( cos_ori_diff > cos(pi/6))
                continue;
            end

            %%%%% extend the approaching cfrag to the intersection
            cfrags{c_id} = [cur_c(e_id, :);  cfrags{c_id}];
            cfrags_idx{c_id} = [cur_c_idx(e_id) cfrags_idx{c_id}];
            
            junction_pts = [junction_pts; cur_c(e_id, :)];
%             introduced_num_junction_points = introduced_num_junction_points + 1; 
            cur_break_e_ids =  [cur_break_e_ids e_id];
            prev_e_id = e_id;
            % update ref maps
            ref_table_start_pts(ref_table_start_pts==c_id) = 0;
        end

        prev_e_id = 1;
        for j = 1:length(unique_end_id)
            %%%%% id of approaching cfrag
            c_id = unique_end_id(j);
            if(c_id==0 || c_id ==i)
                continue;
            end

            %%%%% look for the junction edge, choose the midian
            e_ids = find(end_id_vec==c_id);
            %%%%% find the e_id with min dist to the end pt
            e_locs = cur_c(e_ids, 1:2);
            a_loc = repmat(cfrags{c_id}(end, 1:2), [length(e_ids),1]);

            loc_diff = e_locs - a_loc;
            loc_diff = sqrt(loc_diff(:,1).^2 + loc_diff(:,2).^2);
            [min_dist, min_id] = min(loc_diff);
            e_id = e_ids(min_id);

%             %%%% the junction to insert should be detected in jct_map
%             x_min = max(1, round(cur_c(e_id, 1))+1-jct_nbr_range);
%             x_max = min(w, round(cur_c(e_id, 1))+1+jct_nbr_range);
%             y_min = max(1, round(cur_c(e_id, 2))+1-jct_nbr_range);
%             y_max = min(h, round(cur_c(e_id, 2))+1+jct_nbr_range);
%             if(sum(sum(jct_map(y_min:y_max, x_min:x_max)))==0)
%                continue; 
%             end
            
            %%%% skip locations close to existing junctions
            if(ref_table_jct_pts(round(cur_c(e_id, 2))+1, round(cur_c(e_id, 1))+1))
                continue;
            end
            
            % if there are multiple junctions to insert, they should beyond
            % some distanct
            if(abs(e_id-prev_e_id) < 3)
               continue; 
            end
            
            % when the intersection super close to an end point
            if(e_id<3 )
%                 % simply extend the approaching cfrag to the end
%                 cfrags{c_id} = [cfrags{c_id}; cur_c(1, :) ];
%                 cfrags_idx{c_id} = [ cfrags_idx{c_id} cur_c_idx(1)];
%                 junction_pts = [junction_pts; cur_c(1, :)];
                continue;
            end
            if(e_id > size(cur_c,1)-2)
%                 % simply extend the approaching cfrag to the end
%                 cfrags{c_id} = [cfrags{c_id}; cur_c(end, :) ];
%                 cfrags_idx{c_id} = [ cfrags_idx{c_id} cur_c_idx(end)];
%                 junction_pts = [junction_pts; cur_c(end, :)];
                continue;
            end

            % Break existing T-junction
            if(min_dist == 0)
                junction_pts = [junction_pts; cur_c(e_id, :)];
%                 introduced_num_junction_points = introduced_num_junction_points + 1; 
                cur_break_e_ids =  [cur_break_e_ids e_id];
                continue;
            end

            % only consider junction when abs(cos(ori_diff)) < cos(pi/6);
            if(e_id==1)
                c_ori_vec = cur_c(e_id + 1, 1:2) - cur_c(e_id, 1:2);
            elseif (e_id==size(cur_c,1))
                c_ori_vec = cur_c(e_id, 1:2) - cur_c(e_id-1, 1:2);  
            else
                c_ori_vec = cur_c(e_id + 1, 1:2) - cur_c(e_id-1, 1:2);
            end
            a_ori_vec = cfrags{c_id}(end, 1:2) - cfrags{c_id}(max(1,end-4), 1:2);
            cos_ori_diff = abs(c_ori_vec*a_ori_vec'/norm(c_ori_vec)/norm(a_ori_vec));
            if(cos_ori_diff > cos(pi/6))% || cos_ori_diff < cos(pi/3))
                continue;
            end

            %%%%% extend the approaching cfrag to the intersection
            cfrags{c_id} = [cfrags{c_id}; cur_c(e_id, :) ];
            cfrags_idx{c_id} = [ cfrags_idx{c_id} cur_c_idx(e_id)];
            
            junction_pts = [junction_pts; cur_c(e_id, :)];
%             introduced_num_junction_points = introduced_num_junction_points + 1; 
            cur_break_e_ids =  [cur_break_e_ids e_id];
            prev_e_id = e_id;
            % update ref maps
            ref_table_end_pts(ref_table_end_pts==c_id) = 0;
        end


        %%%%%%%%%%%%%%%%%  break the cfrag at junctions (one or more)
        %%%%%%%%%%%%% iteratively update at cfrags, finally replace new
        if(~isempty(cur_break_e_ids))
            cur_break_e_ids = sort(cur_break_e_ids, 'ascend');

            cfrags{i} = cur_c(cur_break_e_ids(end):end, :);
            cfrags = [cfrags cur_c(1:cur_break_e_ids(1), :)];

            cfrags_idx{i} = cur_c_idx(cur_break_e_ids(end):end);
            cfrags_idx = [cfrags_idx cur_c_idx(1:cur_break_e_ids(1))];
            
%             %%%%%%%  update the ref table of start pts
            ref_table_start_pts(ref_table_start_pts==i) = length(cfrags);
        
        
            %%%%%%%% break at multiple pts
            prev_e_id = cur_break_e_ids(1);
            if(length(cur_break_e_ids)>=2)
                
                for j = 1:length(cur_break_e_ids)-1                    
                    cfrags = [cfrags cur_c(prev_e_id:cur_break_e_ids(j+1), :)];
                    cfrags_idx = [cfrags_idx cur_c_idx(prev_e_id:cur_break_e_ids(j+1))];
                    
                    prev_e_id = cur_break_e_ids(j+1);
                end
            end
        end


    end

end

new_cfrags = cfrags;
new_cfrags_idx = cfrags_idx;
% introduced_num_junction_points    


end