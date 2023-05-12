function [new_cfrags, new_cfrags_idx] = contour_fill_gaps_delaunay(cfrags, cfrags_idx, h, w, edges, edgemap_soft)
bry_conf_th = 0.05;
bry_conf_ceil = 0.5;
ori_diff_th = pi/6;
len_th1 = 10;
len_th2 = 3;

disp('complete large gaps');

%% preparation of reference maps
G_test = construct_fac_graph_from_curve_fragments (cfrags_idx, cfrags);
%%%%%%%%% construct contour id map
EdgeGroupMap = convert_cfrags_to_EdgeGroupMap (cfrags, h, w);
EdgeMap = single(EdgeGroupMap>0);
EdgeMap =  conv2(EdgeMap,ones(2,2),'same');
EdgeMap = EdgeMap>0;
% figure(1);
% imshow(EdgeGroupMap, 'border', 'tight'); 
% hold on;
%%%%%%%%% construct free contour end point map, and correponding orientations
end_point_map = zeros(h,w);
end_dx_map = zeros(h,w);
end_dy_map = zeros(h,w);
jct_map = zeros(h,w);

end_eids = [G_test.var(:).actual_edge_id];
% only conside free end in dense triangulation case now
to_remove = zeros(size(end_eids));
% angle_list = zeros(size(end_eids));
end_dx_list = zeros(size(end_eids));
end_dy_list = zeros(size(end_eids));

for vid =1:length(to_remove)
    if(length(G_test.var(vid).nbrs_fac)==2)
        to_remove(vid)=1;
    elseif(length(G_test.var(vid).nbrs_fac)>2)
        to_remove(vid)=1; 
        cur_x = edges(G_test.var(vid).actual_edge_id,1);
        cur_y = edges(G_test.var(vid).actual_edge_id,2);
        jct_map(round(cur_y)+1, round(cur_x)+1) = 1;
    else
        % compute the orientation at free end
        cid = G_test.var(vid).nbrs_fac;
        cur_c = cfrags{cid};
        e_size = size(cur_c,1);
        
        if(G_test.fac(cid).nbrs_var(1)==vid)
            end_dx_list(vid) =  cur_c(1,1) - cur_c(min(e_size,3),1);
            end_dy_list(vid) =  cur_c(1,2) - cur_c(min(e_size,3),2);
        else
            end_dx_list(vid) =  cur_c(end,1) - cur_c(max(e_size-2, 1),1);
            end_dy_list(vid) =  cur_c(end,2) - cur_c(max(e_size-2, 1),2);
        end

    end
end
end_eids(to_remove>0) = [];
end_dx_list(to_remove>0) = [];
end_dy_list(to_remove>0) = [];
end_Y = round(edges(end_eids', 2))+1;
end_X = round(edges(end_eids', 1))+1;
end_point_map(sub2ind([h,w], end_Y, end_X))=1;
end_dx_map(sub2ind([h,w], end_Y, end_X))=end_dx_list;
end_dy_map(sub2ind([h,w], end_Y, end_X))=end_dy_list;

%%%%%%%%%% construct junction map
EDGEIM = EdgeGroupMap ~= 0;           % make sure image is binary.
EDGEIM = bwmorph(EDGEIM,'clean');     % Remove isolated pixels
EDGEIM = bwmorph(EDGEIM,'skel',Inf);  % and make sure edges are thinned. I
                                      % think using 'skel' is better than 'thin'    
% Find endings and junctions in edge data

[Jct_Y, Jct_X, re, ce] = findendsjunctions(EDGEIM);
jct_map(sub2ind([h,w], Jct_Y, Jct_X))=1;
jct_map = conv2(jct_map,ones(5,5),'same');
%% Propose gap completion between contour end points, based on DT (delaunay) betwen end points

vertices = [G_test.var(:).actual_edge_id];
Y = round(edges(vertices', 2))+1;
X = round(edges(vertices', 1))+1;
is_connect_to_jct = jct_map(sub2ind([h,w], Y, X));
Y(is_connect_to_jct>0)=[];
X(is_connect_to_jct>0)=[];


constraints = [G_test.fac(:).nbrs_var];
constraints = reshape(constraints', [2, length(constraints)/2]);
constraints = constraints';
constraints = sort(constraints,2);

tic
% tri = delaunayTriangulation(X,Y, constraints); 
tri = delaunayTriangulation(X,Y); 
toc
c_links_ids_t1 = tri.edges;
c_links_ids_t1 = sort(c_links_ids_t1, 2);

% [~,index_c,index_r] = intersect(c_links_ids_t1,constraints,'rows');
% c_links_ids_t1(index_c,:)= [];

% %%%%%%%% remove links connected to T-junctions (Y-junction case is ruled out already)
% links_start = [X(c_links_ids_t1(:,1)) Y(c_links_ids_t1(:,1))];
% links_end = [X(c_links_ids_t1(:,2)) Y(c_links_ids_t1(:,2))];
% is_connect_to_jct = jct_map(sub2ind([h,w], links_start(:,2), links_start(:,1)));
% is_connect_to_jct = is_connect_to_jct | jct_map(sub2ind([h,w], links_end(:,2), links_end(:,1)));
% c_links_ids_t1(is_connect_to_jct,:) = [];

%%%%%%%% remove links cross exisiting contours
% links_start = [X(c_links_ids_t1(:,1)) Y(c_links_ids_t1(:,1))];
% links_end = [X(c_links_ids_t1(:,2)) Y(c_links_ids_t1(:,2))];
% links_len = sqrt(sum((links_end-links_start).^2, 2));
to_remove = zeros(size(c_links_ids_t1,1),1);
for i=1:length(to_remove)
    
    % check the link's average confidence
%    [CX, CY, conf_vals] = improfile(EdgeGroupMap>0, X(c_links_ids_t1(i,:)), Y(c_links_ids_t1(i,:)));
    conf_vals = improfile(EdgeMap, X(c_links_ids_t1(i,:)), Y(c_links_ids_t1(i,:)));

   if (sum(conf_vals(2:end-1))>1)
        to_remove(i) = 1;
%    elseif(links_len(i)>50)
%        plot (CX, CY, 'r.');
%        keyboard; 
   end
    
end
c_links_ids_t1(to_remove>0,:) = [];


%%%%%%% keep end-to-end links with good coninuity at both side
links_start = [X(c_links_ids_t1(:,1)) Y(c_links_ids_t1(:,1))];
links_end = [X(c_links_ids_t1(:,2)) Y(c_links_ids_t1(:,2))];
end_dx = end_dx_map(sub2ind([h,w], links_start(:,2), links_start(:,1)));
end_dy = end_dy_map(sub2ind([h,w], links_start(:,2), links_start(:,1)));
% end_dx_at_end = end_dx_map(sub2ind([h,w], links_end(:,2), links_end(:,1)));
% end_dy_at_end = end_dy_map(sub2ind([h,w], links_end(:,2), links_end(:,1)));
dx = links_end(:,1) - links_start(:,1);
dy = links_end(:,2) - links_start(:,2);
% judge at starting side
is_connect_to_end_1 = end_point_map(sub2ind([h,w], links_start(:,2), links_start(:,1)));
is_connect_to_end_1 = is_connect_to_end_1>0;
cos_ori_diff_1 = (end_dx.*dx + end_dy.*dy)./sqrt(dx.*dx + dy.*dy)./sqrt(end_dx.*end_dx + end_dy.*end_dy);
id_to_keep = cos_ori_diff_1 > cos(ori_diff_th);
id_to_keep = id_to_keep & is_connect_to_end_1;

% judge at ending side
is_connect_to_end_2 = end_point_map(sub2ind([h,w], links_end(:,2), links_end(:,1)));
is_connect_to_end_2 = is_connect_to_end_2>0;
dx(is_connect_to_end_2) = -dx(is_connect_to_end_2);
dy(is_connect_to_end_2) = -dy(is_connect_to_end_2);
end_dx = end_dx_map(sub2ind([h,w], links_end(:,2), links_end(:,1)));
end_dy = end_dy_map(sub2ind([h,w], links_end(:,2), links_end(:,1)));
cos_ori_diff_2 = (end_dx.*dx + end_dy.*dy)./sqrt(dx.*dx + dy.*dy)./sqrt(end_dx.*end_dx + end_dy.*end_dy);
id_to_keep = id_to_keep & cos_ori_diff_2 > cos(ori_diff_th);
id_to_keep = id_to_keep & is_connect_to_end_2;

links_len = sqrt(sum((links_end-links_start).^2, 2));
id_to_keep = id_to_keep & (links_len < len_th1);
id_to_keep = id_to_keep | (links_len < len_th2);

%%%%%%% select the most continuous link at each vertex, keep the best if enough bry_conf
% sort by starting index
[c_links_ids_t1, sort_order] =  sortrows(c_links_ids_t1);
id_to_keep = id_to_keep(sort_order);
is_connect_to_end_1 = is_connect_to_end_1(sort_order);
is_connect_to_end_2 = is_connect_to_end_2(sort_order);
cos_ori_diff_1 = cos_ori_diff_1(sort_order);
cos_ori_diff_2 = cos_ori_diff_2(sort_order);

max_continue = -1;
max_i = 0;
prev_end_idx = 0;
for i=1:length(id_to_keep)
    
   if(~is_connect_to_end_1(i))
      continue; 
   end
    
   if(c_links_ids_t1(i,1) == prev_end_idx)
       % only keep the besk link per contour end point
        if(id_to_keep(i)==1) % already found previously enforced link
           max_continue = 2; % asign a value no others can excend 
           max_i = i;
        elseif(cos_ori_diff_1(i)>max_continue)
           max_continue = cos_ori_diff_1(i);
           max_i = i;
       end
   else
       % check the link's average confidence
       if(max_i ~= 0 && max_continue <=1) % exclude the case been enforced
            conf_vals = improfile(edgemap_soft, X(c_links_ids_t1(max_i,:)), Y(c_links_ids_t1(max_i,:)));
            conf_vals(conf_vals>bry_conf_ceil) = bry_conf_th; % aiming for weak curve, take very strong bry as noise
            conf_val = mean(conf_vals(2:end-1));
       
           if(conf_val > bry_conf_th && max_continue > cos(ori_diff_th))
              id_to_keep(max_i)=1; 
           end
       end
       
       % intilize max_conf again
       max_continue = cos_ori_diff_1(i);
       max_i = i;
       prev_end_idx = c_links_ids_t1(i,1);
       
   end
    
end

% sort by ending index
[c_links_ids_t1, sort_order] =  sortrows(c_links_ids_t1, 2);
id_to_keep = id_to_keep(sort_order);
is_connect_to_end_2 = is_connect_to_end_2(sort_order);
cos_ori_diff_2 = cos_ori_diff_2(sort_order);

max_continue = -1;
max_i = 0;
prev_end_idx = 0;
for i=1:length(id_to_keep)
    
   if(~is_connect_to_end_2(i))
      continue; 
   end
    
   if(c_links_ids_t1(i,1) == prev_end_idx)
       % only keep the besk link per contour end point
        if(id_to_keep(i)==1) % already found previously enforced link
           max_continue = 2; % asign a value no others can excend 
           max_i = i;
        elseif(cos_ori_diff_2(i)>max_continue)
           max_continue = cos_ori_diff_2(i);
           max_i = i;
       end
   else
       % check the link's average confidence
       if(max_i ~= 0 && max_continue <=1) % exclude the case been enforced
            conf_vals = improfile(edgemap_soft, X(c_links_ids_t1(max_i,:)), Y(c_links_ids_t1(max_i,:)));
            conf_vals(conf_vals>bry_conf_ceil) = bry_conf_th; % aiming for weak curve, take very strong bry as noise
            conf_val = mean(conf_vals(2:end-1));
       
           if(conf_val > bry_conf_th && max_continue > cos(ori_diff_th))
              id_to_keep(max_i)=1; 
           end
       end
       
       % intilize max_conf again
       max_continue = cos_ori_diff_2(i);
       max_i = i;
       prev_end_idx = c_links_ids_t1(i,1);
       
   end
    
end
c_links_ids_t1 = c_links_ids_t1(id_to_keep, :);


% figure(1);
% imshow(edgemap_soft, 'border', 'tight'); 
% hold on;
draw_contours(cfrags,0,0);
for cid=1:size(c_links_ids_t1,1)
    plot(X(c_links_ids_t1(cid,:)), Y(c_links_ids_t1(cid,:)), '-r');
end

%% Propose T-junction based on DT (delaunay) on all image edges
[Y_all, X_all] = find(EdgeGroupMap>0);
is_connect_to_jct = jct_map(sub2ind([h,w], Y_all, X_all));
Y_all(is_connect_to_jct>0)=[];
X_all(is_connect_to_jct>0)=[];
tri = delaunayTriangulation(X_all,Y_all);
c_links_ids_t2 = tri.edges;

%%%%%%%%% remove links sharing a single contour id
Gid = EdgeGroupMap(sub2ind([h,w], Y_all, X_all));
c_links_ids_t2(Gid(c_links_ids_t2(:,1)) == Gid(c_links_ids_t2(:,2)),:) = [];
%TODO: remove some remaining neibouring links

%%%%%%%%% remove links do not connected to an end point point
links_start = [X_all(c_links_ids_t2(:,1)) Y_all(c_links_ids_t2(:,1))];
links_end = [X_all(c_links_ids_t2(:,2)) Y_all(c_links_ids_t2(:,2))];
is_connect_to_end_1 = end_point_map(sub2ind([h,w], links_start(:,2), links_start(:,1)));
is_connect_to_both_end = is_connect_to_end_1 & end_point_map(sub2ind([h,w], links_end(:,2), links_end(:,1)));
is_connect_to_end_point = is_connect_to_end_1 | end_point_map(sub2ind([h,w], links_end(:,2), links_end(:,1)));
% remove the case of tangential gaps
is_connect_to_end_point = is_connect_to_end_point & (~is_connect_to_both_end);
c_links_ids_t2 = c_links_ids_t2(is_connect_to_end_point,:);

% %%%%%%%% remove links connected to T-junctions (Y-junction case is ruled out already)
% links_start = [X_all(c_links_ids_t2(:,1)) Y_all(c_links_ids_t2(:,1))];
% links_end = [X_all(c_links_ids_t2(:,2)) Y_all(c_links_ids_t2(:,2))];
% is_connect_to_jct = jct_map(sub2ind([h,w], links_start(:,2), links_start(:,1)));
% is_connect_to_jct = is_connect_to_jct | jct_map(sub2ind([h,w], links_end(:,2), links_end(:,1)));
% c_links_ids_t2(is_connect_to_jct,:) = [];

%%%%%%% simple remove those having super weak boundary conf value in middle
links_start = [X_all(c_links_ids_t2(:,1)) Y_all(c_links_ids_t2(:,1))];
links_end = [X_all(c_links_ids_t2(:,2)) Y_all(c_links_ids_t2(:,2))];
links_mid = round((links_start+links_end)/2);
bry_val = edgemap_soft(sub2ind([h,w], links_mid(:,2), links_mid(:,1)));
c_links_ids_t2(bry_val<bry_conf_th, :) = [];

%%%%%%% remove those < 3 pixels in length
links_start = [X_all(c_links_ids_t2(:,1)) Y_all(c_links_ids_t2(:,1))];
links_end = [X_all(c_links_ids_t2(:,2)) Y_all(c_links_ids_t2(:,2))];
links_len = sqrt(sum((links_end-links_start).^2, 2));
c_links_ids_t2(links_len<len_th2, :) = [];

%%%%%%% TODO: remove links between a corner and the middle of a contour 


%%%%%%% remove links not continuous with the stem endpoint
links_start = [X_all(c_links_ids_t2(:,1)) Y_all(c_links_ids_t2(:,1))];
links_end = [X_all(c_links_ids_t2(:,2)) Y_all(c_links_ids_t2(:,2))];
end_dx = end_dx_map(sub2ind([h,w], links_start(:,2), links_start(:,1)));
end_dy = end_dy_map(sub2ind([h,w], links_start(:,2), links_start(:,1)));
end_dx_at_end = end_dx_map(sub2ind([h,w], links_end(:,2), links_end(:,1)));
end_dy_at_end = end_dy_map(sub2ind([h,w], links_end(:,2), links_end(:,1)));
dx = links_end(:,1) - links_start(:,1);
dy = links_end(:,2) - links_start(:,2);
% assume only one contour end point
is_connect_to_end_1 = end_point_map(sub2ind([h,w], links_start(:,2), links_start(:,1)));
dx(~is_connect_to_end_1) = -dx(~is_connect_to_end_1);
dy(~is_connect_to_end_1) = -dy(~is_connect_to_end_1);
end_dx(~is_connect_to_end_1) = end_dx_at_end(~is_connect_to_end_1);
end_dy(~is_connect_to_end_1) = end_dy_at_end(~is_connect_to_end_1);
cos_ori_diff = (end_dx.*dx + end_dy.*dy)./sqrt(dx.*dx + dy.*dy)./sqrt(end_dx.*end_dx + end_dy.*end_dy);
c_links_ids_t2(cos_ori_diff < cos(ori_diff_th), :) = [];

%%%%%%% For each contour end point, keep the best link in avg boundary conf and orientation agreement
% Only one side is of a link contour end point, make it the starting side
links_end = [X_all(c_links_ids_t2(:,2)) Y_all(c_links_ids_t2(:,2))];
is_connect_to_end2 = end_point_map(sub2ind([h,w], links_end(:,2), links_end(:,1)));
c_links_ids_t2(is_connect_to_end2>0,:) = [c_links_ids_t2(is_connect_to_end2>0,2), c_links_ids_t2(is_connect_to_end2>0,1)];
% sort by starting index
c_links_ids_t2 =  sortrows(c_links_ids_t2);
max_conf = 0;
max_i = 0;
prev_end_idx = 0;
to_remove = zeros(size(c_links_ids_t2,1),1);
for i=1:length(to_remove)
    
    % check the link's average confidence
   conf_vals = improfile(edgemap_soft, X_all(c_links_ids_t2(i,:)), Y_all(c_links_ids_t2(i,:)));
   conf_vals(conf_vals>bry_conf_ceil) = bry_conf_th*2; % aiming for weak curve, take very strong bry as noise
   conf_val = mean(conf_vals(2:end-1));
   if(c_links_ids_t2(i,1) == prev_end_idx)
       % only keep the max conf link per contour end point
       if(conf_val>max_conf)
           max_conf = conf_val;
           to_remove(max_i) = 1;
           max_i = i;
       else
           to_remove(i) = 1;
       end
   else
       if(max_conf<bry_conf_th && max_i~=0)
          to_remove(max_i)=1; 
       end
       
       % intilize max_conf again
       max_conf = conf_val;
       max_i = i;
       prev_end_idx = c_links_ids_t2(i,1);
       
   end
    
end
c_links_ids_t2(to_remove>0,:) = [];


for cid=1:size(c_links_ids_t2,1)
    plot(X_all(c_links_ids_t2(cid,:)), Y_all(c_links_ids_t2(cid,:)), '-b');
end
hold off;

new_cfrags = cfrags;
new_cfrags_idx = cfrags_idx;
% only keep the non-empty contours
new_cfrags = new_cfrags(~cellfun('isempty',new_cfrags));
new_cfrags_idx = new_cfrags_idx(~cellfun('isempty',new_cfrags_idx)); 


end