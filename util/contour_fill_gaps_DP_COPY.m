function [new_cfrags, new_cfrags_idx] = contour_fill_gaps_DP(cfrags, cfrags_idx, h, w, edges, edgemap_soft, edgemap, thetamap)
% it maybe better to run DP under increasing range radius, small range DP
% will very efficient. When the small gaps are filled, the remaining number
% of free ends are very small
r_range = 50;
cost_th = 1;
%% preparation of reference maps
G_test = construct_fac_graph_from_curve_fragments (cfrags_idx, cfrags);

% [Ox,Oy]=gradient2(convTri(single(edgemap_soft),4));
% [Oxx,~]=gradient2(Ox); [Oxy,Oyy]=gradient2(Oy);
% O_map=mod(atan(Oyy.*sign(-Oxy)./(Oxx+1e-5)),pi);
% O_map = pi*rand([h,w]);
% O_map(edgemap>0)=thetamap(edgemap>0);
O_map=thetamap;

%%%%%%%%% construct contour id map
EdgeGroupMap = convert_cfrags_to_EdgeGroupMap (cfrags, h, w);

%%%%%%%%% construct free contour end point map, and correponding orientations
end_point_map = zeros(h,w);
end_point_map_org = zeros(h,w);
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
            end_dx_list(vid) =  cur_c(1,1) - cur_c(min(e_size,4),1);
            end_dy_list(vid) =  cur_c(1,2) - cur_c(min(e_size,4),2);
        else
            end_dx_list(vid) =  cur_c(end,1) - cur_c(max(e_size-3, 1),1);
            end_dy_list(vid) =  cur_c(end,2) - cur_c(max(e_size-3, 1),2);
        end

    end
end
end_Y = round(edges(end_eids', 2))+1;
end_X = round(edges(end_eids', 1))+1;
end_point_map_org(sub2ind([h,w], end_Y, end_X))=1; % this map include all the breaking points

end_eids(to_remove>0) = [];
end_dx_list(to_remove>0) = [];
end_dy_list(to_remove>0) = [];
end_Y = round(edges(end_eids', 2))+1;
end_X = round(edges(end_eids', 1))+1;
end_point_map(sub2ind([h,w], end_Y, end_X))=1;
end_dx_map(sub2ind([h,w], end_Y, end_X))=end_dx_list;
end_dy_map(sub2ind([h,w], end_Y, end_X))=end_dy_list;
end_angle_map = atan2(end_dy_map, end_dx_map);
% assigin end orientation to O_map
% O_map(sub2ind([h,w], end_Y, end_X)) = wrapToPi(end_angle_map(sub2ind([h,w], end_Y, end_X))+pi);
% make thick edge map
% DT of binary map
DT_map = bwdist(EdgeGroupMap>0,'euclidean');
DT_map = double(exp(-DT_map.^2));

%%%%%%%%%% construct junction map (This can be removed is T-junction is known)
EDGEIM = EdgeGroupMap ~= 0;           % make sure image is binary.
EDGEIM = bwmorph(EDGEIM,'clean');     % Remove isolated pixels
EDGEIM = bwmorph(EDGEIM,'skel',Inf);  % and make sure edges are thinned. I
                                      % think using 'skel' is better than 'thin'    
% Find endings and junctions in edge data
[Jct_Y, Jct_X, re, ce] = findendsjunctions(EDGEIM);
jct_map(sub2ind([h,w], Jct_Y, Jct_X))=1;
jct_map = conv2(jct_map,ones(3,3),'same');
% remove those T-junctions in free end map
[Jct_Y, Jct_X] = find(jct_map>0);
end_point_map(sub2ind([h,w], Jct_Y, Jct_X))=0;


%% Dynamic programming at each end point and find the optimal path within a range
[Y,X] =  find(end_point_map);
min_cost_all = [];
for i=1:length(X)
    
    start_point = [X(i) Y(i) end_angle_map(Y(i),X(i))];
    
    x_min = max(X(i)-r_range,1);
    x_max = min(X(i)+r_range,w);
    y_min = max(Y(i)-r_range,1);
    y_max = min(Y(i)+r_range,h);    
    
    start_point_ref = start_point; % this is in c++ coordinates
    start_point_ref(1) = start_point_ref(1) - x_min;
    start_point_ref(2) = start_point_ref(2) - y_min;
    h_ref = y_max-y_min+1;
    w_ref = x_max-x_min+1;
    
    EdgeGroupMap_ref = EdgeGroupMap(y_min:y_max, x_min:x_max);
    end_point_map_org_ref = end_point_map_org(y_min:y_max, x_min:x_max);
    [ex_con_pointY, ex_con_pointX] = find(EdgeGroupMap_ref>0 & EdgeGroupMap_ref~= EdgeGroupMap_ref(start_point_ref(2)+1, start_point_ref(1)+1));
 
    %%%%%% A: rule out those not dontinuous with starting contour
    dx = ex_con_pointX - start_point_ref(1)-1;
    dy = ex_con_pointY - start_point_ref(2)-1;
    % More restricted angle may needed for T-junction type gap, The
    % searching range limit should be integrated into DP to accelerate the
    % process
    is_end_point = end_point_map_org_ref(sub2ind([h_ref, w_ref], ex_con_pointY, ex_con_pointX));
    cos_diff = (cos(start_point(3)) * dx + sin(start_point(3)) * dy)./ sqrt(dx.*dx + dy.*dy);
    link_to_remove = (is_end_point==0 & cos_diff<=0.95); % strict rule for T-junction
    link_to_remove = link_to_remove | (is_end_point>0 & cos_diff<=0.5); % allow some turing in end-end gap
    
    ex_con_pointX(link_to_remove)=[];
    ex_con_pointY(link_to_remove)=[];       
    is_end_point(link_to_remove)=[];
    if(isempty(ex_con_pointX))
        continue;
    end    
    
%     edgemap_ref = double(edgemap(y_min:y_max, x_min:x_max)>0);
%      edgemap_ref = double(EdgeGroupMap_ref>0);
%     edgemap_ref(EdgeGroupMap_ref== EdgeGroupMap_ref(start_point_ref(2)+1, start_point_ref(1)+1)) = 0;
    edgemap_ref = DT_map(y_min:y_max, x_min:x_max);

    %% make sure all the matrix input in double format
    [cost_mat, back_p_mat, len_mat] = Live_Wire_2D_DP_gap_cpt(edgemap_ref, edgemap_soft(y_min:y_max, x_min:x_max), ...
                                            O_map(y_min:y_max, x_min:x_max), start_point_ref, h_ref, w_ref);
    
           
    
                                        
    %%%%%% B: rule out those not with avg high cost
% %     cost_mat = cost_mat./len_mat; 
    end_costs = cost_mat(sub2ind([h_ref, w_ref], ex_con_pointY, ex_con_pointX));
    %%%%%% B2: add preference for end-to-end completion
%     end_costs(is_end_point>0) = end_costs(is_end_point>0)*0.8;
    end_costs_norm = end_costs ./ len_mat(sub2ind([h_ref, w_ref], ex_con_pointY, ex_con_pointX));
    end_costs(end_costs_norm>cost_th)=1000;
    
    [min_cost, min_ind] = min(end_costs);
    if(min_cost==1000)
        continue;
    end
    if(end_costs_norm(min_ind) > cost_th/2)
        continue;
    end
    
    % back track to get the optimal path
    end_point_x = ex_con_pointX(min_ind);
    end_point_y = ex_con_pointY(min_ind);
    
    cur_ind = sub2ind([h_ref, w_ref], end_point_y,end_point_x);
    opt_path = cur_ind;
    while(back_p_mat(cur_ind))
        opt_path = [opt_path back_p_mat(cur_ind)+1];
        cur_ind = back_p_mat(cur_ind)+1; % covert to matlab coordinates
    end
    [opt_path_y, opt_path_x] = ind2sub([h_ref,w_ref], opt_path);
    opt_path_y = opt_path_y+y_min-1;
    opt_path_x = opt_path_x+x_min-1;
    opt_path_y = [opt_path_y start_point(2)];
    opt_path_x = [opt_path_x start_point(1)];
    
    
%     %%%%%% C: rule out those gap paths that are longer than its connected contour and are very low avg contrast
%     path_prob = edgemap_soft(sub2ind([h,w], opt_path_y(2:end-1), opt_path_x(2:end-1)));
%     if((mean(path_prob)<0.05 && length(path_prob)> 10)) %|| (mean(path_prob)<0.05 && length(path_prob)<= 5))
%         continue;
%     end
    
    %%%%%% HOW to solve the remaining confict?


    plot(opt_path_x, opt_path_y, 'r-');
    
    min_cost_all = [min_cost_all min_cost];

end
%  keyboard;
new_cfrags = cfrags;
new_cfrags_idx = cfrags_idx;

