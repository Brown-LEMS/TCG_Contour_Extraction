function [new_cfrags, new_cfrags_idx, corner_pts, junction_pts] = contour_breaker_geom_v2(cfrags, cfrags_idx, edgemap, params)
% locate the conners and break the contour
b_extend = 1; % locate the junctions: fill the gap and break the contour
ref_tabel_nbr_range = 2;

diag_ratio= params.diag_ratio;
merge_th_geom = params.merge_th_geom; 
nbr_num_edges = params.nbr_num_edges;
beta_1 = params.beta_1;
fmean_1 = params.fmean_1;
% max_iter = params.max_iter;
max_iter = 0; % do not break

disp('breaking contours at conners and form junctions');
[h,w] = size(edgemap);


%% look for junctions, break the cfrag, extent the approaching cfrag
introduced_num_junction_points = 0;
junction_pts = [];

new_cfrags = cfrags;
num_org_cfrags = length(new_cfrags);

new_cfrags_idx = cfrags_idx;

%%%%%%%%%%%%%%%%%%%%% compute table refering end pts %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% extract start/end pts:  format: x, y, theta
start_pts = zeros(num_org_cfrags, 3);
end_pts = zeros(num_org_cfrags, 3);
c_lens = zeros(num_org_cfrags, 1);
ref_table_start_pts = zeros(size(edgemap));
ref_table_end_pts = zeros(size(edgemap));


for i = 1: num_org_cfrags
    cur_c = new_cfrags{i};
    c_len = contour_length_mex(cur_c');
    c_lens(i) = c_len;
    start_pts(i,:) = new_cfrags{i}(1,1:3);
    end_pts(i,:) = new_cfrags{i}(end,1:3);

    
    %%%%%% skip circles
    if(cur_c(1,:) == cur_c(end,:))
        continue;
    end
    
    %%%%%% only consider length>10 contour as approaching cfrags
    %%%%%%% TODO: should avoid degree 2 nodes
    if(c_len > 5*diag_ratio && size(cur_c, 1) >5)
        x_coords_start = round(start_pts(i,1)+1)-ref_tabel_nbr_range : round(start_pts(i,1)+1)+ref_tabel_nbr_range;
        y_coords_start = round(start_pts(i,2)+1)-ref_tabel_nbr_range : round(start_pts(i,2)+1)+ref_tabel_nbr_range;
        x_coords_start = max(x_coords_start, ones(size(x_coords_start)));
        y_coords_start = max(y_coords_start, ones(size(y_coords_start)));
        x_coords_start = min(x_coords_start, w*ones(size(x_coords_start)));
        y_coords_start = min(y_coords_start, h*ones(size(y_coords_start)));

        ref_table_start_pts(y_coords_start, x_coords_start) = i;

        x_coords_end = round(end_pts(i,1)+1)-ref_tabel_nbr_range : round(end_pts(i,1)+1)+ref_tabel_nbr_range;
        y_coords_end = round(end_pts(i,2)+1)-ref_tabel_nbr_range : round(end_pts(i,2)+1)+ref_tabel_nbr_range;
        x_coords_end = max(x_coords_end, ones(size(x_coords_end)));
        y_coords_end = max(y_coords_end, ones(size(y_coords_end)));
        x_coords_end = min(x_coords_end, w*ones(size(x_coords_end)));
        y_coords_end = min(y_coords_end, h*ones(size(y_coords_end)));

        ref_table_end_pts(y_coords_end, x_coords_end) = i;  
    end
end

%%%%%%%%%%%%%%%%%%%%%  look for junctions, break the cfrag, extent the approaching cfrag
for i = 1: num_org_cfrags
    cur_c = cfrags{i};
    c_len = c_lens(i);

    cur_c_idx = new_cfrags_idx{i};
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

            
            % skip the pts close to ends
            if(e_id<5 || e_id > size(cur_c,1)-4 || (e_id-prev_e_id) < 5)
                continue;
            end

            % if the extend edge is already junction, just to break
            if(min_dist == 0)
                junction_pts = [junction_pts; cur_c(e_id, :)];
                introduced_num_junction_points = introduced_num_junction_points + 1; 
                cur_break_e_ids =  [cur_break_e_ids e_id];
                continue;
            end


            % only consider Y junction when abs(cos(ori_diff)) in [cos(pi/3), cos(pi/6) ];
            c_ori_vec = cur_c(e_id + 1, 1:2) - cur_c(e_id-1, 1:2);

%             try
            a_ori_vec = cfrags{c_id}(5, 1:2) - cfrags{c_id}(1, 1:2);
%             catch
%                 keyboard
%             end
            cos_ori_diff = abs(c_ori_vec*a_ori_vec'/norm(c_ori_vec)/norm(a_ori_vec));
            if( cos_ori_diff > cos(pi/6))% || cos_ori_diff < cos(pi/3))
                continue;
            end


            %%%%% extend the approaching cfrag at the start                
            junction_pts = [junction_pts; cur_c(e_id, :)];
            introduced_num_junction_points = introduced_num_junction_points + 1; 
            cur_break_e_ids =  [cur_break_e_ids e_id];

            if(b_extend)
                cfrags{c_id} = [cur_c(e_id, :);  cfrags{c_id}];
                new_cfrags_idx{c_id} = [cur_c_idx(e_id) new_cfrags_idx{c_id}];
            end
            
%             if(size(cfrags{c_id},1)~=length(new_cfrags_idx{c_id}))
%                 keyboard;
%             end
            
            prev_e_id = e_id;
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


            
            % skip the pts close to ends
            if(e_id<5 || e_id > size(cur_c,1)-4 || (e_id-prev_e_id) < 5)
                continue;
            end

            % if the extend edge is already junction, just to break
            if(min_dist == 0)
                junction_pts = [junction_pts; cur_c(e_id, :)];
                introduced_num_junction_points = introduced_num_junction_points + 1; 
                cur_break_e_ids =  [cur_break_e_ids e_id];
                continue;
            end

            % only consider Y junction when abs(cos(ori_diff)) < cos(pi/4);
            c_ori_vec = cur_c(e_id + 1, 1:2) - cur_c(e_id-1, 1:2);
            a_ori_vec = cfrags{c_id}(end, 1:2) - cfrags{c_id}(end-4, 1:2);
            cos_ori_diff = abs(c_ori_vec*a_ori_vec'/norm(c_ori_vec)/norm(a_ori_vec));
            if(cos_ori_diff > cos(pi/6))% || cos_ori_diff < cos(pi/3))
                continue;
            end

            %%%%% extend the approaching cfrag at the end
            junction_pts = [junction_pts; cur_c(e_id, :)];
            introduced_num_junction_points = introduced_num_junction_points + 1; 
            cur_break_e_ids =  [cur_break_e_ids e_id];

            if(b_extend)
                cfrags{c_id} = [cfrags{c_id}; cur_c(e_id, :) ];
                new_cfrags_idx{c_id} = [ new_cfrags_idx{c_id} cur_c_idx(e_id)];
            end
%             if(size(cfrags{c_id},1)~=length(new_cfrags_idx{c_id}))
%                 keyboard;
%             end
            
            prev_e_id = e_id;

        end


        %%%%%%%%%%%%%%%%%  break the cfrag at junctions (one or more)
        %%%%%%%%%%%%% iteratively update at cfrags, finally replace new
        if(~isempty(cur_break_e_ids))
            cur_break_e_ids = sort(cur_break_e_ids, 'ascend');

            cfrags{i} = cur_c(cur_break_e_ids(end):end, :);
            cfrags = [cfrags cur_c(1:cur_break_e_ids(1), :)];

            new_cfrags_idx{i} = cur_c_idx(cur_break_e_ids(end):end);
            new_cfrags_idx = [new_cfrags_idx cur_c_idx(1:cur_break_e_ids(1))];
            
            
            %%%%%%%  update the ref table of start pts
            x_coords_start = round(cur_c(1,1)+1)-ref_tabel_nbr_range : round(cur_c(1,1)+1)+ref_tabel_nbr_range;
            y_coords_start = round(cur_c(1,2)+1)-ref_tabel_nbr_range : round(cur_c(1,2)+1)+ref_tabel_nbr_range;
            x_coords_start = max(x_coords_start, ones(size(x_coords_start)));
            y_coords_start = max(y_coords_start, ones(size(y_coords_start)));
            x_coords_start = min(x_coords_start, w*ones(size(x_coords_start)));
            y_coords_start = min(y_coords_start, h*ones(size(y_coords_start)));

            ref_table_start_pts(y_coords_start, x_coords_start) = length(cfrags);
        
        
            %%%%%%%% break at multiple pts
            prev_e_id = cur_break_e_ids(1);
            if(length(cur_break_e_ids)>=2)
                
                for j = 1:length(cur_break_e_ids)-1
%                     if((cur_break_e_ids(j+1) - prev_e_id) < nbr_num_edges/2)
%                         continue;
%                     end
                    
                    cfrags = [cfrags cur_c(prev_e_id:cur_break_e_ids(j+1), :)];
                    new_cfrags_idx = [new_cfrags_idx cur_c_idx(prev_e_id:cur_break_e_ids(j+1))];
                    
                    prev_e_id = cur_break_e_ids(j+1);
                end
            end
        end


    end

end

new_cfrags = cfrags;
% introduced_num_junction_points    

%% look for minimum merge prob given geometry (given results from prev)
corner_pts = [];
introduced_num_corner_points = 0;
for iter = 1:max_iter
    num_org_cfrags = length(new_cfrags);  

    for i = 1: num_org_cfrags
        cur_c = new_cfrags{i};
        c_len = contour_length_mex(cur_c');

        cur_c_idx = new_cfrags_idx{i};

        % skip closed circle
        if(cur_c(1,:) == cur_c(end,:))
            continue;
        end

        % only introducing breaking points for curves which are long enough
        if(c_len>(10 * diag_ratio) && size(cur_c,1)>(nbr_num_edges/iter))


            merge_prob_vec = filter_merge_prob_geom_along_cfrag(cur_c, ceil(nbr_num_edges/iter/2), beta_1, fmean_1);
            [min_prob, id] = min(merge_prob_vec);
            id = median(id); % in case there are sequntial id with the same min_prob, due to down sampling

            if(min_prob< merge_th_geom)
                new_cfrags = [new_cfrags cur_c(1:id, :)];
                new_cfrags{i} = cur_c(id:end, :);

                new_cfrags_idx = [new_cfrags_idx cur_c_idx(1:id)];                   
                new_cfrags_idx{i} = cur_c_idx(id:end);
                    
                corner_pts = [corner_pts; cur_c(id,:)];
                introduced_num_corner_points = introduced_num_corner_points+1;
            end

        end

    end

end

% introduced_num_corner_points

end