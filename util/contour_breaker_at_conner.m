function [new_cfrags, new_cfrags_idx, corner_pts] = contour_breaker_at_conner(cfrags, cfrags_idx, params, ori_diff_th)
% locate the conners and break the contour
if(nargin<4)
    ori_diff_th = params.corner_angle_th;
end
if(nargin<5)
    nbr_num_edges = params.nbr_num_edges;
end
% diag_ratio= params.diag_ratio;
disp('breaking contours at conners and form junctions');


%% look for junctions, break the cfrag, extent the approaching cfrag
new_cfrags = cfrags;
new_cfrags_idx = cfrags_idx;
  
%% look for minimum merge prob given geometry (given results from prev)
corner_pts = [];
introduced_num_corner_points = 0;
iter = 1;
% for iter = 2:-1:1 % multi-scale conner detection
    
    num_org_cfrags = length(new_cfrags);  
    i=1;
    while( i <= num_org_cfrags)
        cur_c = new_cfrags{i};
        c_len = contour_length_mex(cur_c');
        cur_c_idx = new_cfrags_idx{i};

        % skip small closed circle
        if(c_len < (2 * nbr_num_edges) && sum(cur_c(1,:) == cur_c(end,:))== 5)
            i = i+1;
            continue;
        end

        % only introducing breaking points for curves which are long enough
        if(c_len>(nbr_num_edges) && size(cur_c,1)>(nbr_num_edges))


%             ori_diff_vec = filter_merge_prob_geom_along_cfrag(cur_c, ceil(nbr_num_edges/iter/2), beta_1, fmean_1);
%             ori_diff_vec = filter_ori_diff_along_cfrag(cur_c, ceil(nbr_num_edges/2));
            cur_nbr_num = min(ceil(nbr_num_edges/2), ceil(size(cur_c,1)/6)) / iter;
            [ori_diff_vec, Dtheta] =  filter_co_circular_along_cfrag(cur_c, cur_nbr_num);
            
            ori_diff_vec(1:2*cur_nbr_num)=0;
            ori_diff_vec(end-2*cur_nbr_num:end)=0;
            
            ori_diff_vec = smooth(ori_diff_vec);
%             ori_diff_vec = hampel(ori_diff_vec);
            % rule out locally extremely smooth spots;
            ori_diff_vec(Dtheta<(ori_diff_th/2)) = Dtheta(Dtheta<(ori_diff_th/2));
            
            [max_ori_diff, id] = max(ori_diff_vec);
            id = median(id); % in case there are sequntial id with the same min_prob, due to down sampling
            
            if(max_ori_diff> ori_diff_th)
                new_cfrags = [new_cfrags cur_c(1:id, :)];
                new_cfrags{i} = cur_c(id:end, :);

                new_cfrags_idx = [new_cfrags_idx cur_c_idx(1:id)];                   
                new_cfrags_idx{i} = cur_c_idx(id:end);
                    
                corner_pts = [corner_pts; cur_c(id,:)];
                introduced_num_corner_points = introduced_num_corner_points+1;
                num_org_cfrags = num_org_cfrags+1;
                i = i-1; % reconsider current contour to check if there are other possible breaking point
            end

        end
        i=i+1;
    end

% end

end