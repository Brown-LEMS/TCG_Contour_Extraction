function [new_cfrags, new_cfrags_idx, break_pts] = contour_breaker_semantic(cfrags, cfrags_idx, hsv_img, tmap, edgemap, params)

beta_0 = params.beta_0;
fmean_0 = params.fmean_0;
max_iter = params.max_iter;
merge_th = params.merge_th;
nbr_num_edges = params.nbr_num_edges;
diag_ratio = params.diag_ratio; 

disp('semantic breaking contours...');
introduced_num_break_points = 0;
break_pts = [];

new_cfrags = cfrags;
new_cfrags_idx = cfrags_idx;

for iter = 1:max_iter
    num_org_cfrags = length(new_cfrags);

%         nbr_range_th_1 = round(nbr_range_th/iter);
    for i = 1: num_org_cfrags
        cur_c = new_cfrags{i};
        cur_c_idx = new_cfrags_idx{i};

        c_len = contour_length_mex(cur_c');

%         % skip closed circle
%         if(cur_c(1,:) == cur_c(end,:))
%             continue;
%         end

        % only introducing breaking points for curves which are long enough
        if(c_len>(30 * diag_ratio) && size(cur_c,1)>(3*nbr_num_edges/iter))

            merge_prob_vec = filter_merge_prob_along_cfrag(cur_c, hsv_img, edgemap, tmap, beta_0, fmean_0, round(nbr_num_edges/iter));

            [min_prob, id] = min(merge_prob_vec);
            id = floor(median(id)); % in case there are sequntial id with the same min_prob, due to down sampling

            % should only introduce when it's quite sure
            if(min_prob< merge_th)
                introduced_num_break_points = introduced_num_break_points + 1;
                new_cfrags = [new_cfrags cur_c(1:id, :)];
                new_cfrags{i} = cur_c(id:end, :);

                new_cfrags_idx = [new_cfrags_idx cur_c_idx(1:id)];                   
                new_cfrags_idx{i} = cur_c_idx(id:end);

                break_pts = [break_pts; cur_c(id,:)];
%                     plot(cur_c(id, 1)+1, cur_c(id, 2)+1, 'yx');
            end

        end
    end
end
% introduced_num_break_points


end