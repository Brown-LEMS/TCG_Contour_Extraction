function [cfrags, cfrags_idx, edges] = break_contours_at_T_junctions(cfrags, cfrags_idx, edges)
% degree of edges are marked, 3 means junctions
% for each contour, look for edges with >=3 degree, and break the contour

for c = 1:length(cfrags_idx)

    cur_c = cfrags{c};
    cur_c_idx = cfrags_idx{c};
    
    edge_degree = edges(cur_c_idx, 5);
    % ignore both end;
    edge_degree(1) = 0;
    edge_degree(end) = 0;
    cur_break_e_ids = find(edge_degree>1);
    
    %%%%%%%%%%%%%  break the cfrag at junctions (one or more)
    %%%%%%%%%%%%% iteratively update at cfrags, finally replace new
    if(~isempty(cur_break_e_ids))
%         cur_break_e_ids = sort(cur_break_e_ids, 'ascend');

        cfrags{c} = cur_c(cur_break_e_ids(end):end, :);
        cfrags = [cfrags cur_c(1:cur_break_e_ids(1), :)];

        cfrags_idx{c} = cur_c_idx(cur_break_e_ids(end):end);
        cfrags_idx = [cfrags_idx cur_c_idx(1:cur_break_e_ids(1))];

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