function  [G, merged_cem, merged_cf_idx] = merge_cfrags_graphical_model_v2(new_cfrags, new_cfrags_idx, edges, hsv_img, tmap, edgemap, params)    
%%%%%% This version merge all the connected curve fragments
%%%%%% for the purpose of evaluating long-range junciton structure
%%%%%%
beta_0 = params.beta_0;
fmean_0 = params.fmean_0;
beta_1 = params.beta_1;
fmean_1 = params.fmean_1;
merge_th = params.merge_th;
merge_th_geom = params.merge_th_geom;
nbr_num_edges = params.nbr_num_edges;
nbr_len_th = params.nbr_len_th;
disp('merge curve fragments under graphical model embeded with trained logistic regressor...');    
    
    
% construct factor graph based on curve fragment candidates
G = construct_fac_graph_from_curve_fragments (edges, new_cfrags_idx, new_cfrags);


% select merge/break nodes just within dim 2 nodes 
dim1_var_count = 0;
dim2_var_count = 0;
dim3_var_count = 0;


merged_cem = new_cfrags;
merged_cf_idx = new_cfrags_idx;

%% decide degree 2 nodes merging in this iteration
for v = 1:length(G.var)
    if(G.var(v).dim ==1)
        dim1_var_count = dim1_var_count +1;
    end

    if(G.var(v).dim ==2) % dim 2 nodes
        dim2_var_count = dim2_var_count +1;

        nbrs_fac = G.var(v).nbrs_fac;

        % do not need to deal with circle anymore
        if(nbrs_fac(1) == nbrs_fac(2))
            continue;
        end

        c1_ids = merged_cf_idx{nbrs_fac(1)};
        c2_ids = merged_cf_idx{nbrs_fac(2)};

        c1 = merged_cem{nbrs_fac(1)};
        c2 = merged_cem{nbrs_fac(2)};



        % only consider portion of the curve fragment within
        % neighbourhood of the connecting node
        if(c1_ids(1) == G.var(v).actual_edge_id)
            c1 = c1(1:min(nbr_num_edges, size(c1,1)), :);
        elseif (c1_ids(end) == G.var(v).actual_edge_id)
            c1 = c1(end-min(nbr_num_edges, size(c1,1))+1:end, :);
        end

        if(c2_ids(1) == G.var(v).actual_edge_id)
            c2 = c2(1:min(nbr_num_edges, size(c2,1)), :);
        elseif (c2_ids(end) == G.var(v).actual_edge_id)
            c2 = c2(end-min(nbr_num_edges, size(c2,1))+1:end, :);
        end            

        % always make directions c1 -> node -> c2
        if(c1_ids(1) == G.var(v).actual_edge_id)
            % reverse the order of edges
            c1 = fliplr(c1');
            c1 = c1';
            c1_ids = fliplr(c1_ids);
        end

        if(c2_ids(end) == G.var(v).actual_edge_id)
            % reverse the order of edges
            c2 = fliplr(c2');
            c2 = c2';    
            c2_ids = fliplr(c2_ids);
        end
%         % compute texture continuity and geometric continuity using
%         % both c1,c2 as inputs                
% 
%         [geom_diff, texture_diff] = degree_2_node_cues(c1,c2, tmap, nbr_num_edges);
% 
% 
%         % compute the cues diff between portion of conneting curve fragments
%         [bg_grad, sat_grad, hue_grad, abs_k, edge_sparsity, wigg, c_len1] = curve_fragment_cues(c1, hsv_img, edgemap);
%         c1_cues = [bg_grad; sat_grad; hue_grad; abs_k; edge_sparsity; wigg];
%         [bg_grad, sat_grad, hue_grad, abs_k, edge_sparsity, wigg, c_len2] = curve_fragment_cues(c2, hsv_img, edgemap);
%         c2_cues = [bg_grad; sat_grad; hue_grad; abs_k; edge_sparsity; wigg];
% 
%         cues_diff = abs(c1_cues-c2_cues);
% 
% 
%         cues_diff = [cues_diff; geom_diff; texture_diff];
        actual_edge = edges(G.var(v).actual_edge_id,:);
% 
%         is_merge = 0;
% %             % for very short curve, just decide merging based on geometry;
%         if(c_len1< min(5, nbr_len_th) || c_len2< min(5, nbr_len_th))
%             p = 1 / (1 + exp(-([1 geom_diff] - fmean_1)*beta_1'));
%             if(p>merge_th_geom)
%                 is_merge = 1;
%             end
%         else
%         % probability to merge at this node, 
%             p = 1 / (1 + exp(-([1 cues_diff']-fmean_0)*beta_0'));
%             if(p > merge_th)
%                 is_merge = 1;
%             end           
%         end
% 
%         G.var(v).p = p;
        G.var(v).p = 1;

        % if p>0.3, merge
%         if(is_merge)
            [G, merged_cem, merged_cf_idx] = merge_at_degree_2_node(G, merged_cem, merged_cf_idx, v, actual_edge, nbrs_fac, c1_ids, c2_ids);
%         end
    end

    %%%%%%%%%%%%% dim_3 node should be dealt with after all dim2 node is solved 
    if(G.var(v).dim ==3)
        dim3_var_count = dim3_var_count +1;
    end
end    

%%  decide dim 3 nodes merging
for v = 1:length(G.var)


    if(G.var(v).dim ==3) % dim 2 nodes
       

        nbrs_fac = G.var(v).nbrs_fac;

        % do not need to deal with circle anymore
        if(nbrs_fac(1) == nbrs_fac(2) || nbrs_fac(1) == nbrs_fac(3) || nbrs_fac(2) == nbrs_fac(3))
            continue;
        end

        c1_ids = merged_cf_idx{nbrs_fac(1)};
        c2_ids = merged_cf_idx{nbrs_fac(2)};
        c3_ids = merged_cf_idx{nbrs_fac(3)};

        c1 = merged_cem{nbrs_fac(1)};
        c2 = merged_cem{nbrs_fac(2)};
        c3 = merged_cem{nbrs_fac(3)};



        % only consider portion of the curve fragment within
        % neighbourhood of the connecting node
        % always make directions c1 -> node <- c2
        %                                ^
        %                                |
        %                               c3
        %
        if(c1_ids(1) == G.var(v).actual_edge_id)
            c1 = c1(1:min(nbr_num_edges, size(c1,1)), :);
            % reverse the order of edges
            c1 = fliplr(c1');
            c1 = c1';
            c1_ids = fliplr(c1_ids);
        elseif (c1_ids(end) == G.var(v).actual_edge_id)
            c1 = c1(end-min(nbr_num_edges, size(c1,1))+1:end, :);
        end

        if(c2_ids(1) == G.var(v).actual_edge_id)
            c2 = c2(1:min(nbr_num_edges, size(c2,1)), :);
            % reverse the order of edges
            c2 = fliplr(c2');
            c2 = c2';    
            c2_ids = fliplr(c2_ids);
        elseif (c2_ids(end) == G.var(v).actual_edge_id)
            c2 = c2(end-min(nbr_num_edges, size(c2,1))+1:end, :);
        end    
        
        if(c3_ids(1) == G.var(v).actual_edge_id)
            c3 = c3(1:min(nbr_num_edges, size(c3,1)), :);
            % reverse the order of edges
            c3 = fliplr(c3');
            c3 = c3';    
            c3_ids = fliplr(c3_ids);
        elseif (c3_ids(end) == G.var(v).actual_edge_id)
            c3 = c3(end-min(nbr_num_edges, size(c3,1))+1:end, :);
        end   


        %%%%%%%%%%%%%%%  compute the pairwise diff and merge probs at
        %%%%%%%%%%%%%%%  degree 3 node

        [cues_diff, probs] = degree_3_node_cues(c1,c2,c3, tmap, hsv_img, edgemap, fmean_0, beta_0, fmean_1, beta_1, nbr_len_th, nbr_num_edges);

        actual_edge = edges(G.var(v).actual_edge_id,:);

        [max_p, max_p_id] = max(probs);

        % probability to merge at this node, 
       
        if(max_p > merge_th)
            is_merge = 1;
        else
            is_merge = 0;
        end           
       
        G.var(v).p = probs;

        if(is_merge)
            [G, merged_cem, merged_cf_idx] = merge_at_degree_3_node(G, merged_cem, merged_cf_idx, v, actual_edge, c1_ids, c2_ids, c3_ids);
        end
    end

    
end

% %% Trace closed cirles and merge all the degree 2 nodes on it.
% for v = 1:length(G.var)
% 
%     % only deal will deg 2 nodes which is not merged yet
%     if(G.var(v).dim ==2 && G.var(v).merged == 0)            
%         nbrs_fac = G.var(v).nbrs_fac;
% 
%         % skip circle
%         if(nbrs_fac(1) == nbrs_fac(2))
%             continue;
%         end
% 
%         % do circle tracing at the node which is not recognized in
%         % circle yet           
%         if (G.var(v).p ~= 1)
%             trace_max_iter = 10;
%             trace_iter = 0;
% 
%             prev_fac_id = nbrs_fac(1);
%             prev_var_id = G.var(v).id;
% 
%             var_id_list = G.var(v).id;
%             is_circle = 0;
%             while (trace_iter < 10)
% 
%                 prev_nbr_var_ids = G.fac(prev_fac_id).nbrs_var;
%                 curr_var_id = prev_nbr_var_ids(find(prev_nbr_var_ids~=prev_var_id));
%                 curr_nbr_fac_ids = G.var(curr_var_id).nbrs_fac;
% 
%                 if(var_id_list(1) == curr_var_id)
%                     is_circle = 1;
%                     break;
%                 elseif (G.var(curr_var_id).dim ~=2)
%                     break;
%                 elseif (G.var(curr_var_id).dim ==2 && var_id_list(1) ~= curr_var_id)
%                     var_id_list = [var_id_list curr_var_id];
%                     curr_fac_id = curr_nbr_fac_ids(find(curr_nbr_fac_ids~=prev_fac_id));
%                     prev_fac_id = curr_fac_id;
%                     prev_var_id = curr_var_id;
%                     trace_iter = trace_iter +1;
%                 end
% 
%             end
% 
%             % when circle is traced, assign the p of each var as 1
%             if(is_circle ==1)
% %                    keyboard;
%                for v2= 1:length(var_id_list)
%                    G.var(var_id_list(v2)).p = 1;
%                end
%             end
% 
% 
%         % merge at the node where is assigned p=1 in circle tracing
%         elseif(G.var(v).p == 1) 
%             c1_ids = merged_cf_idx{nbrs_fac(1)};
%             c2_ids = merged_cf_idx{nbrs_fac(2)};
% 
%             c1 = merged_cem{nbrs_fac(1)};
%             c2 = merged_cem{nbrs_fac(2)};
% 
%             % only consider portion of the curve fragment within
%             % neighbourhood of the connecting node
%             if(c1_ids(1) == G.var(v).actual_edge_id)
%                 c1 = c1(1:min(nbr_num_edges, size(c1,1)), :);
%             elseif (c1_ids(end) == G.var(v).actual_edge_id)
%                 c1 = c1(end-min(nbr_num_edges, size(c1,1))+1:end, :);
%             end
% 
%             if(c2_ids(1) == G.var(v).actual_edge_id)
%                 c2 = c2(1:min(nbr_num_edges, size(c2,1)), :);
%             elseif (c2_ids(end) == G.var(v).actual_edge_id)
%                 c2 = c2(end-min(nbr_num_edges, size(c2,1))+1:end, :);
%             end            
% 
%             % always make directions c1 -> node -> c2
%             if(c1_ids(1) == G.var(v).actual_edge_id)
%                 % reverse the order of edges
%                 c1 = fliplr(c1');
%                 c1 = c1';
%                 c1_ids = fliplr(c1_ids);
%             end
% 
%             if(c2_ids(end) == G.var(v).actual_edge_id)
%                 % reverse the order of edges
%                 c2 = fliplr(c2');
%                 c2 = c2';    
%                 c2_ids = fliplr(c2_ids);
%             end
% 
%             actual_edge = edges(G.var(v).actual_edge_id,:);
% 
%             [G, merged_cem, merged_cf_idx] = merge_at_degree_2_node(G, merged_cem, merged_cf_idx, v, actual_edge, nbrs_fac, c1_ids, c2_ids);
%         end
%     end
% end


% push all the refined cfrags into another cell
% remove those length 2 curve fragments
merged_cem_refined = cell(1,0);
merged_cf_idx_refined = cell(1,0);
for i= 1:size(merged_cem,2)
    if(~isempty(merged_cem{i}))
        merged_cem_refined = [merged_cem_refined merged_cem{i}];
        merged_cf_idx_refined = [merged_cf_idx_refined merged_cf_idx{i}];
    end
end

merged_cem = merged_cem_refined;
merged_cf_idx = merged_cf_idx_refined;