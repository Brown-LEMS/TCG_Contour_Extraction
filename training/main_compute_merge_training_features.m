clear all; close all;
addpath (genpath('../util/'));

%% construct Graph from curve fragments, label degree 2 nodes 0/1 and compute cues difference between connecting curve fragments

cost_thresh = 2; % thresh For pruning out unrelated curve fragments and grouping process 
nbr_num_edges_org = 20; % only consider number of edges close to the connecting points for feature computation % also for matching in building ground truth
diag_of_train = 578.275;

img_src_path = '../Data/CFGD_img/';
edg_src_path = '../Data/gPb_SEL_CFGD/edges/';
cem_src_path = '../Data/gPb_SEL_CFGD/broken_cfrags/';
% vis_dst_path = 'Tests/TO_SEL_cfrags_CFGD2_vis/';
% mkdir(vis_dst_path);
% hist_dst_path = 'hists/TO_SEL_hists_cues_for_merging/';
% mkdir(hist_dst_path);
gt_src_path = '../Data/CFGD/';

prefix = 'gPb_SEL_';


cem_files = dir([cem_src_path '*.cem']);

Y = [];
Features = [];
for c = 1:length(cem_files)
    
    disp([num2str(c) '/'  num2str(length(cem_files))]); 
    % load in curve fragment candidates
    input_file = [cem_src_path cem_files(c).name];
    [CEM, edges, cf_idx] = load_contours(input_file);
    [~, edgemap, thetamap] = load_edg([edg_src_path cem_files(c).name(1:end-3) 'edg']);
    
    img = imread([img_src_path cem_files(c).name(1:end-4) '.jpg']);
    
    % compute hsv space map
    hsv_img = rgb2hsv(img);
    
    % compute texton map
    no = 6;
    ss = 1;
    ns = 2;
    sc = sqrt(2);
    el = 2;
    k = 64;
    fname = sprintf( ...
        'unitex_%.2g_%.2g_%.2g_%.2g_%.2g_%d.mat',no,ss,ns,sc,el,k);
    textonData = load(fname); % defines fb,tex,tsim
    tmap = assignTextons(fbRun(textonData.fb,rgb2gray(img)),textonData.tex);


    % visulize contour fragments
    [h,w,~] = size(img);
    diag = sqrt(h^2+w^2);
    nbr_num_edges = max(round(nbr_num_edges_org*diag/diag_of_train), 10);    
    

    
    % construct factor graph based on curve fragment candidates
    G = construct_fac_graph_from_curve_fragments (edges, cf_idx, CEM{2,1});
    imshow(zeros(h,w), 'border', 'tight');hold on;
    draw_contours(CEM{2,1});    
    for gt_num = 1:3

        % load in ground truth candidates
        gt_file = [gt_src_path cem_files(c).name(1:end-4) '_s' num2str(gt_num) '.cem'];
        CEM_gt = load_contours(gt_file);


        % select merge/break nodes just within dim 2 nodes 
        dim1_var_count = 0;
        dim2_var_count = 0;
        dim3_var_count = 0;

        % selected sample node id vector, actual edge id vector, label vector
        sample_id_vector = [];
        sample_edge_id_vector = [];
        sample_label_vector = [];
        cues_diff_vector = [];


        for v = 1:length(G.var)
%             if(G.var(v).gt_label ~= -1)
%                 continue;
%             end
            
            if(G.var(v).dim ==1)
                dim1_var_count = dim1_var_count +1;
            end

            if(G.var(v).dim ==2) % investigate only on these dim 2 nodes
                dim2_var_count = dim2_var_count +1;

                nbrs_fac = G.var(v).nbrs_fac;

                c1_ids = cf_idx{nbrs_fac(1)};
                c2_ids = cf_idx{nbrs_fac(2)};

                c1 = CEM{2,1}{nbrs_fac(1)};
                c2 = CEM{2,1}{nbrs_fac(2)};

                % do not consider the extremely short curves which do not have
                % curve property
                if(size(c1,1)<=nbr_num_edges || size(c2,1)<=nbr_num_edges)
                    continue;
                end

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

                c1_gt_cost_vector = ones(length(CEM_gt{2,1}),1)*1000;
                c2_gt_cost_vector = ones(length(CEM_gt{2,1}),1)*1000;

                for i=1:length(CEM_gt{2,1})
                    cur_c = CEM_gt{2,1}{i};
                    c1_gt_cost_vector(i) = compute_contours_cost_mex(c1', cur_c', cost_thresh);
                end

                for i=1:length(CEM_gt{2,1})
                    cur_c = CEM_gt{2,1}{i};
                    c2_gt_cost_vector(i) = compute_contours_cost_mex(c2', cur_c', cost_thresh);
                end

                % compute a dist between c1 c2 to judge if they are connected
                % by a very sharp angle, might be matech the same ground truth
                % cf
% 
%                 c1_c2_cost = min(compute_contours_cost_mex(c2', c1', cost_thresh),compute_contours_cost_mex(c1', c2', cost_thresh));

                match_gt_id1 = find(c1_gt_cost_vector<1000);            
                match_gt_id2 = find(c2_gt_cost_vector<1000);

                % CASE: Both side curve fragments get matched to GT. 
                if(length(match_gt_id1)==1&&length(match_gt_id2)==1)


                    % compute the cues diff between portion of conneting curve fragments
                    [bg_grad, sat_grad, hue_grad, abs_k, edge_sparsity, wigg] = curve_fragment_cues(c1, hsv_img, edgemap);
                    c1_cues = [bg_grad; sat_grad; hue_grad; abs_k; edge_sparsity; wigg];
                    [bg_grad, sat_grad, hue_grad, abs_k, edge_sparsity, wigg] = curve_fragment_cues(c2, hsv_img, edgemap);
                    c2_cues = [bg_grad; sat_grad; hue_grad; abs_k; edge_sparsity; wigg];

                    cues_diff = abs(c1_cues-c2_cues);


                    % compute texture continuity and geometric continuity using
                    % both c1,c2 as inputs                
                    % always make directions c1 -> node -> c2
                    if(c1_ids(1) == G.var(v).actual_edge_id)
                        % reverse the order of edges
                        c1 = fliplr(c1');
                        c1 = c1';
                    end

                    if(c2_ids(end) == G.var(v).actual_edge_id)
                        % reverse the order of edges
                        c2 = fliplr(c2');
                        c2 = c2';                    
                    end
                    [geom_diff, texture_diff] = degree_2_node_cues(c1,c2, tmap, nbr_num_edges);
                    cues_diff = [cues_diff; geom_diff; texture_diff];

                    actual_edge = edges(G.var(v).actual_edge_id,1:2);

                    % SUBCASE: if both side cf match to the same GT, + sample
                    if(match_gt_id1==match_gt_id2)
                        if(G.var(v).gt_label == -1)
                            G.var(v).gt_label = 1;
                            sample_label_vector = [sample_label_vector 1];
                            sample_id_vector = [sample_id_vector v];
                            sample_edge_id_vector = [sample_edge_id_vector G.var(v).actual_edge_id];
                            cues_diff_vector = [cues_diff_vector cues_diff];
                            plot(actual_edge(1)+1, actual_edge(2)+1, 'go', 'MarkerSize',10);
                        elseif(G.var(v).gt_label ~= -1 &&  G.var(v).gt_label ~= 1)
                            node_loc = find(sample_edge_id_vector==G.var(v).actual_edge_id);
                            sample_id_vector(node_loc) = [];
                            sample_edge_id_vector(node_loc) = [];
                            sample_label_vector(node_loc) = [];
                            cues_diff_vector(:, node_loc) = [];
                        end
                    else
                        gt_c1 = CEM_gt{2,1}{match_gt_id1};
                        gt1_start_edge = gt_c1(1,1:2);
                        gt1_end_edge = gt_c1(end,1:2);
                        gt_c2 = CEM_gt{2,1}{match_gt_id2};
                        gt2_start_edge = gt_c2(1,1:2);
                        gt2_end_edge = gt_c2(end,1:2);
                        % SUBCASE: if side cfs match to the different GT cfs, and GT cfs connected, and the connected node not far + sample     
                        min_gt_dist = min([pdist2(gt1_start_edge, gt2_start_edge, 'euclidean'), pdist2(gt1_start_edge, gt2_end_edge, 'euclidean') ...
                            pdist2(gt1_end_edge, gt2_start_edge, 'euclidean'), pdist2(gt1_end_edge, gt2_end_edge, 'euclidean')]);
                        min_gt_cp_node = min([pdist2(gt1_start_edge, actual_edge, 'euclidean'), pdist2(gt1_end_edge, actual_edge, 'euclidean')...
                            pdist2(gt2_start_edge, actual_edge, 'euclidean'), pdist2(gt2_end_edge, actual_edge, 'euclidean')]);
                        
                        if(min_gt_dist < 5 && min_gt_cp_node > 3)
                            if(G.var(v).gt_label == -1)
                                G.var(v).gt_label = 1;
                                sample_label_vector = [sample_label_vector 1];
                                sample_id_vector = [sample_id_vector v];
                                sample_edge_id_vector = [sample_edge_id_vector G.var(v).actual_edge_id];
                                cues_diff_vector = [cues_diff_vector cues_diff];
                                plot(actual_edge(1)+1, actual_edge(2)+1, 'gs', 'MarkerSize',10);   
                            elseif(G.var(v).gt_label ~= -1 &&  G.var(v).gt_label ~= 1)
                                node_loc = find(sample_edge_id_vector==G.var(v).actual_edge_id);
                                sample_id_vector(node_loc) = [];
                                sample_edge_id_vector(node_loc) = [];
                                sample_label_vector(node_loc) = [];
                                cues_diff_vector(:, node_loc) = [];
                            end
                        elseif(min_gt_dist < 5 && min_gt_cp_node < 3)
                            if(G.var(v).gt_label == -1)
                                G.var(v).gt_label = 0;
                                sample_label_vector = [sample_label_vector 0];
                                sample_id_vector = [sample_id_vector v];
                                sample_edge_id_vector = [sample_edge_id_vector G.var(v).actual_edge_id];
                                cues_diff_vector = [cues_diff_vector cues_diff];
                                plot(actual_edge(1)+1, actual_edge(2)+1, 'rs', 'MarkerSize',10);
                            elseif(G.var(v).gt_label ~= -1 &&  G.var(v).gt_label ~= 0)
                                node_loc = find(sample_edge_id_vector==G.var(v).actual_edge_id);
                                sample_id_vector(node_loc) = [];
                                sample_edge_id_vector(node_loc) = [];
                                sample_label_vector(node_loc) = [];
                                cues_diff_vector(:, node_loc) = [];
                            end
                        else
                            if(G.var(v).gt_label == -1)
                                G.var(v).gt_label = 0;
                                sample_label_vector = [sample_label_vector 0];
                                sample_id_vector = [sample_id_vector v];
                                sample_edge_id_vector = [sample_edge_id_vector G.var(v).actual_edge_id];
                                cues_diff_vector = [cues_diff_vector cues_diff];
                                plot(actual_edge(1)+1, actual_edge(2)+1, 'ro', 'MarkerSize',10);
                            elseif(G.var(v).gt_label ~= -1 &&  G.var(v).gt_label ~= 0)
                                node_loc = find(sample_edge_id_vector==G.var(v).actual_edge_id);
                                sample_id_vector(node_loc) = [];
                                sample_edge_id_vector(node_loc) = [];
                                sample_label_vector(node_loc) = [];
                                cues_diff_vector(:, node_loc) = [];
                            end
                        end

                    end

    %             % CASE: one match to GT one match to nothing, add to neg sample
    %             elseif((length(match_gt_id1)==1&&length(match_gt_id2)==0) || (length(match_gt_id1)==0&&length(match_gt_id2)==1))
    %                 sample_id_vector = [sample_id_vector v];
    %                 sample_edge_id_vector = [sample_edge_id_vector G.var(v).actual_edge_id];
    % 
    %                 % compute the cues diff between portion of conneting curve fragments
    %                 [bg_grad, sat_grad, hue_grad, abs_k, edge_sparsity, wigg] = curve_fragment_cues(c1, hsv_img, edgemap);
    %                 c1_cues = [bg_grad; sat_grad; hue_grad; abs_k; edge_sparsity; wigg];
    %                 [bg_grad, sat_grad, hue_grad, abs_k, edge_sparsity, wigg] = curve_fragment_cues(c2, hsv_img, edgemap);
    %                 c2_cues = [bg_grad; sat_grad; hue_grad; abs_k; edge_sparsity; wigg];
    %                 
    %                 cues_diff = abs(c1_cues-c2_cues);
    %                 
    %                 
    %                 % compute texture continuity and geometric continuity using
    %                 % both c1,c2 as inputs                
    %                 % always make directions c1 -> node -> c2
    %                 if(c1_ids(1) == G.var(v).actual_edge_id)
    %                     % reverse the order of edges
    %                     c1 = fliplr(c1');
    %                     c1 = c1';
    %                 end
    %                 
    %                 if(c2_ids(end) == G.var(v).actual_edge_id)
    %                     % reverse the order of edges
    %                     c2 = fliplr(c2');
    %                     c2 = c2';                    
    %                 end
    %                 [geom_diff, texture_diff] = degree_2_node_cues(c1,c2, tmap);
    %                 cues_diff = [cues_diff; geom_diff; texture_diff];
    %                 cues_diff_vector = [cues_diff_vector cues_diff];
    %                 
    %                 actual_edge = edges(G.var(v).actual_edge_id,:);   
    %                 G.var(v).gt_label = 0;
    %                 sample_label_vector = [sample_label_vector 0];
    %                 plot(actual_edge(1)+1, actual_edge(2)+1, 'rx', 'MarkerSize',10);
                end

    %             disp(G.var(v).gt_label);
    %             keyboard;

            end
            if(G.var(v).dim ==3)
                dim3_var_count = dim3_var_count +1;
                actual_edge = edges(G.var(v).actual_edge_id,1:2);
%                 plot(actual_edge(1)+1, actual_edge(2)+1, 'co', 'MarkerSize',10);
            end
        end    

    %     print_pdf([vis_dst_path cem_files(c).name(1:end-4) '_select_samples.pdf'])

        onidx = find(sample_label_vector==1);
        offidx = find(sample_label_vector==0); % offidx -> break point is small in number,

        cnt = numel(onidx);
        idx = randperm(cnt);
        onidx = onidx(idx);

        ind = [ offidx onidx ];

        if(numel(onidx)>numel(offidx))
            cnt2 = 2*numel(offidx);
        else
            cnt2 = 2*numel(onidx);
        end
        idx = randperm(cnt2);


        Y = [Y sample_label_vector(ind(idx))];
        Features = [Features cues_diff_vector(:,(ind(idx)))];
%         keyboard;
    end
    hold off;

end

save([prefix 'Features.mat'], 'Features');
save([prefix 'Y.mat'], 'Y');