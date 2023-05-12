clear all; close all;
addpath (genpath('../util/'));

%% construct Graph from curve fragments, label degree 2 nodes 0/1 and compute cues difference between connecting curve fragments

maxDist = 0.01
cost_thresh = 2; % thresh For pruning out unrelated curve fragments and grouping process 
nbr_num_edges = 20; % only consider number of edges close to the connecting points for feature computation % also for matching in building ground truth
diag_of_train = 578.275;

img_src_path = '../Data/VOC2007_img/';
edg_src_path = '../Data/gPb_SEL_VOC2007/edges/';
cem_src_path = '../Data/gPb_SEL_VOC2007/cfrags/';
% vis_dst_path = 'Tests/TO_SEL_cfrags_CFGD2_vis/';
% mkdir(vis_dst_path);
% hist_dst_path = 'hists/TO_SEL_hists_cues_for_merging/';
% mkdir(hist_dst_path);
gt_src_path = '../Data/VOC2007_GT/';

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
    nbr_num_edges = max(round(nbr_num_edges*diag/diag_of_train), 5);    
    

    
    % construct factor graph based on curve fragment candidates
    G = construct_fac_graph_from_curve_fragments (edges, cf_idx, CEM{2,1});
    imshow(zeros(h,w), 'border', 'tight');hold on;
    draw_contours(CEM{2,1});    
    
    load([gt_src_path cem_files(c).name(1:end-4) '.mat' ]);

    gt_edgemap_objects = zeros(size(gt_edgemap{1}));
    disp('construct gt edgemaps of each object silhouette')
    for gt_num = 1:length(gt_edgemap)
        gt_edgemap_objects = gt_edgemap_objects + gt_num*gt_edgemap{gt_num};
    end
    
%     imagesc(gt_edgemap_objects);
%     keyboard;
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

            % find the gt_edges c1, and c2 matches

            cur_Y = b_merge_node(c1, c2, gt_edgemap_objects, maxDist);


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
            [geom_diff, texture_diff] = degree_2_node_cues(c1,c2, tmap);
            cues_diff = [cues_diff; geom_diff; texture_diff];

            actual_edge = edges(G.var(v).actual_edge_id,1:2);

                
                
            if(cur_Y==1)
                G.var(v).gt_label = 1;
                sample_label_vector = [sample_label_vector 1];
                sample_id_vector = [sample_id_vector v];
                sample_edge_id_vector = [sample_edge_id_vector G.var(v).actual_edge_id];
                cues_diff_vector = [cues_diff_vector cues_diff];
                plot(actual_edge(1)+1, actual_edge(2)+1, 'gs', 'MarkerSize',10);
            elseif(cur_Y==0)
                G.var(v).gt_label = 0;
                sample_label_vector = [sample_label_vector 0];
                sample_id_vector = [sample_id_vector v];
                sample_edge_id_vector = [sample_edge_id_vector G.var(v).actual_edge_id];
                cues_diff_vector = [cues_diff_vector cues_diff];
                plot(actual_edge(1)+1, actual_edge(2)+1, 'rs', 'MarkerSize',10);
            end
            
        end
 
        if(G.var(v).dim ==3)
            dim3_var_count = dim3_var_count +1;
            actual_edge = edges(G.var(v).actual_edge_id,1:2);
%                 plot(actual_edge(1)+1, actual_edge(2)+1, 'co', 'MarkerSize',10);
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
    end
    hold off;

%     keyboard;
end

save([prefix '_merge_Features.mat'], 'Features', 'Y');