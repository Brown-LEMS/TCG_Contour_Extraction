clear all; close all;
addpath (genpath('util/'));

%% set paths
edge_src_path = '/home/guoy/Desktop/Chris/silhouettes_hand/edges/';
cfrags_dst_path = '/home/guoy/Desktop/Chris/silhouettes_hand/cfrags/';
mkdir(cfrags_dst_path);
img_src_path = '/home/guoy/Desktop/Chris/silhouettes_hand/img/';
cem_src_path = cfrags_dst_path;
final_cem_path = '/home/guoy/Desktop/Chris/silhouettes_hand/final_curves/';
mkdir(final_cem_path);
beta_src = 'training/';

%% load in beta
prefix = 'TO_SEL_';

input_beta_0 = load([beta_src prefix 'beta_of_cues_for_merging.txt']);
fmean_0 = input_beta_0(1,:);
fstd_0 = input_beta_0(2,:);
beta_0 = input_beta_0(3,:);
beta_0 = beta_0./fstd_0;

input_beta_1 = load([beta_src prefix 'beta_of_geomcon_cue_for_merging.txt']);
fmean_1 = input_beta_1(1,:);
fstd_1 = input_beta_1(2,:);
beta_1 = input_beta_1(3,:);
beta_1 = beta_1./fstd_1;

%% set parameters

nbr_num_edges = 40; % only consider number of edges close to the connecting points for feature computation % also for matching in building ground truth
diag_of_train = 578.275;
% nbr_len_th = 5;
nbr_len_th = 5;

merge_th = 0.4;
merge_th_geom = 0.4;

max_iter = 2;    


% params for computing texton
no = 6;
ss = 1;
ns = 2;
sc = sqrt(2);
el = 2;
k = 64;
fname = sprintf( ...
    'unitex_%.2g_%.2g_%.2g_%.2g_%.2g_%d.mat',no,ss,ns,sc,el,k);
textonData = load(fname); % defines fb,tex,tsim
    
    
%% Extrac Curve Fragments from edges
img_files = dir([img_src_path '*.png']);

for c = 1:length(img_files)
    
    disp([num2str(c) '/'  num2str(length(img_files))]); 
    
    %% compute cfrags from edges
    input = [edge_src_path img_files(c).name(1:end-4) '.edg'];
    output = [cfrags_dst_path img_files(c).name(1:end-4) '.cem'];
    
    cmd = ['!util/dborl_compute_curve_frags ' input ' ' output];
    eval(cmd);
    
    %% merge curve fragments using trained classifier 
    
    disp('merge curve fragments using trained logistic regression classifier...');

    % load in curve fragment candidates
    input_file = [cem_src_path img_files(c).name(1:end-4) '.cem'];
    [CEM, edges, cf_idx] = load_contours(input_file);
    [~, edgemap, thetamap] = load_edg([edge_src_path img_files(c).name(1:end-4) '.edg']);
    
    img = imread([img_src_path img_files(c).name]);
    [h,w,~]= size(img);
    diag = sqrt(h^2+w^2);
    nbr_num_edges = max(round(nbr_num_edges*diag/diag_of_train), 5);
    
    % compute hsv space map
    hsv_img = rgb2hsv(img);
    
    % compute texton map

    tmap = assignTextons(fbRun(textonData.fb,rgb2gray(img)),textonData.tex);
    
    % construct factor graph based on curve fragment candidates
    G = construct_fac_graph_from_curve_fragments (edges, cf_idx, CEM{2,1});
    
    
    % select merge/break nodes just within dim 2 nodes 
    dim1_var_count = 0;
    dim2_var_count = 0;
    dim3_var_count = 0;
    
    
    merged_cem = CEM{2,1};
    merged_cf_idx = cf_idx;
    
    % only solve degree 2 node merging in this iteration
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
            % compute texture continuity and geometric continuity using
            % both c1,c2 as inputs                

            [geom_diff, texture_diff] = degree_2_node_cues(c1,c2, tmap);

            
            % compute the cues diff between portion of conneting curve fragments
            [bg_grad, sat_grad, hue_grad, abs_k, edge_sparsity, wigg, c_len1] = curve_fragment_cues(c1, hsv_img, edgemap);
            c1_cues = [bg_grad; sat_grad; hue_grad; abs_k; edge_sparsity; wigg];
            [bg_grad, sat_grad, hue_grad, abs_k, edge_sparsity, wigg, c_len2] = curve_fragment_cues(c2, hsv_img, edgemap);
            c2_cues = [bg_grad; sat_grad; hue_grad; abs_k; edge_sparsity; wigg];

            cues_diff = abs(c1_cues-c2_cues);

           
            cues_diff = [cues_diff; geom_diff; texture_diff];
            actual_edge = edges(G.var(v).actual_edge_id,:);

            is_merge = 0;
            % for very short curve, just decide merging based on geometry;
            if(c_len1< nbr_len_th || c_len2< nbr_len_th)
                p = 1 / (1 + exp(-([1 geom_diff] - fmean_1)*beta_1'));
                if(p>merge_th_geom)
                    is_merge = 1;
                end
            else
            % probability to merge at this node, 
                p = 1 / (1 + exp(-([1 cues_diff']-fmean_0)*beta_0'));
                if(p>merge_th)
                    is_merge = 1;
                end           
            end
            
            if(is_merge)
                [G, merged_cem, merged_cf_idx] = merge_at_degree_2_node(G, merged_cem, merged_cf_idx, v, actual_edge, nbrs_fac, c1_ids, c2_ids);

            end
        end
        
        % dim_3 node should be dealt with after all dim2 node is solved 
        if(G.var(v).dim ==3)
            dim3_var_count = dim3_var_count +1;
        end
    end    
    
    % Trace closed cirles and merge all the degree 2 nodes on it.
    for v = 1:length(G.var)
        
        % only deal will deg 2 nodes which is not merged yet
        if(G.var(v).dim ==2 && G.var(v).merged == 0)            
            nbrs_fac = G.var(v).nbrs_fac;
            
            % skip circle
            if(nbrs_fac(1) == nbrs_fac(2))
                continue;
            end
            
            % do circle tracing at the node which is not recognized in
            % circle yet           
            if (G.var(v).p ~= 1)
                trace_max_iter = 10;
                trace_iter = 0;
                
                prev_fac_id = nbrs_fac(1);
                prev_var_id = G.var(v).id;
                
                var_id_list = G.var(v).id;
                is_circle = 0;
                while (trace_iter < 10)
                
                    prev_nbr_var_ids = G.fac(prev_fac_id).nbrs_var;
                    curr_var_id = prev_nbr_var_ids(find(prev_nbr_var_ids~=prev_var_id));
                    curr_nbr_fac_ids = G.var(curr_var_id).nbrs_fac;
                    
                    if(var_id_list(1) == curr_var_id)
                        is_circle = 1;
                        break;
                    elseif (G.var(curr_var_id).dim ~=2)
                        break;
                    elseif (G.var(curr_var_id).dim ==2 && var_id_list(1) ~= curr_var_id)
                        var_id_list = [var_id_list curr_var_id];
                        curr_fac_id = curr_nbr_fac_ids(find(curr_nbr_fac_ids~=prev_fac_id));
                        prev_fac_id = curr_fac_id;
                        prev_var_id = curr_var_id;
                        trace_iter = trace_iter +1;
                    end
                    
                end
                
                % when circle is traced, assign the p of each var as 1
                if(is_circle ==1)
%                    keyboard;
                   for v2= 1:length(var_id_list)
                       G.var(var_id_list(v2)).p = 1;
                   end
                end
                
                
            % merge at the node where is assigned p=1 in circle tracing
            elseif(G.var(v).p == 1) 
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

                actual_edge = edges(G.var(v).actual_edge_id,:);
            
                [G, merged_cem, merged_cf_idx] = merge_at_degree_2_node(G, merged_cem, merged_cf_idx, v, actual_edge, nbrs_fac, c1_ids, c2_ids);
            end
        end
    end 
    
    % push all the refined cfrags into another cell
    merged_cem_refine = cell(1,0);
    for i= 1:size(merged_cem,2)
        if(~isempty(merged_cem{i}))
            merged_cem_refine = [merged_cem_refine merged_cem{i}];
        end
    end
    
%     %% post process to break long curve fragment at local mixima
%     
%     disp('post process to insert breaking points...');
%     introduced_num_break_points = 0;
% 
    new_cfrags = merged_cem_refine;
%     for iter = 1:max_iter
%         num_org_cfrags = length(new_cfrags);
%         
% %         nbr_range_th_1 = round(nbr_range_th/iter);
%         for i = 1: num_org_cfrags
%             cur_c = new_cfrags{i};
%             c_len = contour_length_mex(cur_c');
% 
%             % skip closed circle
%             if(cur_c(1,:) == cur_c(end,:))
%                 continue;
%             end
% 
%             % only introducing breaking points for curves which are long enough
%             if(c_len>(30 * diag/ diag_of_train) && size(cur_c,1)>(3*nbr_num_edges/iter))
%                 merge_prob_vec = filter_merge_prob_along_cfrag(cur_c, hsv_img, edgemap, tmap, beta_0, fmean_0, round(nbr_num_edges/iter));
% 
%                 [min_prob, id] = min(merge_prob_vec);
%                 
%                 % should only introduce when it's quite sure
%                 if(min_prob< merge_th-0.1)
%                     introduced_num_break_points = introduced_num_break_points + 1;
%                     new_cfrags = [new_cfrags cur_c(1:id, :)];
%                     new_cfrags{i} = cur_c(id:end, :);
% %                     plot(cur_c(id, 1)+1, cur_c(id, 2)+1, 'yx');
%                 end
% 
%             end
%         end
%     end
%     introduced_num_break_points
    
    % display and save results
    imshow(img, 'border', 'tight'); hold on;
    draw_contours(new_cfrags);
    hold off;
%     print_pdf([final_cem_path img_files(c).name(1:end-4) '.pdf'])
    
    h = figure(1);
    set(h, 'PaperPositionMode','auto')
    print(h,'-dpng','-r0',[final_cem_path img_files(c).name(1:end-4) '.png']);
%     write_cem([final_cem_path img_files(c).name(1:end-4) '.cem'], new_cfrags, h, w)
    write_cem_fixed_edge_id([final_cem_path img_files(c).name(1:end-4) '.cem'], new_cfrags, edges, h, w)

    det_save_cemv([final_cem_path img_files(c).name(1:end-4) '.cemv'], new_cfrags);
    
end