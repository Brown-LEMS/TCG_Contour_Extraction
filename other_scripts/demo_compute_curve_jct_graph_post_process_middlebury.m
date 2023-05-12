clear all; close all;
addpath (genpath('util'))
img_src_path = '/media/guoy/Research/Datasets/MiddEval/trainingQ/';
edge_src_path = '/media/guoy/Research/Datasets/MiddEval/trainingQ_SE/';
frags_dst_path = '/media/guoy/Research/Datasets/MiddEval/trainingQ_SE/';
final_cem_path = frags_dst_path;
mkdir(frags_dst_path);
beta_src = 'training/';

%% load in beta
prefix = 'SE_SEL_';

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

input_beta_2 = load([beta_src prefix 'beta_of_cues_for_seletion.txt']);
fmean_2 = input_beta_2(1,:);
fstd_2 = input_beta_2(2,:);
beta_2 = input_beta_2(3,:);
beta_2 = beta_2./fstd_2;

params.beta_0 = beta_0;
params.fmean_0 = fmean_0;
params.beta_1 = beta_1;
params.fmean_1 = fmean_1;
params.beta_2 = beta_2;
params.fmean_2 = fmean_2;

%% set parameters

nbr_num_edges = 20; % only consider number of edges close to the connecting points for feature computation % also for matching in building ground truth
params.diag_of_train = 578.275;
params.nbr_len_th = 5; % short curve under this length will be grouped due to geometry.
params.merge_th = 0.4;
params.merge_th_geom = 0.4;
params.max_iter = 2;    

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

%% processing data
image_names{1} = 'Adirondack';
image_names{2} = 'ArtL';
image_names{3} = 'Jadeplant';
image_names{4} = 'Motorcycle';
image_names{5} = 'MotorcycleE';
image_names{6} = 'Piano';
image_names{7} = 'PianoL';
image_names{8} = 'Pipes';
image_names{9} = 'Playroom';
image_names{10} = 'Playtable';
image_names{11} = 'PlaytableP';
image_names{12} = 'Recycle';
image_names{13} = 'Shelves';
image_names{14} = 'Teddy';
image_names{15} = 'Vintage';

% img_files = dir([img_src_path img_ext]);

for c = 1:length(image_names)
    
    disp([num2str(c) '/'  num2str(length(image_names))]); 
    
    %% im0:  Compute cfrags from edges: curvelets, unambiguous cfrags, HPG cfrags
    input = [edge_src_path image_names{c} '/im0.edg'];
    output = [frags_dst_path image_names{c} '/im0.cem'];
    
    cmd = ['!util/dborl_compute_curve_frags ' input ' ' output];
%     cmd = ['!util/dborl_compute_curve_frags_pos_uncertain05 ' input ' ' output];
    eval(cmd);
    
    
    %%%%%%%%%%%%%%%%%% Post Process introduce more junctions
    img = imread([img_src_path image_names{c} '/im0.png']);
    [h,w,~]= size(img);
    diag = sqrt(h^2+w^2);
    params.nbr_num_edges = max(round(nbr_num_edges*diag/params.diag_of_train), 5);
    params.diag_ratio = diag/ params.diag_of_train;
    [~, edgemap, ~] = load_edg([edge_src_path image_names{c} '/im0.edg']);
    
    
    [CEM, edges, cfrags_idx] = load_contours(output);    
    tic;
    [new_cfrags, new_cfrags_idx, T_junctions] = contour_introduce_jct(CEM{2}, cfrags_idx, h,w, edges);
    toc;
    %%%%%%%%%%%%%%%%%% Merge curve fragments under graphical model
    tic
    [G_test, new_cfrags, new_cfrags_idx, T_junctions] = merge_cfrags_graphical_model_geom(new_cfrags, new_cfrags_idx, edges, params);
    toc 
    %%%%%%%%%%%%%%%%%%%%%% break curve fragments at conners
    tic
    [new_cfrags, new_cfrags_idx, corner_pts] = contour_breaker_at_conner(new_cfrags, new_cfrags_idx, params);
    toc
    %%%%%%%%%%%%%%%%%%%%%% classify junction types and merge at T-junction
    tic
    [G_test, new_cfrags, new_cfrags_idx, T_junctions] = classify_junction_type_wrt_graph(new_cfrags, new_cfrags_idx, edges, params);
    toc    
       
    %%%%%%%%%%%%%%%%%  visulization %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    G_test = construct_fac_graph_from_curve_fragments (new_cfrags_idx, new_cfrags);
    H = figure(1);
    imshow(img, 'border', 'tight'); hold on;
    draw_contours(new_cfrags);
    if(~isempty(T_junctions))
        plot(T_junctions(:,1)+1, T_junctions(:,2)+1, 'rx');
    end
%     if(~isempty(corner_pts))
%         plot(corner_pts(:,1)+1, corner_pts(:,2)+1, 'yo');
%     end
    for v =1:length(G_test.var)
        if(length(G_test.var(v).nbrs_fac)>=3)
            c_x = edges(G_test.var(v).actual_edge_id, 1)+1;
            c_y = edges(G_test.var(v).actual_edge_id, 2)+1;
            plot(c_x, c_y, 'gx');
        end
    end
    hold off;
    set(H, 'PaperPositionMode','auto')
    print(H,'-dpdf','-r0',[frags_dst_path image_names{c} 'im0.pdf']);
    cmd = ['!pdfcrop ' frags_dst_path image_names{c} 'im0.pdf ' frags_dst_path image_names{c} 'im0.pdf'];
    eval(cmd);
    
    %%%%%%%%%%%%%%%%%%  save results  %%%%%%%%%%%%%%%%%%%%%%%%%%%
    write_cem_fixed_edge_id(output, new_cfrags, edges, h, w, new_cfrags_idx)
    
    %% im1: Compute cfrags from edges: curvelets, unambiguous cfrags, HPG cfrags
    input = [edge_src_path image_names{c} '/im1.edg'];
    output = [frags_dst_path image_names{c} '/im1.cem'];
    
    cmd = ['!util/dborl_compute_curve_frags ' input ' ' output];
    eval(cmd);
    
    %%%%%%%%%%%%%%%%%% Post Process introduce more junctions
    img = imread([img_src_path image_names{c} '/im1.png']);
    [h,w,~]= size(img);
    
    [~, edgemap, ~] = load_edg([edge_src_path image_names{c} '/im1.edg']);
    
    
    [CEM, edges, cfrags_idx] = load_contours(output);    
    tic;
    [new_cfrags, new_cfrags_idx, T_junctions] = contour_introduce_jct(CEM{2}, cfrags_idx, h,w, edges);
    toc;
    %%%%%%%%%%%%%%%%%% Merge curve fragments under graphical model
    tic
    [G_test, new_cfrags, new_cfrags_idx, T_junctions] = merge_cfrags_graphical_model_geom(new_cfrags, new_cfrags_idx, edges, params);
    toc        
    %%%%%%%%%%%%%%%%%%%%%% break curve fragments at conners
    tic
    [new_cfrags, new_cfrags_idx, corner_pts] = contour_breaker_at_conner(new_cfrags, new_cfrags_idx, params);
    toc    
    %%%%%%%%%%%%%%%%%%%%%% classify junction types and merge at T-junction
    tic
    [G_test, new_cfrags, new_cfrags_idx, T_junctions] = classify_junction_type_wrt_graph(new_cfrags, new_cfrags_idx, edges, params);
    toc
    
    %%%%%%%%%%%%%%%%%  visulization %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    G_test = construct_fac_graph_from_curve_fragments (new_cfrags_idx, new_cfrags);
    H = figure(2);
    imshow(img, 'border', 'tight'); hold on;
    draw_contours(new_cfrags);
    if(~isempty(T_junctions))
        plot(T_junctions(:,1)+1, T_junctions(:,2)+1, 'rx');
    end
%     if(~isempty(corner_pts))
%         plot(corner_pts(:,1)+1, corner_pts(:,2)+1, 'yo');
%     end
    for v =1:length(G_test.var)
        if(length(G_test.var(v).nbrs_fac)>=3)
            c_x = edges(G_test.var(v).actual_edge_id, 1)+1;
            c_y = edges(G_test.var(v).actual_edge_id, 2)+1;
            plot(c_x, c_y, 'gx');
        end
    end
    hold off;
    set(H, 'PaperPositionMode','auto')
    print(H,'-dpdf','-r0',[frags_dst_path image_names{c} 'im1.pdf']);
    cmd = ['!pdfcrop ' frags_dst_path image_names{c} 'im1.pdf ' frags_dst_path image_names{c} 'im1.pdf'];
    eval(cmd);
    
    %%%%%%%%%%%%%%%%%%  save results  %%%%%%%%%%%%%%%%%%%%%%%%%%%
    write_cem_fixed_edge_id(output, new_cfrags, edges, h, w, new_cfrags_idx)
    
end