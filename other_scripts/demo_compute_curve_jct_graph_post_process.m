clear all; close all;
addpath (genpath('util/'));

%% set paths
base_src_path = '/media/guoy/Research/Datasets/Drapery_Datasets/Cloth1_COB/';
edge_src_path = [base_src_path 'edges/'];
% jct_src_path = '/media/guoy/Research/Datasets/Junction_curve_structure_Dataset/Data/BSDS500_eval/test_SETO/edges/';
cfrags_dst_path = [base_src_path 'MSEL/'];
mkdir(cfrags_dst_path);
% img_src_path = base_src_path;
img_src_path = '/media/guoy/Research/Datasets/Drapery_Datasets/Cloth1/';
cem_src_path = cfrags_dst_path;
% broken_cfrag_path = 'Data/TO_SEL_CFGD/broken_cfrags/';
% mkdir(broken_cfrag_path);
% final_cem_path = [cfrags_dst_path 'final_curves2/'];
final_cem_path = cfrags_dst_path;
mkdir(final_cem_path);
img_ext = '*.png';

    
%% Extrac Curve Fragments from edges
img_files = dir([img_src_path img_ext]);

for c = 1:length(img_files)
    
    disp([num2str(c) '/'  num2str(length(img_files))]); 
    
    %% Compute cfrags from edges: curvelets, unambiguous cfrags, HPG cfrags
    input = [edge_src_path img_files(c).name(1:end-4) '.edg'];
    output = [cfrags_dst_path img_files(c).name(1:end-4) '.cem'];
    
    cmd = ['!util/dborl_compute_curve_frags ' input ' ' output];
    eval(cmd);
    
    %% Post Process introduce more junctions
    img = imread([img_src_path img_files(c).name]);
    [h,w,~]= size(img);
    
    input_file = [cem_src_path img_files(c).name(1:end-4) '.cem'];
    [CEM, edges, cfrags_idx] = load_contours(input_file);

        
%     test_file2 = [jct_src_path  img_files(c).name(1:end-4) '_pj.bmp'];
%     jct_map = im2double(imread(test_file2));
%     %%%%%%%%%%%%%%%%%%%%%% break curve fragments at conners and junctions
%     tic;
%     [new_cfrags, new_cfrags_idx, junction_pts] = contour_introduce_jct_third_party(CEM{2}, cfrags_idx, jct_map, edges, params);
%     toc;

    tic;
    [new_cfrags, new_cfrags_idx, junction_pts] = contour_introduce_jct(CEM{2}, cfrags_idx, h,w, edges);
    toc;
    
    %%%%%%%%%%%%%%%%%%  save results  %%%%%%%%%%%%%%%%%%%%%%%%%%%
    write_cem_fixed_edge_id([final_cem_path img_files(c).name(1:end-4) '.cem'], new_cfrags, edges, h, w, new_cfrags_idx)

    
    %%%%%%%%%%%%%%%%%  visulization %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     [CEM, edges, new_cfrags_idx] = load_contours([final_cem_path img_files(c).name(1:end-4) '.cem']);
%     new_cfrags = CEM{2};
    % display and save results
    imshow(img, 'border', 'tight'); hold on;
    draw_contours(new_cfrags);
    G_test = construct_fac_graph_from_curve_fragments (new_cfrags_idx, new_cfrags);

    for v =1:length(G_test.var)
        if(length(G_test.var(v).nbrs_fac)>=3)
            c_x = round(edges(G_test.var(v).actual_edge_id, 1)+1);
            c_y = round(edges(G_test.var(v).actual_edge_id, 2)+1);
            plot(c_x, c_y, 'gx');
        end
    end
    hold off;
    H = figure(1);
    set(H, 'PaperPositionMode','auto')
    print(H,'-dpdf','-r0',[final_cem_path img_files(c).name(1:end-4) '.pdf']);
    cmd = ['!pdfcrop ' final_cem_path img_files(c).name(1:end-4) '.pdf ' final_cem_path img_files(c).name(1:end-4) '.pdf'];
    eval(cmd);

        

end
