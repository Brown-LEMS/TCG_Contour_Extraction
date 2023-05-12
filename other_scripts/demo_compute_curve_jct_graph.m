clear all; close all;
img_src_path = '/media/guoy/Research/Datasets/Drapery_Datasets/Cloth1/';
edge_src_path = '/media/guoy/Research/Datasets/Drapery_Datasets/Cloth1/GE/edges/';
frags_dst_path = '/media/guoy/Research/Datasets/Drapery_Datasets/Cloth1/GE/MSEL/';
mkdir(frags_dst_path);
img_ext = '*.png';

img_files = dir([img_src_path img_ext]);

for c = 1:length(img_files)
    
    disp([num2str(c) '/'  num2str(length(img_files))]); 
    
    %% Compute cfrags from edges: curvelets, unambiguous cfrags, HPG cfrags
    input = [edge_src_path img_files(c).name(1:end-4) '.edg'];
    output = [frags_dst_path img_files(c).name(1:end-4) '.cem'];
    
    cmd = ['!util/dborl_compute_curve_frags ' input ' ' output];
    eval(cmd);
    
    
    %% visualization
    [CEM, edges, cfrags_idx] = load_contours(output);
    G_test = construct_fac_graph_from_curve_fragments (cfrags_idx, CEM{2});

    H = figure(1);
    img = imread([img_src_path img_files(c).name]);
    imshow(img, 'border', 'tight'); hold on;
    draw_contours(CEM{2});
    for v =1:length(G_test.var)
        if(length(G_test.var(v).nbrs_fac)>=3)
            c_x = round(edges(G_test.var(v).actual_edge_id, 1)+1);
            c_y = round(edges(G_test.var(v).actual_edge_id, 2)+1);
            plot(c_x, c_y, 'gx');
        end
    end
    hold off;
    set(H, 'PaperPositionMode','auto')
    print(H,'-dpdf','-r0',[frags_dst_path img_files(c).name(1:end-4) '.pdf']);
    cmd = ['!pdfcrop ' frags_dst_path img_files(c).name(1:end-4) '.pdf ' frags_dst_path img_files(c).name(1:end-4) '.pdf'];
    eval(cmd);
end