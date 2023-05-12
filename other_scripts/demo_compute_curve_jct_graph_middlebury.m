clear all; close all;
addpath (genpath('util'))
img_src_path = '/media/guoy/Research/Datasets/MiddEval/trainingQ/';
edge_src_path = '/media/guoy/Research/Datasets/MiddEval/trainingQ_SE/';
frags_dst_path = '/media/guoy/Research/Datasets/MiddEval/trainingQ_SE/';
mkdir(frags_dst_path);

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
    
    %% Compute cfrags from edges: curvelets, unambiguous cfrags, HPG cfrags
    input = [edge_src_path image_names{c} '/im0.edg'];
    output = [frags_dst_path image_names{c} '/im0.cem'];
    
    cmd = ['!util/dborl_compute_curve_frags ' input ' ' output];
%     cmd = ['!util/dborl_compute_curve_frags_pos_uncertain05 ' input ' ' output];
    eval(cmd);
    
    
    %% visualization
    [CEM, edges, cfrags_idx] = load_contours(output);
    G_test = construct_fac_graph_from_curve_fragments (cfrags_idx, CEM{2});

    H = figure(1);
    img = imread([img_src_path image_names{c} '/im0.png']);
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
    print(H,'-dpdf','-r0',[frags_dst_path image_names{c} 'im0.pdf']);
    cmd = ['!pdfcrop ' frags_dst_path image_names{c} 'im0.pdf ' frags_dst_path image_names{c} 'im0.pdf'];
    eval(cmd);
    
    %% Compute cfrags from edges: curvelets, unambiguous cfrags, HPG cfrags
    input = [edge_src_path image_names{c} '/im1.edg'];
    output = [frags_dst_path image_names{c} '/im1.cem'];
    
    cmd = ['!util/dborl_compute_curve_frags ' input ' ' output];
%     cmd = ['!util/dborl_compute_curve_frags_pos_uncertain05 ' input ' ' output];
    eval(cmd);
    
    
    %% visualization
    [CEM, edges, cfrags_idx] = load_contours(output);
    G_test = construct_fac_graph_from_curve_fragments (cfrags_idx, CEM{2});
    H = figure(1);
    img = imread([img_src_path image_names{c} '/im1.png']);
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
    print(H,'-dpdf','-r0',[frags_dst_path image_names{c} 'im1.pdf']);
    cmd = ['!pdfcrop ' frags_dst_path image_names{c} 'im1.pdf ' frags_dst_path image_names{c} 'im1.pdf'];
    eval(cmd);
    
end