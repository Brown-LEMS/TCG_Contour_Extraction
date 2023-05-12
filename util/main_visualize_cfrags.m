clear all; close all;
addpath (genpath('./'));

src_path = '../Data/TO_Kokkinos_VOC2007/';
edge_src_path = '../Data/TO_SEL_VOC2007/edges/';
image_path = '../Data/VOC2007_img/';
output_path = '/media/New_Volume/Research/My_Papers/Yuliang_Guo_working_paper_contour/figs/visualization_voc/TO_Kokkinos/';
mkdir(output_path);
beta_src = '../training/';
prefix = 'TO_SEL_';
top_K = 100;


input_beta_2 = load([beta_src prefix 'beta_of_cues_for_seletion.txt']);
fmean_2 = input_beta_2(1,:);
fstd_2 = input_beta_2(2,:);
beta_2 = input_beta_2(3,:);
beta_2 = beta_2./fstd_2;

files = dir([output_path '*png']);
for i = 1:length(files)
    i
    name=[files(i).name];
    src=[src_path,name(1:end-4), '.cem'];
    [CEM, edges, cf_idx]=load_contours(src);
%     load([src(1:end-4) '.mat']);
%     [edg, edgemap, thetamap] = load_edg([edge_path, name(1:end-4), '.edg']);
%     I = zeros(d{1}(2),d{1}(1),3);
    I = imread([image_path, name(1:end-4) '.jpg']);
    hsv_img = rgb2hsv(I);
    [~, edgemap, thetamap] = load_edg([edge_src_path name(1:end-4) '.edg']);


%     imshow(I,'border', 'tight'); hold on;
    imshow(rgb2gray(I),'border', 'tight');hold on;

    
%%%%%%%%%%%%% Rank of results cfrags and prune using logistic regressor        
    disp('rank curve fragments');
    new_cfrags = CEM{2};
    P_vec = [];
    for c = 1:length(new_cfrags)
%         [bg_grad, sat_grad, hue_grad, abs_k, edge_sparsity, wigg, len, mean_conf] = curve_fragment_cues(new_cfrags{i}, hsv_img, edgemap);
% 
%         p = 1 / (1 + exp(-([1, bg_grad, sat_grad, hue_grad, abs_k, edge_sparsity, wigg, len, mean_conf]-fmean_2)*beta_2'));
%         
        p = contour_length_mex(new_cfrags{c}');
        P_vec = [P_vec p];
    end 
    
    % rank the cfrags
    [~, sort_id] = sort(P_vec, 2, 'descend');
    new_cfrags = new_cfrags(sort_id); 
    
    top_cfrags = new_cfrags(1: min(top_K, length(new_cfrags)));
            
    
    draw_contours(top_cfrags);
    
% % write cem has problem that save identical edges multiple times so that the graph structure is not saved    
%     G = construct_fac_graph_from_curve_fragments (edges, cf_idx, CEM{2,1});
%     
%     % only ploting the connecting nodes
%     for v = 1:length(G.var)
%         if(G.var(v).dim >=2)
%             actual_edge = edges(G.var(v).actual_edge_id,:);
%             plot(actual_edge(1,1)+1, actual_edge+1, 'y.');
%         end
%     end
    
    hold off;
%     disp_edg(edg);
    h = figure(1);
%     keyboard;
    print_pdf([output_path name(1:end-4), '.pdf'], h);
%     set(h, 'PaperPositionMode','auto')
%     print(h,'-djpeg95','-r0',[src(1:end-4), '_color.jpg']);
end