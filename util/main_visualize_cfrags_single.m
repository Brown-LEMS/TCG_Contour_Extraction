clear all; close all;
addpath (genpath('./'));

src_path = '../Data/TO_Kovesi_CFGD/';
edge_src_path = '../Data/TO_SEL_CFGD/edges/';
image_path = '../Data/CFGD_img/';
output_path = '/media/guoy/Research/My_Papers/Yuliang_Guo_working_paper_contour/figs/visualization_problem/';
% mkdir(output_path);
% beta_src = '../training/';
% prefix = 'TO_SEL_';
% top_K = 100;
% 
% 
% input_beta_2 = load([beta_src prefix 'beta_of_cues_for_seletion.txt']);
% fmean_2 = input_beta_2(1,:);
% fstd_2 = input_beta_2(2,:);
% beta_2 = input_beta_2(3,:);
% beta_2 = beta_2./fstd_2;

files = dir([output_path '*10081.jpg']);
for i = 1:length(files)
    i
    name=[files(i).name];
    src=[src_path,name(1:end-4), '.cem'];
    [CEM, edges, cf_idx]=load_contours(src);
    I = imread([image_path, name(1:end-4) '.jpg']);
    [~, edgemap, thetamap] = load_edg([edge_src_path name(1:end-4) '.edg']);
    output_file1 = [output_path name(1:end-4), 'TO.png'];
    output_file2 = [output_path name(1:end-4), 'TO_Kovesi.pdf'];

    
    h1 = figure(1);
    imshow(1-edgemap,'border', 'tight');
    set(h1, 'PaperPositionMode','auto')
    print(h1,'-dpng',output_file1);
    
    h2 = figure(2);
    imshow(I,'border', 'tight'); hold on;
% %     load([src(1:end-4) '.mat']);
% %     [edg, edgemap, thetamap] = load_edg([edge_src_path, name(1:end-4), '.edg']);
% %     I = zeros(d{1}(2),d{1}(1),3);
    
%     hsv_img = rgb2hsv(I);
% 
% 
%     imshow(rgb2gray(I),'border', 'tight');hold on;
% 
%     
% %%%%%%%%%%%%% Rank of results cfrags and prune using logistic regressor        
%     disp('rank curve fragments');
%     new_cfrags = CEM{2};
%     P_vec = [];
%     for c = 1:length(new_cfrags)
% %         [bg_grad, sat_grad, hue_grad, abs_k, edge_sparsity, wigg, len, mean_conf] = curve_fragment_cues(new_cfrags{i}, hsv_img, edgemap);
% % 
% %         p = 1 / (1 + exp(-([1, bg_grad, sat_grad, hue_grad, abs_k, edge_sparsity, wigg, len, mean_conf]-fmean_2)*beta_2'));
% %         
%         p = contour_length_mex(new_cfrags{c}');
%         P_vec = [P_vec p];
%     end 
%     
%     % rank the cfrags
%     [~, sort_id] = sort(P_vec, 2, 'descend');
%     new_cfrags = new_cfrags(sort_id); 
%     
%     top_cfrags = new_cfrags(1: min(top_K, length(new_cfrags)));
%             
%     
%     draw_contours(top_cfrags);
    
    
    draw_contours(CEM{2});
    
    hold off;
    set(h2, 'PaperPositionMode','auto')
    print(h2,'-dpdf',output_file2);
    cmd = ['!pdfcrop ' output_file2 ' ' output_file2];
    eval(cmd);
end