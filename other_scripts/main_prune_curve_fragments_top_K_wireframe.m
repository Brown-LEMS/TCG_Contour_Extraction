clear all; close all;
addpath (genpath('../util/'));

%% set paths
top_K = 150;   
edge_src_path = '/media/guoy/Research/Project_contour/MSEL_contour_extraction_full/Data/wireframe/wireframe_SE/edges/';
img_src_path = '/media/guoy/Research/Project_contour/MSEL_contour_extraction_full/Data/wireframe/wireframe_SE/SE_TCG_wireframe_subset/';
cem_src_path = '/media/guoy/Research/Project_contour/MSEL_contour_extraction_full/Data/wireframe/wireframe_SE/SE_TCG_wireframe/';
final_cem_path = '/media/guoy/Research/Project_contour/MSEL_contour_extraction_full/Data/wireframe/wireframe_SE/SE_TCG_wireframe_subset/Top150/';
mkdir(final_cem_path);
beta_src = '../training/';

%% load in beta
prefix = 'SE_SEL_';

input_beta_2 = load([beta_src prefix 'beta_of_cues_for_seletion.txt']);
fmean_2 = input_beta_2(1,:);
fstd_2 = input_beta_2(2,:);
beta_2 = input_beta_2(3,:);
beta_2 = beta_2./fstd_2;



%% set parameters
p_th_vec = [0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.96 0.99];

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
img_files = dir([img_src_path '*.jpg']);

for c = 1:length(img_files)
    
    disp([num2str(c) '/'  num2str(length(img_files))]); 
    

    %% merge curve fragments using trained classifier 
    

    % load in curve fragment candidates
    input_file = [cem_src_path img_files(c).name(1:end-4) '.cem'];
    [CEM] = load_contours(input_file);
    [~, edgemap, thetamap] = load_edg([edge_src_path img_files(c).name(1:end-4) '.edg']);
    
    img = imread([img_src_path img_files(c).name]);
    [h,w,~]= size(img);

    
    % compute hsv space map
    hsv_img = rgb2hsv(img);
    
    % compute texton map

%     tmap = assignTextons(fbRun(textonData.fb,rgb2gray(img)),textonData.tex);
    

    new_cfrags = CEM{2};
    
%     con_cnt = length(new_cfrags);
%     col_all = hsv(con_cnt);    % HSV colour map with con_cnt entries
%     col_all = col_all(randperm(con_cnt),:); 
    
    
    
    %% prune curve fragments using trained logistic regression classifier
    disp('prune curve fragments using trained logistic regression classifier')
    P_vec = [];
    for j = 1:length( new_cfrags)
%         [bg_grad, sat_grad, hue_grad, abs_k, edge_sparsity, wigg, len, mean_conf] = curve_fragment_cues( new_cfrags{j}, hsv_img, edgemap);
% 
%         p = 1 / (1 + exp(-([1, bg_grad, sat_grad, hue_grad, abs_k, edge_sparsity, wigg, len, mean_conf]-fmean_2)*beta_2'));
       p = contour_length_mex(new_cfrags{j}');
       P_vec = [P_vec p];
    end 
    
    % rank the cfrags
    [~, sort_id] = sort(P_vec, 2, 'descend');
    new_cfrags =  new_cfrags(sort_id); 
    new_cfrags =  new_cfrags(1:min(top_K, length(new_cfrags)));
                
        
%         write_cem([dst_path img_files(c).name(1:end-4) '.cem'], cfrags, h, w);
    det_save_cemv([final_cem_path img_files(c).name(1:end-4) '.cemv'], new_cfrags);
%     write_cem_fixed_edge_id([final_cem_path img_files(c).name(1:end-4) '.cem'], cfrags, edges, h, w);

    H = figure(1);
    imshow(img, 'border', 'tight'); 
    hold on;
    draw_contours(new_cfrags);
    hold off;
    set(H, 'PaperPositionMode','auto')
    outimg = [final_cem_path img_files(c).name(1:end-4) '.pdf'];
    print(H,'-dpdf','-r0', outimg);
    cmd = ['!pdfcrop ' outimg ' ' outimg];
    eval(cmd);

    close all;
%     keyboard;
    
end