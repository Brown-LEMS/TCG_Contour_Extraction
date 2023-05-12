clear all; close all;
addpath (genpath('util/'));

%% set paths
edge_src_path = 'Data/gPb_SEL_CFGD/edges/';
img_src_path = 'Data/CFGD_img/';
final_cem_path = 'Data/gPb_SEL_CFGD/final_curves/';
mkdir(final_cem_path);
cem_src_path = final_cem_path;
beta_src = 'training/';

%% load in beta
prefix = 'gPb_SEL_';

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
    [CEM, edges, cf_idx] = load_contours(input_file);
    [~, edgemap, thetamap] = load_edg([edge_src_path img_files(c).name(1:end-4) '.edg']);
    
    img = imread([img_src_path img_files(c).name]);
    [h,w,~]= size(img);

    
    % compute hsv space map
    hsv_img = rgb2hsv(img);
    
    % compute texton map

    tmap = assignTextons(fbRun(textonData.fb,rgb2gray(img)),textonData.tex);
    

    new_cfrags = CEM{2};
    
    con_cnt = length(new_cfrags);
    col_all = hsv(con_cnt);    % HSV colour map with con_cnt entries
    col_all = col_all(randperm(con_cnt),:); 
    
    imshow(img, 'border', 'tight'); hold on;
    draw_contours_sorted_col(new_cfrags, col_all);
    print_pdf([cem_src_path img_files(c).name(1:end-4) '.pdf']);
    hold off;
    
    %% prune curve fragments using trained logistic regression classifier
    disp('prune curve fragments using trained logistic regression classifier')
    P_vec = [];
    for i = 1:length(new_cfrags)
        [bg_grad, sat_grad, hue_grad, abs_k, edge_sparsity, wigg, len] = curve_fragment_cues(new_cfrags{i}, hsv_img, edgemap);

        p = 1 / (1 + exp(-([1, bg_grad, sat_grad, hue_grad, abs_k, edge_sparsity, wigg, len]-fmean_2)*beta_2'));
        P_vec = [P_vec p];
    end
    
    for j = 1:length(p_th_vec)
        
        dst_path = [final_cem_path 'p_' num2str(p_th_vec(j)*100) '/'];
        mkdir(dst_path);
        
        p_th = p_th_vec(j);
        disp(['  prob th: ' num2str(p_th)]);
        
        cfrags = cell(1,0);
        col = [];
        imshow(img, 'border', 'tight'); hold on;
        for i = 1:length(new_cfrags)
            if (P_vec(i) >p_th)
                cfrags = [cfrags new_cfrags{i}];
                col = [col; col_all(i,:)];
            end
        end
                
        
        draw_contours_sorted_col(cfrags, col);
%         write_cem([dst_path img_files(c).name(1:end-4) '.cem'], cfrags, h, w);
%         det_save_cemv([dst_path img_files(c).name(1:end-4) '.cemv'], cfrags);
        write_cem_fixed_edge_id([dst_path img_files(c).name(1:end-4) '.cem'], cfrags, edges, h, w);

        print_pdf([dst_path img_files(c).name(1:end-4) '.pdf']);
        hold off;
    end
    
    
    


%     keyboard;
    
end