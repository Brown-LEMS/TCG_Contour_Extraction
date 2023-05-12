clear all; close all;
addpath (genpath('util/'));

src_path = 'Data/CUB_200_2011_TG/';
img_path = 'Data/CUB14_GT_FINAL/test/';
folder_path = dir(src_path);
k = 0;
for i = 1 : length(folder_path)
    if strcmp(folder_path(i).name, '.') | strcmp(folder_path(i).name, '..')| ...
            strcmp(folder_path(i).name, '.DS_Store')
        continue;
    else
        k = k + 1;
        subfolder{k} = folder_path(i).name;
    end
end

for g = 1:length(subfolder)
    
%% set paths
edge_src_path = [src_path subfolder{g} '/edges/'];
img_src_path = [img_path subfolder{g}(5:end) '/'];
final_cem_path = [src_path subfolder{g} '/test/'];
cem_src_path = [src_path subfolder{g} '/final_curves/'];

% edge_src_path = 'Data/TO_SEL_CFGD/edges/';
% img_src_path = 'Data/CFGD_img/';
% final_cem_path = 'Data/TO_SEL_CFGD/final_curves/';
rmdir(final_cem_path, 's')
mkdir(final_cem_path);
% cem_src_path = final_cem_path;
beta_src = 'training_for_bird_silhouette/';

%% load in beta
prefix = 'texture_silhou_';

input_beta_2 = load([beta_src prefix 'beta_of_cues_for_seletion.txt']);
fmean_2 = input_beta_2(1,:);
fstd_2 = input_beta_2(2,:);
beta_2 = input_beta_2(3,:);
beta_2 = beta_2./fstd_2;



%% set parameters
K = 40;
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
    %% prune curve fragments using trained logistic regression classifier
    disp('prune curve fragments using trained logistic regression classifier')
    P_vec = [];
    for i = 1:length(new_cfrags)
        [bg_grad, sat_grad, hue_grad, abs_k, edge_sparsity, wigg, len] = curve_fragment_cues(new_cfrags{i}, hsv_img, edgemap);

        p = 1 / (1 + exp(-([1, bg_grad, sat_grad, hue_grad, abs_k, edge_sparsity, wigg, len]-fmean_2)*beta_2'));
        P_vec = [P_vec p];
    end
    
    [~,p_index] =  sort(P_vec, 2, 'descend');
    new_cfrags = new_cfrags(p_index);
        
    dst_path = [final_cem_path 'p_top_K/'];
    mkdir(dst_path);


    cfrags = cell(1,0);
    imshow(img, 'border', 'tight'); hold on;
    for i = 1:K
            cfrags = [cfrags new_cfrags{i}];
    end

    draw_contours(cfrags);
    if ~isempty(cfrags)
    write_cem([dst_path img_files(c).name(1:end-4) '.cem'], cfrags, h, w);
%         det_save_cemv([dst_path img_files(c).name(1:end-4) '.cemv'], cfrags);
    print_pdf([dst_path img_files(c).name(1:end-4) '.pdf']);
    end
    hold off;

%     keyboard;
    
end
end
