% Demo for Structured Edge Detector (please see readme.txt first).
addpath(genpath('toolbox-master/'))
addpath(genpath('/media/guoy/Research/Project_contour/MSEL_contour_extraction_full/util'))
addpath(genpath('/media/guoy/Research/Project_contour/COB/'))

%% set opts for training (see edgesTrain.m)
opts=edgesTrain();                % default options (good settings)
opts.modelDir='models/';          % model will be in models/forest
opts.modelFnm='modelBsds';        % model name
opts.nPos=5e5; opts.nNeg=5e5;     % decrease to speedup training
opts.useParfor=0;                 % parallelize if sufficient memory

%% train edge detector (~20m/8Gb per tree, proportional to nPos/nNeg)
tic, model=edgesTrain(opts); toc; % will load model if already trained

%% set detection parameters (can set after training)
model.opts.multiscale=0;          % for top accuracy set multiscale=1
model.opts.sharpen=2;             % for top speed set sharpen=0
model.opts.nTreesEval=4;          % for top speed set nTreesEval=1
model.opts.nThreads=4;            % max number threads for evaluation
model.opts.nms=true;                 % set to true to enable nms

%% evaluate edge detector on BSDS500 (see edgesEval.m)
if(0), edgesEval( model, 'show',1, 'name','' ); end

%% detect edge and visualize results
th = 0.07;
sigma = 1;
% src_path = '/media/guoy/Research/Project_contour/MSEL_contour_extraction_full/Data/CFGD/img/';
% src_path = '/media/guoy/Research/Datasets/3D_Drawing/pavilion-midday-mcs-work/';
% src_path = '/media/New_Volume/Research/Datasets/capital_dataset_resized/jpeg/';
% src_path = '/media/New_Volume/Research/Datasets/Illumination_selected/SET003/';
% src_path = '/media/guoy/Research/Datasets/random_noise/sigma5/';
% src_path = '/media/guoy/Research/Datasets/Drapery_Datasets/Cloth1/';
% src_path = 'Data/CUB14/all_imgs/';
src_path = '/media/guoy/Research/Project_contour/MSEL_contour_extraction_full/Data/trifocal/img/';
img_files = dir([src_path '*.png']);
% out_path = 'results/CUB14/';
% out_path = 'results/capital_dataset_resized_test/';
% out_path = 'results/Illumination_selected/SET003/';
% out_path = '/media/guoy/Research/Datasets/Drapery_Datasets/Cloth1_SE/edges/';
% out_path = '/media/guoy/Research/Project_contour/MSEL_contour_extraction_full/Data/CFGD/CFGD_SE/edges/';
% out_path = '/media/guoy/Research/Datasets/3D_Drawing/pavilion-midday-mcs-work/SE_edges/';
out_path = '/media/guoy/Research/Project_contour/MSEL_contour_extraction_full/Data/trifocal/SE_edges/';
mkdir(out_path);
for i = 1:length(img_files)
    i
    imgFile = [src_path img_files(i).name];
    outFile = [out_path img_files(i).name(1:end-4) '.edg'];
    outFile2 = [out_path img_files(i).name(1:end-4) '.mat'];
    I = imread(imgFile);
%     tic, [E,O,~,~]=edgesDetect(I,model); toc
    tic,[E,~,edginfo,inds,segs, E_org, O_org] = edgesDetect_TO(I,model, th, sigma); toc

    % change edges to cxx coodinates
    edginfo(:,1) = edginfo(:,1)-1;
    edginfo(:,2) = edginfo(:,2)-1;
    save_edg(outFile, edginfo, [size(I,2) size(I,1)]);
        
    outImg = [out_path img_files(i).name(1:end-4) '_SE.png'];
    imwrite(1-E, outImg, 'PNG');
    outImg = [out_path img_files(i).name(1:end-4) '_bry.png'];
    imwrite(1-E_org, outImg, 'PNG');
    
    % compute ucm
    outImg = [out_path img_files(i).name(1:end-4) '_ucm.png'];
    ucms = edge2ucm(E_org,O_org);
    ucms = ucms/max(ucms(:));
    imwrite(1-ucms, outImg, 'PNG');
    
    E = E_org;
    O.angle = O_org;
    save(outFile2, 'ucms', 'E', 'O');
   
end