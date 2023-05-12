% Demo for Structured Edge Detector (please see readme.txt first).
addpath(genpath('toolbox-master/'))
addpath(genpath('../'))
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
th = 0.05;
sigma = 2;
src_path = '/gpfs/scratch/yg13/Datasets/KITTI/data_stereo_flow/training/colored_1/';
img_files = dir([src_path '*.png']);
out_path = '/gpfs/scratch/yg13/Datasets/KITTI/data_stereo_flow/training/colored_1_SE/';
mkdir(out_path);
for i = 1:length(img_files)
    i
    %% img 0
    imgFile = [src_path img_files(i).name];
    outFile = [out_path img_files(i).name(1:end-4) '.edg'];
    I = imread(imgFile);
    tic,[E,O,edginfo,inds,segs] = edgesDetect_TO(I,model, th, sigma); toc
    
    outImg = [out_path img_files(i).name(1:end-4) '_SE.png'];
    imwrite(1-E, outImg, 'PNG');
    
    % change edges to cxx coodinates
    edginfo(:,1) = edginfo(:,1)-1;
    edginfo(:,2) = edginfo(:,2)-1;
    save_edg(outFile, edginfo, [size(I,2) size(I,1)]);
   
end