% Demo for Structured Edge Detector (please see readme.txt first).
addpath(genpath('toolbox-master/'))
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
% src_path = 'Data/BSDS500/data/images/test/';
src_path = '/media/New_Volume/Research/Project_contour/SEL_release/Data/CFGD_img/';
img_files = dir([src_path '*.jpg']);
out_path = '/media/New_Volume/Research/Project_contour/SEL_release/Data/SE_no_TO_SEL_CFGD/edges/';
mkdir(out_path);

th = 0.07;
for i = 1:length(img_files)
    i
    imgFile = [src_path img_files(i).name(1:end-4) '.jpg'];
    outFile = [out_path img_files(i).name(1:end-4) '.edg'];
    I = imread(imgFile);
    tic, [E,O,~,~] = edgesDetect(I,model); toc
    
    [h,w,~] = size(I);
    [Y,X] = find(E>th);
    
    mag = E(sub2ind([h,w], Y,X));
    Angle = O(sub2ind([h,w], Y,X))+pi/2;
    Angle = wrapTo2Pi(Angle);
    
    edge_map  = [X-1 Y-1 Angle mag];
    
    save_edg(outFile, edge_map, [w,h]);
%     keyboard
%     figure(1); im(I); figure(2); im(1-E);
%     outImg = [out_path img_files(i).name(1:end-4) '.png'];
%     imwrite(E, outImg, 'PNG');
   
end