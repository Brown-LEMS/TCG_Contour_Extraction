% Demo for Structured Edge Detector (please see readme.txt first).
addpath(genpath('toolbox-master/'))
addpath(genpath('/media/New_Volume/Research/Yuliang_Tracking_Toolbox'))
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

edg_ext = '_TO_edges.png';
orient_ext = '_TO_orient.txt';
th = 0.03;
src_path = '/media/New_Volume/Research/BMS_mice_pose_tracking/Seqs_full_length/';
mat_files = dir([src_path '*bg.mat']);
% img_files = dir([src_path '*.png']);
% out_path = 'results/CUB14/';
% out_path = '/gpfs/scratch/yg13/cornell/101815M1-Social/edge_org/';
out_path = src_path;
mkdir(out_path);
for i = 1:length(mat_files)
    i
    bg_file = [src_path mat_files(i).name]
    load(bg_file);
%     outFile = [out_path img_files(i).name(1:end-4) '.edg'];
    if(size(bg,3)==1)
        th = 0.03;
        I = uint8(repmat(bg, [1,1,3]));
    else
        th = 0.05;
        I = uint8(bg);
    end
%     tic, [E,O,~,~]=edgesDetect(I,model); toc
    tic,[E,O,TO_edges,inds,segs] = edgesDetect_TO(I,model,th,2); toc
%     [xx,yy,mag] = find(E~=0);
%     Angle = zeros(length(xx),1); 
%     for j=1:length(xx)
%         Angle(j) = O(xx(j),yy(j));
%     end
%     edge_map  = [xx yy Angle mag];
%     figure(1); im(I); figure(2); im(1-E);
    
%     outImg = [out_path img_files(i).name(1:end-4) '.png'];
%     imwrite(E>th, outImg, 'PNG');
    
   
    
    edgemap = zeros(size(I,1), size(I,2));
    thetamap = zeros(size(I,1), size(I,2));
    [m,k] = size(TO_edges);
    
    for k = 1:m
        x = int32(TO_edges(k,1))+1;
        y = int32(TO_edges(k,2))+1;
        edgemap(y, x) = 255;
        thetamap(y,x) = TO_edges(k,3);
    end
%     imshow(edgemap)
    
    dst_edg_file = [bg_file(1:end-4) edg_ext];
    dst_ori_file = [bg_file(1:end-4) orient_ext];
    
    imwrite(edgemap, dst_edg_file);
    dlmwrite(dst_ori_file, thetamap, ' ');
    
    TO_edges(:,1) = TO_edges(:,1)-1;
    TO_edges(:,2) = TO_edges(:,2)-1;
    save_edg([bg_file(1:end-4) '.edg'], TO_edges, [size(I,2) size(I,1)]);
    
%     keyboard;
end