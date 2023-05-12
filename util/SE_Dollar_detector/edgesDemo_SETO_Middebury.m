% Demo for Structured Edge Detector (please see readme.txt first).
addpath(genpath('toolbox-master/'))
% addpath('/media/guoy/Research/Print_PDF/')
addpath(genpath('/media/guoy/Research/Project_contour/MSEL_contour_extraction_full/util'))
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
th = 0.03;
sigma = 1.5;
src_path = '/media/guoy/Research/Datasets/MiddEval/trainingQ/';
% img_files = dir([src_path '*.png']);
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

out_path = '/media/guoy/Research/Datasets/MiddEval/trainingQ_SE/';
mkdir(out_path);
for i = 1:length(image_names)
    i
    mkdir([out_path image_names{i}])
    %% img 0
    imgFile = [src_path image_names{i} '/im0.png'];
    outFile = [out_path image_names{i} '/im0.edg'];
    I = imread(imgFile);
    tic,[E,O,edginfo,inds,segs, E_org] = edgesDetect_TO(I,model, th, sigma); toc
    
    outImg = [out_path image_names{i} '/im0_SE.png'];
    imwrite(1-E, outImg, 'PNG');
    outImg = [out_path image_names{i} '/im0_bry.png'];
    imwrite(1-E_org, outImg, 'PNG');
    
    % change edges to cxx coodinates
    edginfo(:,1) = edginfo(:,1)-1;
    edginfo(:,2) = edginfo(:,2)-1;
    save_edg(outFile, edginfo, [size(I,2) size(I,1)]);
    %% img 1
    imgFile = [src_path image_names{i} '/im1.png'];
    outFile = [out_path image_names{i} '/im1.edg'];
    I = imread(imgFile);
    tic,[E,O,edginfo,inds,segs, E_org] = edgesDetect_TO(I,model, th, sigma); toc
    
    outImg = [out_path image_names{i} '/im1_SE.png'];
    imwrite(1-E, outImg, 'PNG');

    outImg = [out_path image_names{i} '/im1_bry.png'];
    imwrite(1-E_org, outImg, 'PNG');
    
    % change edges to cxx coodinates
    edginfo(:,1) = edginfo(:,1)-1;
    edginfo(:,2) = edginfo(:,2)-1;
    save_edg(outFile, edginfo, [size(I,2) size(I,1)]);
   
end