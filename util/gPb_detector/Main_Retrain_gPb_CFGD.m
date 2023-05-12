addpath(genpath('../../third_order'));
addpath(fullfile(pwd,'lib'));
clear all; close all;

n=1000000;
buffer=2;
load('beta_mPb_CFGD.txt');
mbeta = beta_mPb_CFGD(2,:)./beta_mPb_CFGD(1,:);
dim =1;

pres = 'gPb_CFGD';
train_path = 'CFGD_train/';
img_path = 'CFGD/';
% test_path = 'test_4/';
img_files = dir([img_path,'*.jpg']);

% number of samples per image
nPer = ceil(n/length(img_files));

y = zeros(1,0);
f = [];


for i = 1:length(img_files),
    tic;
    % read the image
    iid = img_files(i).name;
    fprintf(2,'Processing image %d/%d (iid=%s)...\n',i,length(img_files),iid);
    %   im = double(imread(imgFilename(iid))) / 255;
    im = double(imread([img_path, iid])) / 255;
    load([train_path iid(1:end-4) '.mat']);    
    % run the detector to get feature vectors
    fprintf(2,'  Running detector...\n');
    features = feval(@gPbdetector_new,im, mbeta);
    % load the segmentations and union the boundary maps
    bmap = zeros(size(groundTruth{1}.Boundaries(:,:,1)));
    bmap = bmap | groundTruth{1}.Boundaries(:,:,1);
    bmap = bmap | groundTruth{2}.Boundaries(:,:,1);
    bmap = bmap | groundTruth{3}.Boundaries(:,:,1);

    %   bmap = logical(imread([train_path, iid(1:end-4), '.png']));
    bmap = bmap(:,:,1);
    dmap = bwdist(bmap);

    % sample 
    fprintf(2,'  Sampling...\n');
    onidx = find(bmap)';
    offidx = find(dmap>buffer)';
    ind = [ onidx offidx ];

    % here when the positive potion is small, use equal number of negative samples 
    %   cnt = numel(ind);
    %   idx = randperm(cnt);
    %   idx = idx(1:min(cnt,nPer));
    cnt = numel(offidx);
    idx = randperm(cnt);
    offidx = offidx(idx);

    ind = [ onidx offidx ];

    cnt2 = 2*min(numel(onidx), numel(offidx));
    idx = randperm(cnt2);


    y = [ y bmap(ind(idx)) ];
    f = [ f features(:,ind(idx)) ];
    fprintf(2,'  %d samples.\n',numel(idx));
    toc;
end

save(['features_' pres, '.mat'], 'y', 'f');
load(['features_' pres, '.mat'], 'y', 'f');

f=f'; y=y';
% normalize features to unit variance
fstd = std(f);
fstd = fstd + (fstd==0);
f = f ./ repmat(fstd,size(f,1),1);

if(dim==1)
    f=f(:, [1:4, 11:end]);
    std(5:10)=1;
end
% fit the model
fprintf(2,'Fitting model...\n');
beta = logist2(y,f);
beta = beta';

if(dim==1)
   beta_new = zeros(1,14);
   beta_new(1:4) = beta(1:4);
   beta_new(11:14) = beta(5:8);
   beta = beta_new;
end
beta
% save the result
save_name = ['beta_', pres, '.txt'];
save(save_name,'fstd','beta','-ascii');


%% testing
% saved_beta = load(save_name);
% fstd = saved_beta(1,:);
% beta = saved_beta(2,:);
% beta = beta ./ fstd;
% 
% file=dir([test_path,'*.jpg']);
% 
% for j=1:length(file)
%     name=file(j).name;
%     imgFile=[test_path,name];
%     fprintf(['filename ', imgFile, '\n']);
%     outFile=[test_path,name(1:(end-4)), '_TO_pb.edg'];
%     im = imread (imgFile);
%     im = im2double(im);
%     edg = new_pb_subpixel_TO(im, 0.1, 4, beta);
%     save_edg(outFile, edg, [size(im,2) size(im,1)]);
% 
% end



