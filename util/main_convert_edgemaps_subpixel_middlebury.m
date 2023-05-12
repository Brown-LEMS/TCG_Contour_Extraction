clear all; close all;
addpath(genpath('./'))
% addpath('SE_Dollar_detector/private/')
img_path = '/media/guoy/Research/Datasets/MiddEval/trainingQ/';
src_path = '/media/guoy/Research/Datasets/MiddEval/trainingQ_COB/';
out_path = '/media/guoy/Research/Datasets/MiddEval/trainingQ_COB/';
mkdir(out_path);

input_files = dir([src_path '*.mat']);
sigma = 3;
th = 0.1;

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

for i = 1:length(image_names)
    
    %%%%%%%%%%%%%%%%%%%%%%%%  im0  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    img = imread([img_path image_names{i} '/im0.png']);
    [h,w,~] = size(img);
    load([src_path image_names{i} '/im0.mat']);
    O_input = O;
    
    E = imresize(E{1}, [h,w]);
    outImg = [out_path image_names{i} '/im0_bry.png'];
    imwrite(1-E, outImg, 'PNG');
    O = imresize(single(O.angle), [h,w]);
    E=convTri(single(E),3);
%     [Ox,Oy]=gradient2(convTri(E,4));
%     [Oxx,~]=gradient2(Ox); [Oxy,Oyy]=gradient2(Oy);
%     O=mod(atan(Oyy.*sign(-Oxy)./(Oxx+1e-5)),pi);
    
    
    E_thin = edgesNmsMex(E,O,2,5,1.01, 1);
%     imshow(E_thin);
    [ TO_edgemap, TO_orientation, edginfo ] = subpix_TO_correction( E_thin,O,th,sigma, 0);
    figure(1)
    imshow(1-TO_edgemap)
    outImg = [out_path image_names{i} '/im0_bry_thin.png'];
    imwrite(1-TO_edgemap, outImg, 'PNG');
    % change edges to cxx coodinates
%     edginfo(:,3) = O_input.angle(sub2ind([h,w], round(edginfo(:,2)), round(edginfo(:,1))))+pi/2;
%     edginfo(:,3) =  wrapToPi(edginfo(:,3));
    edginfo(:,1) = edginfo(:,1)-1;
    edginfo(:,2) = edginfo(:,2)-1;
    outFile = [out_path image_names{i} '/im0.edg'];
    save_edg(outFile, edginfo, [w h]);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%  im1  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    img = imread([img_path image_names{i} '/im1.png']);
    [h,w,~] = size(img);
    load([src_path image_names{i} '/im1.mat']);
    O_input = O;
    
    E = imresize(E{1}, [h,w]);
    outImg = [out_path image_names{i} '/im1_bry.png'];
    imwrite(1-E, outImg, 'PNG');
    O = imresize(single(O.angle), [h,w]);
    E=convTri(single(E),3);
%     [Ox,Oy]=gradient2(convTri(E,4));
%     [Oxx,~]=gradient2(Ox); [Oxy,Oyy]=gradient2(Oy);
%     O=mod(atan(Oyy.*sign(-Oxy)./(Oxx+1e-5)),pi);
    
    
    E_thin = edgesNmsMex(E,O,2,5,1.01, 1);
%     imshow(E_thin);
    [ TO_edgemap, TO_orientation, edginfo ] = subpix_TO_correction( E_thin,O,th,sigma, 0);
    figure(2)
    imshow(1-TO_edgemap)
    outImg = [out_path image_names{i} '/im1_bry_thin.png'];
    imwrite(1-TO_edgemap, outImg, 'PNG');
    % change edges to cxx coodinates
%     edginfo(:,3) = O_input.angle(sub2ind([h,w], round(edginfo(:,2)), round(edginfo(:,1))))+pi/2;
%     edginfo(:,3) =  wrapToPi(edginfo(:,3));
    edginfo(:,1) = edginfo(:,1)-1;
    edginfo(:,2) = edginfo(:,2)-1;
    outFile = [out_path image_names{i} '/im1.edg'];
    save_edg(outFile, edginfo, [w h]);
    
end