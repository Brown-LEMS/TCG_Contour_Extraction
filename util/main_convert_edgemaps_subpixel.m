clear all; close all;
addpath(genpath('./'))
% addpath('SE_Dollar_detector/private/')
img_path = '/media/guoy/Research/Datasets/Drapery_Datasets/Cloth1/';
src_path = '/media/guoy/Research/Datasets/Drapery_Datasets/Cloth1_COB/';
out_path = '/media/guoy/Research/Datasets/Drapery_Datasets/Cloth1_COB/edges/';
mkdir(out_path);

input_files = dir([src_path '*.mat']);
sigma = 2;
th = 0.05;

for i = 1:length(input_files)
    
    img = imread([img_path input_files(i).name(1:end-4) '.png']);
    [h,w,~] = size(img);
    load([src_path input_files(i).name]);
    O_input = O;
    
    E = imresize(E{1}, [h,w]);
    outImg0 = [out_path input_files(i).name(1:end-4) '_bry.png'];
    imwrite(1-E, outImg0, 'PNG');
    O = imresize(single(O.angle), [h,w]);
    E=convTri(single(E),1);
%     [Ox,Oy]=gradient2(convTri(E,4));
%     [Oxx,~]=gradient2(Ox); [Oxy,Oyy]=gradient2(Oy);
%     O=mod(atan(Oyy.*sign(-Oxy)./(Oxx+1e-5)),pi);
    
    E_thin = edgesNmsMex(E,O,2,5,1.01, 1);
%     imshow(E_thin);
    [ TO_edgemap, TO_orientation, edginfo ] = subpix_TO_correction( E_thin,O,th,sigma, 0);
    imshow(1-TO_edgemap)
    outImg1 = [out_path input_files(i).name(1:end-4) '_bry_thin.png'];
    imwrite(1-TO_edgemap, outImg1, 'PNG');
    % change edges to cxx coodinates
%     edginfo(:,3) = O_input.angle(sub2ind([h,w], round(edginfo(:,2)), round(edginfo(:,1))))+pi/2;
%     edginfo(:,3) =  wrapToPi(edginfo(:,3));
    edginfo(:,1) = edginfo(:,1)-1;
    edginfo(:,2) = edginfo(:,2)-1;
    outFile = [out_path input_files(i).name(1:end-4) '.edg'];
    save_edg(outFile, edginfo, [w h]);
end