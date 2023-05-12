clear all; close all;

addpath(genpath('/media/guoy/Research/Project_contour/MSEL_contour_extraction_full/util'))
img_path = '/media/guoy/Research/Datasets/MiddEval/trainingQ/';
out_path = '/media/guoy/Research/Datasets/MiddEval/trainingQ/';
mkdir(out_path);

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

    i
    
    GT_disparity0 = readpfm([img_path image_names{i} '/disp0GT.pfm']);
    GT_disparity1 = readpfm([img_path image_names{i} '/disp1GT.pfm']);
    GT_BW0 = edge(GT_disparity0, 'canny');
    GT_BW1 = edge(GT_disparity1, 'canny');
    
    outfile0 = [out_path image_names{i} '/im_0_disp_bry.png'];
    imwrite(GT_BW0, outfile0, 'PNG');
    outfile1 = [out_path image_names{i} '/im_1_disp_bry.png'];
    imwrite(GT_BW1, outfile1, 'PNG');    
end