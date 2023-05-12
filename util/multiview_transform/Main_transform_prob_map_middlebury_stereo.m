%Load the CEM file
addpath(genpath('/media/guoy/Research/Project_contour/MSEL_contour_extraction_full/util'))
img_path = '/media/guoy/Research/Datasets/MiddEval/trainingQ/';
bry_path = '/media/guoy/Research/Datasets/MiddEval/trainingQ_COB/';

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
    
    %%%%%%%%%%%%% load ucm maps and GT disparity
    bry_map0 = imread([bry_path image_names{i} '/im0.png']);
    im0 = imread([img_path image_names{i} '/im0.png']);
    bry_map1 = imread([bry_path image_names{i} '/im1.png']);
%     im1 = imread([img_path image_names{i} '/im1.png']);
    
    [h, w, ~] = size(im0);
    bry_map0 = imresize(bry_map0, [h,w]);
    bry_map0 = 1-im2double(bry_map0);

    bry_map1 = imresize(bry_map1, [h,w]);
    bry_map1 = 1-im2double(bry_map1);
    
    GT_disparity0 = readpfm([img_path image_names{i} '/disp0GT.pfm']);
    GT_disparity1 = readpfm([img_path image_names{i} '/disp1GT.pfm']);

    %%%%%%%%%%%%% transform ucm maps given disparity
    bry_map_trans01 = transform_bry_map_given_GT_disparity( bry_map0, GT_disparity0, true);
    bry_map_trans10 = transform_bry_map_given_GT_disparity( bry_map1, GT_disparity1, false);
    
    %%%%%%%%%%%% visulize im1 to im0

    H=figure(1);
    imshow(1-bry_map_trans01);
    set(H, 'PaperPositionMode','auto')
    out_name = [bry_path image_names{i} 'im_0to1.pdf'];
    print(H,'-dpdf','-r0', out_name);
    cmd = ['!pdfcrop ' out_name ' ' out_name];
    eval(cmd);

    
    %%%%%%%%%%%% visulize im1 to im0
    H=figure(2);
    imshow(1-bry_map_trans10);
    set(H, 'PaperPositionMode','auto')
    out_name = [bry_path image_names{i} 'im_1to0.pdf'];
    print(H,'-dpdf','-r0', out_name);
    cmd = ['!pdfcrop ' out_name ' ' out_name];
    eval(cmd);
end