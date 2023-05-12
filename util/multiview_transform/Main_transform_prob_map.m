%Load the CEM file
addpath(genpath('/media/guoy/Research/Project_contour/MSEL_contour_extraction_full/util'))
img_path = '/media/guoy/Research/Datasets/Drapery_Datasets/Cloth1/';
bry_path = '/media/guoy/Research/Datasets/Drapery_Datasets/Cloth1_COB/';


for i = 1:23
    figure(1)
    id0 = i-1;
    id1 = i;
    
    im0 = imread([img_path 'Camera' num2str(id0) '.png']);
    [h, w, ~] = size(im0);
    bry_map0 = imread([bry_path 'Camera' num2str(id0) '.png']);
    bry_map0 = imresize(bry_map0, [h,w]);
    figure(1)
    imshow(bry_map0, 'border', 'tight');
    bry_map0 = 1-im2double(bry_map0);

    %Get Corresponding Ground Truth Contours
    %Arguments: Image ID1 
    %			Image ID2
    %			CEM data just loaded
    %			The path of dataset
    %Output: A Cell with all the corresponding curves in the second image.
    bry_map_trans = transform_bry_map_given_GT_depth(id0,id1, bry_map0, img_path);
    bry_map_trans = 1- bry_map_trans;
    
    H=figure(2);
    imshow(bry_map_trans, 'border', 'tight');
    
    set(H, 'PaperPositionMode','auto')
    out_name = [bry_path 'Camera_' num2str(id0) 'to' num2str(id1) '.pdf'];
    print(H,'-dpdf','-r0', out_name);
    cmd = ['!pdfcrop ' out_name ' ' out_name];
    eval(cmd);

    
end