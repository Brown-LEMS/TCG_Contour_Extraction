clear all; close all;

voc_path = '/media/yuliangguo/NewVolume_1/VOCdevkit/VOC2007/';
voc_seg_img_list = textread([voc_path 'ImageSets/Segmentation/trainval.txt'], '%s');
voc_image_path = [voc_path 'JPEGImages/'];
voc_edge_GT_src_path = [voc_path 'EdgGT_union/'];


edge_map_src_path = '/media/yuliangguo/NewVolume_1/EdgeBoxDollar/results/VOC2007/';
edge_map_dst_path = '../Data/SE_SEL_VOC2007/edges/';

dst_image_path = '../Data/VOC2007_img/';
dst_edge_GT_path = '../Data/VOC2007_GT/';

for i = 1:length(voc_seg_img_list)
    i
    name = voc_seg_img_list{i};
    copyfile([voc_image_path name '.jpg'], [dst_image_path name '.jpg'])
    copyfile([voc_edge_GT_src_path name '.mat'], [dst_edge_GT_path name '.mat'])
    copyfile([edge_map_src_path name '.edg'], [edge_map_dst_path name '.edg'])
    
end
