clear all; close all;

edg_src_path = '../Data/gPb_SEL_VOC2007/edges/';

edg_files = dir([edg_src_path '*.edg']);

adjust = -1;

for i = 1:length(edg_files)
    i
    edg_file_in = [edg_src_path edg_files(i).name];
    edg_file_out = [edg_src_path edg_files(i).name];
    
    adjust_cood_edg (edg_file_in, edg_file_out, adjust);
end