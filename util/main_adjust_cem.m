clear all; close all;

edg_src_path = '../Data/gPb_KGS_CFGD/';

cem_files = dir([edg_src_path '*.cem']);

adjust = -1;

for i = 1:length(cem_files)
    i
    cem_file_in = [edg_src_path cem_files(i).name];
    cem_file_out = [edg_src_path cem_files(i).name];
    
    adjust_cood_cem (cem_file_in, cem_file_out, adjust);
end