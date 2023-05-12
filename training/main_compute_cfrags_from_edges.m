clear all; close all;
addpath (genpath('../util/'));

edge_src_path = '../Data/gPb_SEL_CFGD/edges/';
cfrags_dst_path = '../Data/gPb_SEL_CFGD/cfrags/';
mkdir(cfrags_dst_path);

edge_files = dir([edge_src_path '*.edg']);

for i = 1:length(edge_files)
   
    input = [edge_src_path edge_files(i).name];
    output = [cfrags_dst_path edge_files(i).name(1:end-4) '.cem'];
    
    cmd = ['!../util/dborl_compute_curve_frags ' input ' ' output];
    eval(cmd);
    
end