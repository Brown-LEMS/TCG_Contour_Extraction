clear all; close all;

%> define source and output save directories
dataset_dir = "/gpfs/data/bkimia/cchien3/test_shock_code/shock-test-suite/imagenette/";
edg_src_path = strcat(dataset_dir, "edges");
cem_dst_path = strcat(dataset_dir, "contours");

if ~exist(cem_dst_path, "dir")
    error(strcat(cem_dst_path, " does not exist!"));
end

%> get the current code absolute path
current_script_path = fileparts(mfilename('fullpath'));

%> collect all edge files, specified by the edge type (including category subfolders)
edge_files = dir(fullfile(edg_src_path, "**", "*_to.edg"));

for i = 1:length(edge_files)
    fprintf("Processing %d/%d: %s\n", i, length(edge_files), edge_files(i).name);

    edge_path = fullfile(edge_files(i).folder, edge_files(i).name);
    edge_category = extractAfter(string(edge_files(i).folder), "imagenette/edges/");
    [~, edg_stem, ~] = fileparts(edge_files(i).name);
    export_path = fullfile(cem_dst_path, edge_category, strcat(edg_stem, "_dborl.cem"));
    
    exec_file = [current_script_path '/dborl_compute_curve_frags '];
    cmd = ['!' exec_file char(edge_path) ' ' char(export_path)];
    eval(cmd);
end

fprintf("Done. Wrote %d .cem files to %s\n", length(edge_files), cem_dst_path);








    


