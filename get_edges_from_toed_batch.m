clear all; close all;

%> define source and output save directories
dataset_dir = "/gpfs/data/bkimia/cchien3/test_shock_code/shock-test-suite/imagenette/";
img_src_path = strcat(dataset_dir, "images");
edg_dst_path = strcat(dataset_dir, "edges");

if ~exist(edg_dst_path, "dir")
    error(strcat(edg_dst_path, " does not exist!"));
end

%> collect all JPEG images (including category subfolders)
img_files = dir(fullfile(img_src_path, "**", "*.JPEG"));

%> Detect third-order edges parameters
sigma = 2;
n = 1; % interpolate 2-fold
threshold = 2;
QaD_flag = 1;

for i = 1:length(img_files)
    fprintf("Processing %d/%d: %s\n", i, length(img_files), img_files(i).name);

    img_path = fullfile(img_files(i).folder, img_files(i).name);
    img_category = extractAfter(string(img_files(i).folder), "imagenette/images/");
    [~, img_stem, ~] = fileparts(img_files(i).name);
    export_path = fullfile(edg_dst_path, img_category, strcat(img_stem, "_to.edg"));

    %> read the image
    img = imread(img_path);
    opts.w = size(img, 2);
    opts.h = size(img, 1);

    %> Detect third-order edges
    [TO_edges, ~] = third_order_edge_detector(img, sigma, n, threshold, QaD_flag);

    %> Check that TO edges are within the image boundaries
    invalid_indices_x = find(TO_edges(:,1) < 2);
    TO_edges(invalid_indices_x,:) = [];
    invalid_indices_y = find(TO_edges(:,2) < 2);
    TO_edges(invalid_indices_y,:) = [];

    invalid_indices_x = find(TO_edges(:,1) > (opts.w-2));
    TO_edges(invalid_indices_x,:) = [];
    invalid_indices_y = find(TO_edges(:,2) > (opts.h-2));
    TO_edges(invalid_indices_y,:) = [];

    %> change to cxx coordinates
    TO_edges(:,1) = TO_edges(:,1)-1;
    TO_edges(:,2) = TO_edges(:,2)-1;
    save_edg(export_path, TO_edges, [opts.w opts.h]);
end

fprintf("Done. Wrote %d .edg files to %s\n", length(img_files), edg_dst_path);
