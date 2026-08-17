clear all; close all;

%> Load two contour files with CEM file extension
current_script_path = fileparts(mfilename('fullpath'));
cem_file_1 = fullfile(current_script_path, "outputs", "n03425413_14351_tcg_matlab.cem");
cem_file_2 = fullfile(current_script_path, "outputs", "n03425413_14351_tcg_cpp.cem");

img_path = fullfile(current_script_path, "../example_data/images", "n03425413_14351.JPEG");
img = imread(img_path);
if length(size(img)) > 2
    img = rgb2gray(img);
end

[CEM_1, edges_1, cfrags_idx_1] = load_contours(cem_file_1);
[CEM_2, edges_2, cfrags_idx_2] = load_contours(cem_file_2);

contour_min_length = 0;
rand_color = 1;

fig = figure(1);
imshow(img, 'border', 'tight'); 
hold on;
draw_contours(CEM_1{2}, contour_min_length, rand_color);
fig.Name = "Contour fragments from MATLAB overlayed on the image";
hold off;
pause(0.5);

fig = figure(2);
imshow(img, 'border', 'tight'); 
hold on;
draw_contours(CEM_2{2}, contour_min_length, rand_color);
fig.Name = "Contour fragments from CPP overlayed on the image";
