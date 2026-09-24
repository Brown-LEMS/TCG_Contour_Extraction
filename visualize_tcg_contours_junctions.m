% visualize_tcg_contours_junctions.m
% Independent script: overlay TCG contours (.cem) and junctions (.jct) on an image.
%
% Coordinate note: .cem / .jct use 0-based coords, so plot with +1 in matlab

clear; close all;

script_dir = fileparts(mfilename('fullpath'));
img_path = fullfile(script_dir, 'example_data', 'n03425413_14351.JPEG');
cem_path = fullfile(script_dir, 'example_data', 'n03425413_14351_to_tcg_cpp.cem');
jct_path = fullfile(script_dir, 'example_data', 'n03425413_14351_to_tcg_cpp.jct');

addpath(genpath(fullfile(script_dir, 'util')));

if ~isfile(img_path), error('Image not found: %s', img_path); end
if ~isfile(cem_path), error('CEM not found: %s', cem_path); end
if ~isfile(jct_path)
    warning(['JCT not found: %s\n' ...
             'Re-run the C++ TCG binary (rebuild after junction support) to produce .jct.'], ...
            jct_path);
end

img = imread(img_path);
[CEM, ~, ~] = load_contours(cem_path);
cfrags = CEM{2};

[T_jcts, Y_jcts] = load_tcg_jct(jct_path);

figure('Name', 'TCG contours + junctions');
imshow(img, 'border', 'tight'); hold on;
draw_contours(cfrags, 0, 1);

h_leg = gobjects(0);
leg_str = {};
if ~isempty(T_jcts)
    hT = plot(T_jcts(:,1) + 1, T_jcts(:,2) + 1, 'go', ...
              'MarkerSize', 8, 'LineWidth', 1.5, 'MarkerFaceColor', 'w');
    h_leg(end+1) = hT;
    leg_str{end+1} = sprintf('T junctions (%d)', size(T_jcts, 1));
end
if ~isempty(Y_jcts)
    hY = plot(Y_jcts(:,1) + 1, Y_jcts(:,2) + 1, 'rs', ...
              'MarkerSize', 8, 'LineWidth', 1.5, 'MarkerFaceColor', 'w');
    h_leg(end+1) = hY;
    leg_str{end+1} = sprintf('Y junctions (%d)', size(Y_jcts, 1));
end
if ~isempty(h_leg)
    lgd = legend(h_leg, leg_str, 'Location', 'northwest');
    lgd.Color = [0.9 0.9 0.9];
end
hold off;
title(sprintf('Contours=%d | T=%d | Y=%d', numel(cfrags), size(T_jcts,1), size(Y_jcts,1)));

%% ---- local loader for the C++ .jct format ----
function [T_jcts, Y_jcts] = load_tcg_jct(path)
%LOAD_TCG_JCT Read [T_junctions] / [Y_junctions] blocks written by CPP write_jct.
% Each row: x y dir conf  (0-based x,y). Returns Nx4 (possibly empty).
    T_jcts = zeros(0, 4);
    Y_jcts = zeros(0, 4);
    if ~isfile(path), return; end

    fid = fopen(path, 'r');
    if fid < 0, error('Cannot open %s', path); end
    cleaner = onCleanup(@() fclose(fid));

    section = '';
    while true
        line = fgetl(fid);
        if ~ischar(line), break; end
        line = strtrim(line);
        if isempty(line) || startsWith(line, '#'), continue; end
        if startsWith(line, 'size='), continue; end

        if strcmp(line, '[T_junctions]')
            section = 'T'; continue;
        elseif strcmp(line, '[Y_junctions]')
            section = 'Y'; continue;
        end

        if startsWith(line, 'count=')
            continue;
        end

        vals = sscanf(line, '%f');
        if numel(vals) < 2, continue; end
        row = [vals(1), vals(2), 0, 0];
        if numel(vals) >= 3, row(3) = vals(3); end
        if numel(vals) >= 4, row(4) = vals(4); end

        if strcmp(section, 'T')
            T_jcts(end+1, :) = row;
        elseif strcmp(section, 'Y')
            Y_jcts(end+1, :) = row;
        end
    end
end
