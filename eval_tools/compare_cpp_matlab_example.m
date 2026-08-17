% Headless MATLAB vs C++ comparison on the *_to / *_to_dborl example inputs
% (same pairing as example_data/get_cem_from_tcg_matlab.m).

repo = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(repo, 'util')));

img_stem = 'n03425413_14351';
edg_stem = [img_stem '_to'];
img_path = fullfile(repo, 'example_data', [img_stem '.JPEG']);
edg_path = fullfile(repo, 'example_data', [edg_stem '.edg']);
cem_path = fullfile(repo, 'example_data', [edg_stem '_dborl.cem']);
ref_cemv = fullfile(repo, 'example_data', [edg_stem '_tcg_matlab.cemv']);

out_dir = fullfile(repo, 'outputs');
if ~exist(out_dir, 'dir'), mkdir(out_dir); end
matlab_cem = fullfile(out_dir, [edg_stem '_tcg_matlab.cem']);
matlab_cemv = fullfile(out_dir, [edg_stem '_tcg_matlab.cemv']);
cpp_cem = fullfile(out_dir, [edg_stem '_tcg_cpp.cem']);
cpp_cemv = fullfile(out_dir, [edg_stem '_tcg_cpp.cemv']);

params.vis = 0;
params.DP_gap_range = 15;
params.DP_angle_th = pi/4;
params.DP_contrast_th = 0.1;
params.shape_gap_range = 8;
params.shape_ori_range = pi/9;
params.BP_clen_th = 15;
params.BP_merge_angle_th = pi/9;
params.BP_nbr_num_edges = 20;
params.geom_merge_angle_th = pi/6;
params.nbr_num_edges = 20;
params.corner_angle_th = pi/6;
params.noise_prob_th = 0.05;
params.noise_len_th = 5;

img = imread(img_path);
[h, w, ~] = size(img);

[~, edgemap, ~] = load_edg(edg_path);
edgemap_soft0 = edgemap;
[CEM, edges, cfrags_idx] = load_contours(cem_path);

fprintf('[MATLAB] inputs: %s | %s\n', edg_path, cem_path);
fprintf('[MATLAB] initial contours = %d\n', numel(CEM{2}));

[new_cfrags, new_cfrags_idx, corner_pts] = contour_breaker_at_conner(CEM{2}, cfrags_idx, params, pi/18);
fprintf('[MATLAB Step 3] fragments = %d / corners = %d\n', numel(new_cfrags), size(corner_pts,1));

if size(img, 3) == 3
  img_gray = rgb2gray(img);
else
  img_gray = img(:, :, 1);
end
[edgemap_soft, thetamap] = imgradient(img_gray);
edgemap_soft = edgemap_soft/max(edgemap_soft(:));
thetamap = thetamap/180*pi;
thetamap = wrapToPi(-thetamap+pi/2);
t0 = tic;
[new_cfrags, new_cfrags_idx, edges] = contour_fill_gaps_DP(new_cfrags, new_cfrags_idx, h, w, edges, edgemap_soft, thetamap, params);
fprintf('[MATLAB Step 4] fragments = %d / edges = %d / time = %.3f ms\n', numel(new_cfrags), size(edges,1), toc(t0)*1000);

t0 = tic;
[new_cfrags, new_cfrags_idx, edges] = break_contours_at_T_junctions(new_cfrags, new_cfrags_idx, edges);
fprintf('[MATLAB T-break] fragments = %d\n', numel(new_cfrags));

[new_cfrags, new_cfrags_idx] = prune_noise_curves(new_cfrags, new_cfrags_idx, edgemap_soft0, params);
fprintf('[MATLAB Prune] fragments = %d\n', numel(new_cfrags));

[~, new_cfrags, new_cfrags_idx] = merge_cfrags_graphical_model_geom(new_cfrags, new_cfrags_idx, edges, params);
fprintf('[MATLAB Merge-geom] fragments = %d\n', numel(new_cfrags));

[~, new_cfrags, new_cfrags_idx, T_junctions] = classify_junction_type_wrt_graph_BP(new_cfrags, new_cfrags_idx, edges, params);
fprintf('[MATLAB Classify-BP] fragments = %d / T junctions = %d\n', numel(new_cfrags), size(T_junctions,1));

[new_cfrags, new_cfrags_idx, corner_pts] = contour_breaker_at_conner(new_cfrags, new_cfrags_idx, params);
fprintf('[MATLAB Corner2] fragments = %d / corners = %d\n', numel(new_cfrags), size(corner_pts,1));

[new_cfrags, new_cfrags_idx] = prune_noise_curves(new_cfrags, new_cfrags_idx, edgemap_soft0, params);
fprintf('[MATLAB Final] fragments = %d / elapsed_post = %.3f ms\n', numel(new_cfrags), toc(t0)*1000);

% Same writers as get_cem_from_tcg_matlab.m
write_cem_fixed_edge_id(matlab_cem, new_cfrags, edges, h, w, new_cfrags_idx);
det_save_cemv(matlab_cemv, new_cfrags);
fprintf('[MATLAB Write] %s\n', matlab_cem);
fprintf('[MATLAB Write] %s\n', matlab_cemv);

compare_cemv_pair(ref_cemv, matlab_cemv, 'your saved cemv', 'regenerated MATLAB cemv');
if exist(cpp_cemv, 'file')
  compare_cemv_pair(ref_cemv, cpp_cemv, 'your saved cemv', 'C++ cemv');
  compare_cemv_pair(matlab_cemv, cpp_cemv, 'regenerated MATLAB cemv', 'C++ cemv');
else
  fprintf('[Compare] missing C++ cemv: %s\n', cpp_cemv);
end

function compare_cemv_pair(path_a, path_b, name_a, name_b)
  ca = load_cemv_contours(path_a);
  cb = load_cemv_contours(path_b);
  fprintf('\n=== %s vs %s ===\n', name_a, name_b);
  fprintf('%s: %d contours, %d pts\n', name_a, numel(ca), sum(cellfun(@(c)size(c,1), ca)));
  fprintf('%s: %d contours, %d pts\n', name_b, numel(cb), sum(cellfun(@(c)size(c,1), cb)));

  eps_a = zeros(0,2); eps_b = zeros(0,2);
  for i = 1:numel(ca)
    if isempty(ca{i}), continue; end
    eps_a(end+1,:) = ca{i}(1,1:2); %#ok<AGROW>
    eps_a(end+1,:) = ca{i}(end,1:2);
  end
  for i = 1:numel(cb)
    if isempty(cb{i}), continue; end
    eps_b(end+1,:) = cb{i}(1,1:2); %#ok<AGROW>
    eps_b(end+1,:) = cb{i}(end,1:2);
  end
  key_a = unique(round(eps_a*100), 'rows');
  key_b = unique(round(eps_b*100), 'rows');
  matched = 0;
  used = false(size(key_b,1),1);
  for i = 1:size(key_a,1)
    d = sum((double(key_b) - double(key_a(i,:))).^2, 2);
    [mind, j] = min(d);
    if ~isempty(mind) && mind <= 100^2 && ~used(j)
      used(j) = true;
      matched = matched + 1;
    end
  end
  fprintf('Unique endpoints: %d / %d\n', size(key_a,1), size(key_b,1));
  fprintf('Endpoint matches within 1px: %d (%.1f%% of %s)\n', matched, ...
          100*matched/max(1,size(key_a,1)), name_a);
end

function contours = load_cemv_contours(path)
  txt = fileread(path);
  blocks = regexp(txt, '\[BEGIN CONTOUR\]\nEDGE_COUNT=(\d+)\n(.*?)\[END CONTOUR\]', 'tokens');
  contours = cell(numel(blocks), 1);
  for i = 1:numel(blocks)
    n = str2double(blocks{i}{1});
    body = blocks{i}{2};
    lines = regexp(body, '\n', 'split');
    pts = zeros(0,2);
    for k = 1:numel(lines)
      tok = regexp(lines{k}, '\[0,\s*0\]\s+0\.0+\s+0\.0+\s+\[([-\d\.eE\+]+),\s*([-\d\.eE\+]+)\]', 'tokens', 'once');
      if isempty(tok), continue; end
      pts(end+1,:) = [str2double(tok{1}), str2double(tok{2})]; %#ok<AGROW>
    end
    if size(pts,1) ~= n
      error('parse mismatch in %s contour %d: expected %d got %d', path, i, n, size(pts,1));
    end
    contours{i} = pts;
  end
end
