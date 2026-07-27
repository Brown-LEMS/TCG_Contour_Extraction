## Topological Contour Graph (TCG) Extraction

This repository implements the topological contour graph (TCG) extractor. Originally developed by Yuliang Guo; documented by Hongyi Fan and Chiang-Heng Chien. Used, edited, and updated by Chiang-Heng Chien.

## Prerequisite data
1. Edges detected on the image, _e.g._ SE / third-order edge detector. The edges are formatted in a `.edg` file.
2. Compute initial curve fragments from edges using `dborl_compute_curve_frags` to produce a `.cem` file.

## MATLAB usage

The main script is `main_TCG.m` which undergoes several steps from curve fragements to a topological contour graph:
  - loads the image, `.edg`, and `.cem`
  - breaks contours at corners
  - fills large gaps (DP)
  - breaks at T-junctions, prunes noise, merges geometrically, classifies junctions (BP)
  - breaks at corners again and prunes

## C++ port (`CPP/`)

The C++ code converts the MATLAB code in `main_TCG.m` (Step 3 and onward). It takes the original image, the edge file (.edg), and the contour fragement file (.cem) as inputs, and returns a new `.cem` file containing the new contour fragments of the topological contour graph.

### Dependencies

- CMake ≥ 3.14
- C++17 compiler
- OpenCV 4 and above built with `contib` modules

### Build
Follow the standard build and compile process to produce the executable `TCG`.
```bash
cd CPP
mkdir -p build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
cmake --build . -j
```
Once built is successful, provide specific `.edg`, `.cem`, and the original image file with optional threshold `ori_diff` (corner angle threshold in radians used to break long curves) and output file name `output.cem`. By default, the output file name is `<image_name>_tcg_cpp.cem` in the current working directory. 
```bash
./TCG <input.edg> <input.cem> <input.image> [ori_diff_th] [output.cem]
```
Using the files under `example_data`,
```bash
./build/TCG ./example_data/edges/n03425413_14351.edg \
            ./example_data/contours/n03425413_14351.cem \
            ./example_data/images/n03425413_14351.JPEG \
            ./outputs/n03425413_14351_tcg_cpp.cem
```
Console output reports fragment counts and timings for each stage.

### Visualize C++ results in MATLAB

The output `.cem` is readable by `load_contours` / `draw_contours`:
```matlab
addpath(genpath('util'));
img = imread('../example_data/images/cable.png');
[CEM_cpp, ~, ~] = load_contours('outputs/cable_tcg_cpp.cem');

figure; imshow(img, 'border', 'tight'); hold on;
draw_contours(CEM_cpp{2}, 0, 1);
title('C++ final contours');
```

### Notes / known differences

- "Corner break" (step 3) matches MATLAB closely on tested examples, _e.g._ fragment/corner counts.
- "Gap fill" (step 4) differ slightly from MATLAB because morphological skeleton (`bwmorph('skel')` used in MATLAB versus OpenCV Guo–Hall thinning), distance transform, and `imgradient` / Sobel are not bit-identical. Later stages inherit that discrepancy.
- Indices inside C++ are 0-based and `.cem` files keep the usual 0-based edge IDs. MATLAB loaders convert to 1-based when needed.
