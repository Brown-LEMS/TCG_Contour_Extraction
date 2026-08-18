# Curve fragments comparison

Given a set of curve fragments, the goal is to compare another set of curve fragments and evaluate the consistency/similarity of the two sets. The comparison is asymmetric, meaning that the first set of curves are treated as the reference (historically it is referred to as the "ground-truth (GT)"), while the second is treated as the test map (historically it is referred to as the "computed (CP)"). Precision and recall rates are used as the evaluation metrics, and both are defined with respect to the GT-CP pair assignment.

## Core Files

| File | Role |
|---|---|
| `util/cfgd_eval/Compare_Curve_Fragment_Maps.m` | Compare two maps given file paths pointing to the GT and CP cem files |
| `util/cfgd_eval/Compare_Curve_Fragment_Maps_2.m` | Compare two maps given loaded CEM structs |
| `util/cfgd_eval/Compare_Curve_Fragment_Maps_v2_new_score.m` | To be determined |
| `Eval/Main_eval_strict_CFGD_new_prune.m` | To be determined |
| `Eval/Main_eval_loose_CFGD_new_prune.m` | To be determined |

## CEM Contour representation

A `.cem` file can be loaded with `load_contours` function in `util/io/load_contours.m`. For example,
```matlab
[CEM, edges, cfrags_idx] = load_contours('image.cem');
```
where `CEM` is a cell array:
- `CEM{1}` — image size `[width, height]`
- `CEM{2}` — cell array of curve fragments. Fragment `k` is an `N x 5` matrix of edges `(x, y, orientation, confidence, d2f)`
- `CEM{3}` — optional per-contour properties

## Evaluation Metrics

All variants return four length (or pixel-count) quantities:
| Output | Meaning |
|---|---|
| `TP_GT_L` | Matched portion of the reference map (true positives (TP) of the "ground-truth" curve fragments) |
| `TP_CP_L` | Matched portion of the test map (true positives (TP) of the "computed" curve fragements) |
| `GT_L` | Total length / pixel count of the reference map |
| `CP_L` | Total length / pixel count of the test map |

From those we can calculate the precision and recall rates, _i.e._,
```
Recall    = TP_GT_L / GT_L     % fraction of the reference that was recovered
Precision = TP_CP_L / CP_L     % fraction of the test map that matches the reference
F         = 2 * P * R / (P + R)
```

Additionally, there are optional extra outputs
- `match_GT_frags` / `match_CP_frags`: fragments that were accepted as matches
- `miss_GT_frags`: reference fragments with no match
- `extra_CP_frags`: test fragments with no match
- `match_select_frags`: for each matched group, the side with fewer fragments (or, on a tie, the longer side)

## Usage

```matlab
[TP_GT_L, TP_CP_L, GT_L, CP_L, match_GT_frags, match_CP_frags, miss_GT_frags, extra_CP_frags, match_select_frags] = ...
    Compare_Curve_Fragment_Maps(GT_contour_path, CP_contour_path, prune_len_thresh, ...
                                display_image, bgroup_only, cost_thresh, edit_thresh, local_dist);
```
`Compare_Curve_Fragment_Maps_2` function does basically the same thing, except that it takes the `CEM` cell array loaded from the paths.

### Inputs

| Argument | Default | Description |
|---|---|---|
| `GT_contour_path` | required | Path to the reference `.cem` or `.mat` |
| `CP_contour_path` | required | Path to the test `.cem` or `.mat` |
| `prune_len_thresh` | `0` | Test fragments shorter than this are excluded from matching. `0` means no prune. |
| `display_image` | `0` | If `> 0`, draw three figures: refined fragments, matched groups, and misses/extras. |
| `bgroup_only` | `0` | If `> 0`, skip edit-distance refinement and accept whole spatial groups as matches. |
| `cost_thresh` | `2` | Cost above which two fragments are treated as unrelated (sentinel `1000` in the cost matrix). Also used while grouping. |
| `edit_thresh` | `1.2` | Maximum transform / edit cost for an accepted fragment combination. |
| `local_dist` | `2` | Allowed localization error (pixels) when chaining fragments into combinations. |

Note that _(i)_ `cost_thresh`, `edit_thresh`, and `local_dist` are all initialized if `nargin < 5`, _(ii)_ the default thresholds in the `Compare_Curve_Fragment_Maps_2` function is looser than `Compare_Curve_Fragment_Maps`.

### Minimal example
Refer to `Eval/minimal_eval_example.m` for evaluating the CPP and MATLAB version using the example data.