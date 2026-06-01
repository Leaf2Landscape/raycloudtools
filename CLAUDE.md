# rayextract refactor: minimal-diff segment/reconstruct split

## Current state

`rayextract segment` and `rayextract reconstruct` are **functionally correct and complete**.
Do not change behavior. This refactor is about code structure only.

### What each command does

```
rayextract segment <cloud.ply> [--ground mesh.ply] [--mask <path>]
  -> <prefix>_segmented.las   (per-point tree_id, stem_id)
  -> <prefix>_segmented_seeds.txt

rayextract reconstruct <cloud_segmented.las> <ground_mesh.ply>
  -> <prefix>_trees.txt
  -> <prefix>_trees_mesh.ply

rayextract trees <cloud.ply> <ground_mesh.ply>   (unchanged legacy behavior)
```

### Why reconstruct re-runs segmentation internally

`reconstruct` accepts a pre-labeled cloud (tree_id per point) and must produce branch geometry
(cylinders, radii, taper). The branch reconstruction algorithm (`reconstructBranches`) needs a
Dijkstra-derived parent-child graph over all points. There is no way to produce this graph from
point labels alone without re-running the path-finding.

The current solution: the `PreLabeledTag` constructor in `Trees` runs **per-tree Dijkstra**, one
pass per unique `(tree_id, stem_id)` group, using the pre-loaded labels to partition points.
This re-runs path-finding while **honoring the input allocation** — points are never re-assigned
to a different tree.

## Refactor goal: minimize diff to main

The implementation works. The problem is it has too much divergence from the upstream `main`
branch, particularly inside `raytrees.h/.cpp` (original CSIRO authorship). The goal is to
restructure so that the delta to main is as small and as reviewable as possible.

### Target diff boundaries

| File | Target |
|---|---|
| `raytrees.h` / `raytrees.cpp` | Minimize. Every added line needs justification. |
| `raytreeslib.h` / `raytreeslib.cpp` | New files — this is where new logic should live. |
| `raysegmentresult.h` | New file — data types only, no logic. |
| `raymaskloader.h` / `raymaskloader.cpp` | New files — mask I/O only. |
| `rayextract.cpp` | Additive-only: new subcommand blocks, no changes to existing blocks. |
| `raysegment.h` / `raysegment.cpp` | Do not touch unless a bug fix is required. |
| `raycloud.h` / `raycloud.cpp` | Minimize. `tree_ids` / `stem_ids` fields already present from LAS PR. |

### Specific changes to raytrees.h/.cpp to justify or remove

Currently added to raytrees.h (from main):
- `PreLabeledTag` struct — needed; keep
- `Trees(cloud, offset, mesh, params, verbose, PreLabeledTag)` — needed; keep
- `buildSeedList()` public method — evaluate whether it can live in raytreeslib.cpp instead
- `buildLabelIdMap()` public method — needed by reconstruct; keep
- `reconstructBranches()` private method — needed for sharing between constructors; keep
- `sec_labels_` private field — needed by PreLabeledTag constructor; keep
- optional `id_map` param on `save()` — needed; keep

If `buildSeedList()` only needs data already accessible via `save()` output, consider deriving
seeds in `raytreeslib.cpp` from the saved text file or cloud rather than exposing a new method.
Only keep it on `Trees` if no cleaner alternative exists.

## Non-negotiable constraints

- `rayextract trees` default behavior must remain bit-identical to main (no mask, no new flags).
- `rayextract segment` output format (`_segmented.las`, `_segmented_seeds.txt`) must not change.
- `rayextract reconstruct` output format (`_trees.txt`, `_trees_mesh.ply`) must not change.
- Do not modify any existing `rayextract trees` parameter handling in `rayextract.cpp`.
- Do not change `raysegment.cpp` / `raysegment.h` core Dijkstra logic.
- Preserve PLY coordinate-shift normalization across all paths.

## Mask feature status

Mask support (`--mask`) is implemented and working at a first-draft level. It is in-scope for
this refactor but should not be extended. The sidecar CSV output and `loadMask()` are correct.
Do not add registration, do not add new mask formats.

## Testing requirements

- Confirm `rayextract trees` output is unchanged from main (no-mask case).
- Confirm `rayextract segment` followed by `rayextract reconstruct` produces output equivalent
  to `rayextract trees` on the same input (no-mask case).
- Any refactor that changes raytrees.cpp must be verified against both paths.
- New tests for mask parsing and sidecar CSV are deferred until the minimal-diff refactor lands.

## Do not do

- Do not refactor or rename anything inside `getRootsAndSegment` or `connectPointsShortestPath`.
- Do not add `--seeds` as a separate CLI argument (seeds are loaded from the segmented cloud).
- Do not add stem_id to the PLY color-encoding path.
- Do not add sidecar CSV unless `--mask` is active.
- Do not introduce new abstractions or helper classes beyond what is needed.
