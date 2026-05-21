# rayextract trees refactor: segment/reconstruct split + multistem support

## Branch goal
Refactor `rayextract trees` to (1) cleanly separate segmentation from
reconstruction at the raylib level, and (2) treat (tree_id, stem_id) as
the reconstructable unit, not tree_id alone. Add a `--mask` input to the
segmentation step for carrying labels forward from a previous collection.

## Domain context (NON-NEGOTIABLE — read before designing anything)
- tree_id identifies one inventory-matched tree.
- stem_id identifies one reconstructable stem within that tree.
- A tree may have multiple stems. Stems may diverge at ground (e.g.
  mallee, lignotuberous Eucalyptus) or may split off the side of a main
  stem below girth height. In raycloudtools' current default behaviour,
  the latter would be merged into a single tree with branches; for our
  data, the stem separation has already been done upstream and must be
  preserved.
- A `reconstructable unit` is therefore a (tree_id, stem_id) pair, not
  a tree_id. Reconstruction must operate per-unit, then group outputs
  by tree_id for multistem-aware reporting.
- When stem_id is absent or all stems within a tree share one stem_id,
  behaviour must reduce exactly to the current rayextract trees output.

## Invariants
- The existing `rayextract trees` CLI must continue to work and produce
  byte-identical (or behaviourally identical, modulo intentional fixes)
  output when called without --mask and without per-stem input. We are
  adding capability, not breaking it.
- Segmentation outputs (tree_ids, stem_ids, trunk metadata) must be
  sufficient inputs for reconstruction. No hidden state crosses the
  boundary.
- The PLY loader's coordinate shift normalisation must be preserved
  through any in-memory path we add — large UTM coordinates must not
  reach float32 internals unshifted.

## --mask semantics
- `--mask` accepts either:
  (a) a directory of PLY files named `{tree_id}_{stem_id}.ply`
      (stem_id optional in filename if not present in data), OR
  (b) a single .las/.laz file with extra-bytes fields `tree_id` and
      optionally `stem_id`.
- Each unique (tree_id, stem_id) becomes a trunk seed for segmentation.
- Output tree_ids and stem_ids in the new cloud match those in the mask
  for matched stems. ID-collision policy for stems present in the new
  cloud but not the mask is configurable (default: new IDs starting
  beyond max mask ID).
- Mask cloud is assumed co-registered with input cloud. We will not
  attempt registration here.

## CLI shape (NON-NEGOTIABLE)

Two new subcommands are added under `rayextract`, mirroring the existing
`trunks → forest --trunks` composition pattern:

  rayextract segment cloud.ply [--mask <path>] [--ground cloud_mesh.ply]
    → cloud_segmented.ply  (labelled cloud, tree_id + stem_id fields)
    → cloud_seeds.txt      (per-stem trunk metadata)

  rayextract reconstruct cloud_segmented.ply cloud_mesh.ply
                         [--seeds cloud_seeds.txt]
    → cloud_trees.txt      (branch structures, treetools-compatible)
    → cloud_trees_mesh.ply (branch mesh)

  rayextract trees cloud.ply cloud_mesh.ply       # unchanged
    ≡ segment then reconstruct, in-memory, no intermediate files

Constraints:
- `rayextract trees` behaviour is preserved bit-for-bit in the no-mask,
  no-pre-segmented case. It is the back-compat entry point.
- The seeds file format extends `cloud_trunks.txt` by adding `tree_id`
  and `stem_id` columns. Existing tools reading `cloud_trunks.txt`
  (notably `rayextract forest --trunks`) must continue to work
  unchanged when handed our seeds file (they ignore extra columns).
- The library API exposes `ray::segment()`, `ray::reconstruct()`, and
  `ray::trees()` (the composition) as first-class functions. Nothing
  is reachable only via the CLI; nothing is reachable only via the
  library. CLI subcommands are thin `main()` wrappers.
- `rayextract reconstruct` accepts a labelled cloud with no seeds file
  and derives seeds by clustering low points per (tree_id, stem_id)
  and fitting trunks. A warning is printed in this mode.
