# Data model — segment/reconstruct split and multistem support

Companion to `docs/upstream_map.md`. Specifies every type, file format, and
policy decision that crosses a component boundary. No code is written here.

---

## 1. In-memory boundary between `ray::segment()` and `ray::reconstruct()`

### 1.1 Per-point label vectors

Two parallel `int32_t` vectors, one entry per **bounded ray** (i.e. per point
satisfying `cloud.rayBounded(i)`). This is the same indexing convention as
`points_[]` inside the current `Trees` constructor and as `cloud.tree_ids[]`
after `segmentCloud` (which also uses `cloud.ends`-parallel storage with -1
sentinels for unbounded rays).

```
std::vector<int32_t>  tree_ids   // parallel to bounded-ray subset of cloud.ends
std::vector<int32_t>  stem_ids   // parallel to bounded-ray subset of cloud.ends
```

When stored in `Cloud` for serialisation (sections 2, 8), both vectors are
**`cloud.ends`-parallel** (all rays), with sentinel **-1** for unbounded rays
and for points that were not assigned to any tree or stem. This matches the
existing `cloud.tree_ids` convention (`segmentCloud`, raytrees.cpp:1127).

**Rationale for two separate vectors**: `tree_ids` already exists in `Cloud`
(`raycloud.h:41`); adding a parallel `stem_ids` at `raycloud.h` is the minimal
extension (see upstream_map §H). A packed `int64_t` combining both would save
one vector allocation but would break the existing `cloud.tree_ids` interface
and the LAS extra-bytes layout.

| Vector | Type | Units | Valid range | Sentinel |
|--------|------|-------|-------------|---------|
| `tree_ids` | `int32_t` | — | 0 … INT32\_MAX | -1 (unassigned) |
| `stem_ids` | `int32_t` | — | 0 … INT32\_MAX | -1 (unassigned, or point not in any stem) |

`stem_id = 0` is the canonical single-stem value. When stem_ids are absent from
all inputs and the single-stem code path runs, every assigned point carries
`stem_id = 0`. This is the backwards-compatible default (see section 7).

---

### 1.2 `StemSeed` and `SeedList`

`SeedList` is the per-stem metadata that `ray::reconstruct()` needs from
`ray::segment()`. It is also serialised to `cloud_seeds.txt` (section 3).

```
struct StemSeed {
    int32_t          tree_id;       // inventory tree ID; matches tree_ids[] above
    int32_t          stem_id;       // stem within tree; matches stem_ids[] above
    Eigen::Vector3d  base;          // stem base position, local coords
    double           radius;        // trunk radius at girth-height window
    double           tree_height;   // height of tallest stem point above base
};

using SeedList = std::vector<StemSeed>;
```

Field documentation:

| Field | C++ type | Units | Valid range | Sentinel / default |
|-------|----------|-------|-------------|-------------------|
| `tree_id` | `int32_t` | — | 0 … INT32\_MAX | no sentinel; every StemSeed has a valid tree_id |
| `stem_id` | `int32_t` | — | 0 … INT32\_MAX | 0 = single-stem default |
| `base` | `Eigen::Vector3d` | metres | any finite float64 | — |
| `radius` | `double` | metres | > 0 | no sentinel; a radius of 0 is invalid and must be rejected |
| `tree_height` | `double` | metres | > 0 | no sentinel; a height of 0 is invalid |

**Coordinate frame for `base`**: local coordinates after `cloud.removeStartPos()`
has been applied (same frame as `points_[].pos` inside the `Trees` constructor).
When written to `cloud_seeds.txt` the `offset` returned by `removeStartPos()` is
re-added, matching the convention used by the existing trees.txt writer
(`raytrees.cpp:1241`). `ray::reconstruct()` receives a cloud that has already
had `removeStartPos()` applied, and the seeds file is read back with the same
offset subtracted.

**Why `tree_height` is in SeedList**: `girth_height_ratio × tree_height` sets the
vertical window for trunk-radius estimation (upstream_map §E, `raytrees.cpp:81`).
Storing `tree_height` in the seed lets `ray::reconstruct()` measure at exactly
the same window height that `ray::segment()` used, rather than re-deriving it
from the (now per-stem) point set and potentially getting a slightly different
value if the tallest point is at a different path distance.

**What is NOT in SeedList**: `distance_to_end`, `parent` links, `children[]`. These
are re-derived by `ray::reconstruct()` via a local Dijkstra seeded from
`StemSeed::base` within each `(tree_id, stem_id)` point set (upstream_map §D).
No intermediate path state crosses the segment/reconstruct boundary.

---

### 1.3 Full boundary type

```
struct SegmentResult {
    SeedList                 seeds;       // one entry per (tree_id, stem_id)
    std::vector<int32_t>     tree_ids;    // bounded-ray-parallel
    std::vector<int32_t>     stem_ids;    // bounded-ray-parallel
};
```

`ray::reconstruct()` signature (library level; CLI wrappers are thin):

```
ForestStructure reconstruct(
    Cloud&                  cloud,        // labelled; cloud.tree_ids and
                                          // cloud.stem_ids populated
    const Eigen::Vector3d&  offset,
    const Mesh&             mesh,
    const SeedList&         seeds,
    const TreesParams&      params,
    bool                    verbose
);
```

`ray::segment()` returns a `SegmentResult`. `ray::trees()` calls `segment()` then
passes its output directly to `reconstruct()` without serialisation, preserving
bit-identical floating-point state for the `rayextract trees` backwards-compat path.

---

## 2. `cloud_segmented` file format

### 2.1 PLY — changes required

The existing PLY writer (`rayply.cpp:165–188`) has a fixed 9-property binary
layout (`RayPlyEntry`, `rayply.h:15–21`). Neither `tree_id` nor `stem_id` is
stored today, and there is no dynamic property mechanism
(upstream_map §H, `rayply.cpp:37–80`).

**Conceptual PLY property names** (for when PLY support is added):

```
property int tree_id
property int stem_id
```

Both `int` (32-bit signed, PLY type `int`), appended after the existing nine
properties so old readers that parse the header to skip unknown properties will
ignore them gracefully.

**Minimum writer extension for PLY**:

| Component | Change |
|-----------|--------|
| `rayply.h:15–21` | Add `int32_t tree_id, stem_id` fields to `RayPlyEntry` |
| `rayply.cpp:37–80` | Append `property int tree_id` and `property int stem_id` to the PLY ASCII header, conditionally on the caller supplying non-empty vectors |
| `rayply.cpp:165–188` | Write the two `int32_t` values per point in `writePlyRayCloud()`; add `tree_ids` and `stem_ids` parameters |
| `rayply.cpp:330–770` | Detect `"property int tree_id"` and `"property int stem_id"` by name in the header parse loop; record their byte offsets; extract them into output vectors |
| `raycloud.cpp:40–41` | Thread `cloud.tree_ids` and `cloud.stem_ids` into `writePlyRayCloud()` |

This is a non-trivial change affecting four functions in two files.

### 2.2 LAS/LAZ — minimal extension (recommended for Level 1)

`tree_id` is already stored as int32 at extra-bytes byte offset 12
(upstream_map §H, `raylaz.cpp:673`). Adding `stem_id` follows the same pattern
at the next available slot:

**Extra-bytes byte layout with stem_id added**:

| Bytes | Field | Type | Notes |
|-------|-------|------|-------|
| 0–3 | sx | float32 | ray origin X − end X |
| 4–7 | sy | float32 | ray origin Y − end Y |
| 8–11 | sz | float32 | ray origin Z − end Z |
| 12–15 | tree_id | int32 | existing; -1 = unassigned |
| 16–19 | stem_id | int32 | **new**; -1 = unassigned |
| 20 | alpha | uint8 | moves from byte 16 to byte 20 |
| 21… | original sensor extra-bytes | — | preserved unchanged |

**Minimum writer extension for LAS** (~20 lines total, upstream_map §H):

| File | Change |
|------|--------|
| `raycloud.h:41` | Add `std::vector<int32_t> stem_ids;` after `tree_ids` |
| `raylaz.cpp:673` | Register `"stem_id"` as int32 EXTRA_BYTES attribute after `"tree_id"` |
| `raylaz.cpp:850–854` | Write `stem_ids[i]` at `point_->extra_bytes + 16`; move alpha write to `+20` |
| `raylaz.cpp:122–123` | Add `has_stem_id_attr` flag to VLR scan; detect `"stem_id"` by name |
| `raylaz.cpp:183–188` | Read back `int32_t` from byte 16 into `stem_ids_out` |
| `raylaz.cpp:268` | Extend alpha-offset logic: `has_stem_id_attr ? 20u : (has_tree_id_attr ? 16u : 12u)` |
| `raylaz.cpp:95` | Append `"stem_id"` to the known-names list |
| `raycloud.cpp:38–39` | Thread `cloud.stem_ids` through `writeLasRayCloud()` |

**Level 1 policy**: `rayextract segment` always writes `.las`. The output
filename matches the input extension if the input is `.las/.laz`; otherwise
it is forced to `.las` with a warning. PLY segmented output is deferred.

### 2.3 Field values

| Field | Value when assigned | Value when unassigned |
|-------|--------------------|-----------------------|
| `tree_id` | 0 … INT32\_MAX | -1 |
| `stem_id` | 0 … INT32\_MAX; 0 = single-stem default | -1 (or 0 if stem detection was not run) |
| RGB colour | `convertIntToColour(tree_id)` — same as today | (0, 0, 0) black — same as today |

The existing colour encoding is preserved for visual compatibility. It encodes
`tree_id` (not stem_id), matching current `segmentCloud` behaviour.

---

## 3. `cloud_seeds.txt` schema

### 3.1 Column order

```
# tree seeds file:
x,y,z,radius,tree_id,stem_id,tree_height
<x>, <y>, <z>, <radius>, <tree_id>, <stem_id>, <tree_height>
```

Columns in order:

| Position | Name | Type | Units | Notes |
|----------|------|------|-------|-------|
| 1 | x | float64 | metres | Trunk base easting, world coords (offset re-added) |
| 2 | y | float64 | metres | Trunk base northing, world coords |
| 3 | z | float64 | metres | Trunk base height, world coords |
| 4 | radius | float64 | metres | Trunk radius at girth-height measurement window; > 0 |
| 5 | tree_id | int32 (written as float64) | — | Inventory tree ID; 0-based |
| 6 | stem_id | int32 (written as float64) | — | Stem within tree; 0 = single-stem default |
| 7 | tree_height | float64 | metres | `StemSeed::tree_height`; drives girth window in reconstruct |

- **Comment line**: `# tree seeds file:` (mirrors `Trunks::save` at `raytrunks.cpp:686`)
- **Delimiter**: comma-space `", "` (matches `Trunks::save` at `raytrunks.cpp:695`)
- **Header**: no spaces between column name tokens (`x,y,z,radius,tree_id,stem_id,tree_height`)
- **One data line per stem**
- `tree_id` and `stem_id` are written as integers but the `ForestStructure::load`
  parser reads them as `double` (stored in `segment.attributes[0]`, `[1]`, `[2]`);
  integer round-trip through `double` is exact for `int32_t` values up to 2^53.
- `tree_height` is stored as float64; it is required by `ray::reconstruct()` so
  the girth-height measurement window (`girth_height = girth_height_ratio ×
  tree_height`, `raytrees.cpp:81`) is reproducible without re-deriving height
  from the point cloud.

### 3.2 Compatibility with `rayextract forest --trunks`

Per upstream_map §G, `rayextract forest --trunks` at `rayextract.cpp:323` calls
`ForestStructure::load()`, not `Trunks::load()`. `ForestStructure::load` stores
post-mandatory columns in `segment.attributes[]` (`rayforeststructure.cpp:250`)
and the forest extractor uses only `segment.tip` and `segment.radius`
(`rayextract.cpp:329`). The seeds file is therefore a strict superset of
`cloud_trunks.txt` for this use case. No change to the forest reader is needed.

The strict `Trunks::load()` (`raytrunks.cpp:717–739`) counts commas and rejects
any line that does not have exactly 3; it would break on a seeds file. However,
`Trunks::load()` is not called by `rayextract forest --trunks` and is not in the
pipeline described by CLAUDE.md. No change to `Trunks::load()` is needed.

### 3.3 `ForestStructure::load` parse walk-through for seeds format

The critical detail: both `commas_per_segment` (from the header) and
`num_commas` (from each data line) are computed with the formula
`1 + std::count(str.begin(), str.end(), ',')` (`rayforeststructure.cpp:155,
198`). Both therefore count **fields**, not separators. They are consistent.

Given header `x,y,z,radius,tree_id,stem_id,tree_height` (6 commas):
- `mandatory_text = "x,y,z,radius"` found at position 0 → `commas_per_tree = 0`
  (`rayforeststructure.cpp:120–151`)
- `lines[1]` = `"x,y,z,radius,tree_id,stem_id,tree_height"` → 6 commas →
  `commas_per_segment = 1 + 6 = 7` (`rayforeststructure.cpp:155`,
  computed **before** the `substr` trim)
- After `substr(mandatory_text.length())`: `lines[1]` = `",tree_id,stem_id,tree_height"`;
  attribute loop runs for `i = 4, 5, 6`: `attributes[1]` = `["tree_id", "stem_id", "tree_height"]`
- `has_parent_id = false`

For data line `"394123.1, 6012044.2, 18.5, 0.142, 7, 0, 22.3"` (6 commas):
- `num_commas = 1 + 6 = 7` (`rayforeststructure.cpp:198`)
- Validation: `(7 − 0) % 7 = 0` → **passes** ✓

`segment.attributes[0] = 7.0` (tree_id), `segment.attributes[1] = 0.0` (stem_id),
`segment.attributes[2] = 22.3` (tree_height). `rayextract forest --trunks` reads
only `segment.tip` and `segment.radius` (`rayextract.cpp:329`) and discards all
attributes. ✓

The seeds file format at §3.1 is fully compatible with `ForestStructure::load`
and with `rayextract forest --trunks` without any changes to either reader.

---

## 4. `cloud_trees.txt` schema extension

### 4.1 Choice: option (b)

`tree_id` and `stem_id` are added as **per-line (per-tree) attributes**, placed
before `x,y,z,radius` in the header. Each data line represents one
reconstructable unit — one `(tree_id, stem_id)` pair — and begins with those two
values.

**Extended header**:
```
# Tree file. tree_id,stem_id per stem; x,y,z,radius and optional per-segment attributes follow.
tree_id,stem_id, x,y,z,radius,parent_id,section_id[,weight,len,accuracy,junction_weight]
```

**Extended data line**:
```
<tree_id>, <stem_id>, <root_x>,<root_y>,<root_z>,<root_r>,-1,0, <c1_x>,...
```

### 4.2 Rationale over option (a)

| Criterion | Option (a): stem_id per-segment after section_id | Option (b): tree_id,stem_id per-line (chosen) |
|-----------|--------------------------------------------------|----------------------------------------------|
| Grouping stems to trees | Must scan all lines and aggregate by segment.attributes[-2] | Immediate: read first two fields of each line |
| Parser support | `segment.attributes[1]` (extra column is double) | `tree.treeAttributes()[0]` and `[1]` (per-tree doubles) |
| Single-stem compat | Writes stem_id=0 per-segment; parsers ignoring attributes unaffected | Writes tree_id, 0 per-line; parsers ignoring treeAttributes() unaffected |
| treetools compatibility | treetools can ignore `segment.attributes` | treetools can ignore `treeAttributes()` |
| tree_id in file | Would require adding a second per-segment column | Naturally encoded per-line |

Option (b) explicitly groups all branches of a stem under one line in the way
the format was designed: one line = one logical tree/stem. Downstream tools
that aggregate multi-stem trees by `tree_id` need only read the first field of
each data line, not walk per-segment attributes. When all stems share
`stem_id = 0`, the file differs from today's only in the first two fields of
each data line and the first two tokens of the header — backward-compatible with
any reader that does not inspect `treeAttributes()`.

### 4.3 Parser walk-through for extended format

`ForestStructure::load` (`rayforeststructure.cpp:93–284`), using the same
`1 + count(',')` field-counting convention throughout:

Header: `"tree_id,stem_id, x,y,z,radius,parent_id,section_id"` (note the
mandatory `, ` space separator between the per-tree prefix and `x,y,z,radius`,
matching the existing `trees.txt` writer convention at `raytrees.cpp:1228`).

- `mandatory_text = "x,y,z,radius"` found at position > 0
  (`rayforeststructure.cpp:120–126`)
- `lines[0]` = `header.substr(0, found-2)` = `"tree_id,stem_id"` (1 comma) →
  `commas_per_tree = 1 + 1 = 2`
- `lines[1]` = `header.substr(found)` = `"x,y,z,radius,parent_id,section_id"`
  (5 commas) → `commas_per_segment = 1 + 5 = 6`; after trim:
  `has_parent_id = true`; `attributes[1] = ["section_id"]`

For a **single-segment** data line (trunk-only stem):
`"42, 0, 394123.1,6012044.2,18.5,0.142,-1,0"` — comma count:
- 1 (between `42` and `0`) + 1 (between `0` and first coordinate) + 5
  (between the 6 segment fields) = **7 commas**
- `num_commas = 1 + 7 = 8`
- Validation: `(8 − 2) % 6 = 6 % 6 = 0` → **passes** ✓

For a **two-segment** data line:
`"42, 0, x,y,z,r,-1,0, cx,cy,cz,cr,0,1"` — 7 + 6 = **13 commas**:
- `num_commas = 1 + 13 = 14`
- Validation: `(14 − 2) % 6 = 12 % 6 = 0` → **passes** ✓

The format is fully compatible with `ForestStructure::load` for any number
of segments per stem.

### 4.4 Column table (extended format)

Per-line attributes (before `x,y,z,radius`):

| Column | Type written | Value |
|--------|-------------|-------|
| tree_id | int (written as-is) | Inventory tree ID; -1 never written here (valid stems only) |
| stem_id | int (written as-is) | Stem within tree; 0 = single-stem default |

Per-segment attributes (existing, unchanged):

| Column | Type | Value |
|--------|------|-------|
| x, y, z | float64 | Branch tip, world coords |
| radius | float64 | `Trees::radius(section)`, metres |
| parent_id | int | -1 for root; per-tree 0-based parent section id |
| section_id | int | `contiguous_section_ids_[sec]`, 0-based per stem |
| weight, len, accuracy, junction_weight | float64 | Verbose only |

### 4.5 Single-stem round-trip

A file written by the new pipeline with all `stem_id = 0`:
```
tree_id,stem_id, x,y,z,radius,parent_id,section_id
0, 0, 10.1,20.2,0.5,0.25,-1,0, 10.2,20.3,2.1,0.18,0,1, ...
1, 0, 15.3,22.1,0.4,0.31,-1,0, ...
```

This round-trips through `ForestStructure::load` with
`tree.treeAttributes() = [0.0, 0.0]` (tree_id=0, stem_id=0) for each tree.
Existing treetools code that reads segment geometry and ignores
`tree.treeAttributes()` is unaffected. Existing tests that compare moment
statistics (`getMoments()`) will see the same segment geometry; the per-tree
attributes do not affect geometry moments.

---

## 5. Mask loader output type

### 5.1 Both input forms produce the same type

```
struct MaskResult {
    SeedList   seeds;           // one StemSeed per unique (tree_id, stem_id)
    ray::Cloud labelled_cloud;  // empty Cloud{} in Level 1;
                                // populated in Level 2+ for spatial-prior use
    std::vector<int32_t> labelled_tree_ids;  // empty in Level 1
    std::vector<int32_t> labelled_stem_ids;  // empty in Level 1
};
```

`ray::Cloud{}` is a valid default-constructed empty cloud. The Level 1
implementation sets `labelled_cloud = {}` and leaves both label vectors empty.
The type is forward-compatible: Level 2 populates the cloud fields without
changing the interface.

### 5.2 Input form (a) — directory of per-stem PLY files

**Filename pattern**: `{tree_id}_{stem_id}.ply` or `{tree_id}.ply`.
- `stem_id` absent in filename → `stem_id = 0`.
- `tree_id` and `stem_id` are non-negative integers; the filename is the
  canonical source of (tree_id, stem_id) — not any in-file metadata.
- The filename parser must handle leading zeros (treat as decimal, not octal).
- Files not matching either pattern are skipped with a warning.

**StemSeed derivation from PLY points**:
1. Load each PLY as a bare point cloud (positions only; `cloud.ends`).
2. `base`: centroid of the lowest 10 % of points by height (z coordinate),
   after `removeStartPos()` relative to the *input cloud*, not the mask cloud.
   The mask cloud must first be translated to the same coordinate frame as
   the new cloud (they are assumed co-registered per CLAUDE.md).
3. `radius`: fit a cylinder to the lowest-`girth_height_ratio × tree_height`
   points; if fewer than 5 points, use half the bounding-box XY diameter.
4. `tree_height`: max(z) − base.z over all points in the file.

### 5.3 Input form (b) — single LAS/LAZ file

**Extra-bytes requirements**:
- `tree_id` extra-bytes field (int32) — required.
- `stem_id` extra-bytes field (int32) — optional; if absent, `stem_ids` all default to 0.
- Detection is by name (upstream_map §H, `raylaz.cpp:122–123`).

**StemSeed derivation from LAS**:
1. Load the LAS file; `cloud.tree_ids` and `cloud.stem_ids` are populated.
2. Group bounded rays by `(cloud.tree_ids[i], cloud.stem_ids[i])`.
3. For each group, derive `base`, `radius`, `tree_height` as in form (a) step 2–4,
   using points in the group.

### 5.4 Shared invariants for both forms

- The coordinate frame of all StemSeed fields must be local (after
  `removeStartPos()` applied to the **new input cloud**), not the mask cloud's
  own frame. The mask loader applies the same offset subtraction.
- If a file contains no points, a warning is printed and that (tree_id, stem_id)
  is excluded from the SeedList.
- Duplicate (tree_id, stem_id) pairs within the mask input are an error; the
  loader prints a fatal error and returns empty.

---

## 6. ID-collision policy when `--mask` is supplied

### 6.1 Matching mask stems to new-cloud stems

After `ray::segment()` produces its SeedList without a mask (or after the
constrained segmentation), each new-cloud stem is matched to a mask stem by
nearest-neighbour search on base XY position, subject to a radius threshold.

Match criterion: a new-cloud stem `S` matches mask stem `M` if:
- `dist_XY(S.base, M.base) ≤ match_radius` (default 2 × `max_diameter`)
- `|S.tree_height − M.tree_height| / M.tree_height ≤ 0.5`
  (height tolerance; prevents matching seedlings to mature trees)

One-to-one assignment: if multiple new-cloud stems map to the same mask stem,
the closest one wins; others are treated as unmatched.

**CLI flag**: `--match_radius <metres>` (default: `2 × params.max_diameter`)

### 6.2 Case 1 — matched stems

New-cloud stem inherits the mask's `(tree_id, stem_id)` exactly.
Output `tree_ids[]` and `stem_ids[]` use the mask values for all points
of this stem.

### 6.3 Case 2 — stems in new cloud not in mask (unmatched new)

Default: assign new tree_ids starting from `max(mask_tree_id) + 1`, incrementing
by 1 per new tree. `stem_id = 0` for unmatched new stems.

**CLI flag**: `--new_id_start <n>` — overrides the starting value.
If `n ≤ max(mask_tree_id)`, a warning is printed and the default
(`max(mask_tree_id) + 1`) is used.

### 6.4 Case 3 — mask stems with no matching points in new cloud (missing mask)

Default: omit from output. A warning summary is printed to stderr listing
each missing `(tree_id, stem_id)`.

**CLI flag**: `--keep_missing_stems` — if set, missing mask stems are written
to `cloud_trees.txt` with zero branches and a note in the comment header. They
appear in `cloud_seeds.txt` with their mask radius and base, so a downstream
merge tool can re-insert them.

**CLI flag**: `--missing_stems_file <path>` — writes the missing-stem summary
to a file (one `tree_id,stem_id` per line) instead of (or in addition to) stderr.
Useful for automated pipelines.

### 6.5 Collision-free guarantee

Because new IDs start from `max(mask_tree_id) + 1`, no new-cloud stem ID
can collide with any mask-carried ID. This is an invariant that must be
checked and asserted at the start of the output phase.

---

## 7. Backwards-compatibility proof

**Claim**: when called without `--mask` and without pre-segmented input, the
new pipeline produces output that is behaviourally identical to current
`rayextract trees`. "Behaviourally identical" means the geometry, radii, and
topology of `cloud_trees.txt` and `cloud_trees_mesh.ply` are the same; pixel
colours in `cloud_segmented.las` encode the same tree assignment; and
`cloud_trees.txt` round-trips through the same `ForestStructure::load` path.

### 7.1 How `rayextract trees` reduces to the new pipeline

`ray::trees()` is defined as:

```
result = ray::segment(cloud, offset, mesh, params, verbose);
return ray::reconstruct(cloud, offset, mesh, result.seeds, params, verbose);
```

No file I/O occurs between segment and reconstruct. The `Cloud` object and the
`SegmentResult` are passed by reference/value in memory, preserving the exact
floating-point state that the current monolithic `Trees` constructor maintains.

### 7.2 Phase-by-phase reduction (no-mask, no pre-segmented input)

| New phase | Equivalent old phase | Equivalence argument |
|-----------|---------------------|----------------------|
| `ray::segment` calls `getRootsAndSegment` | Phase 5 | Call is identical: same arguments, same function, same point data. Result: identical `points_[]` and `roots_list`. |
| stem assignment: `tree_ids[j] = root_to_tree_index[points_[j].root]`, `stem_ids[j] = 0` | Phase 15 `segmentCloud` sets `cloud.tree_ids[i] = contiguous_section_ids_[seg]` | In the no-mask case, each root_index maps 1-to-1 to one contiguous_section_id. The `tree_ids` assignment is a relabelling of the same Dijkstra-root grouping. ✓ |
| SeedList built: one StemSeed per root section, base and radius from trunk estimation loop | Phases 7–9 (generateRootSections + trunk estimation) | `ray::segment` runs phases 5–9 internally. StemSeed fields come from the same calculations: `sections_[i].tip` → base, `radius(sections_[i])` → radius, `sections_[i].tree_height` → tree_height. ✓ |
| `ray::reconstruct` re-runs local Dijkstra per (tree_id, stem_id=0) group seeded from StemSeed.base | Not applicable in old code — Dijkstra runs globally once in phase 5 | In the `ray::trees()` fast path (no serialisation), reconstruct receives the unmodified `points_[]` from segment, with parent links already set by the global Dijkstra. No re-run occurs. The local re-Dijkstra path is used only when reading back from files. For `rayextract trees`, the fast path ensures bit-identity. ✓ |
| Phase 10 branch reconstruction loop | Phase 10 unchanged | `ray::reconstruct` runs the same BranchSection loop with the same `params`. No changes to bifurcation, taper, or radius estimation code for the no-mask path. ✓ |
| Phase 15 `segmentCloud` | Phase 15 unchanged | Called with the same `sections_[]` and `section_ids[]`. `cloud.tree_ids` and `cloud.colours` are set identically. ✓ |
| trees.txt output: new header has `tree_id,stem_id` per-tree attributes prepended | Old header: `x,y,z,radius,parent_id,section_id` | The geometry (x,y,z,radius,parent_id,section_id) is identical. The header adds two new fields. Existing `ForestStructure::load` readers that do not inspect `treeAttributes()` are unaffected. Old treetools is unaffected. The segment geometry round-trips correctly. **Not byte-identical** to old output, but **behaviourally identical** per CLAUDE.md. This divergence is documented in CHANGELOG as "trees.txt now includes tree_id and stem_id per-line attributes." |
| cloud_segmented: new format is LAS with stem_ids all = 0 | Old format: same extension as input, no stem_ids | If input was `.ply`, old output was `_segmented.ply`; new output is `_segmented.las`. This is a **known behavioural difference** for `.ply` input users. Document in CHANGELOG: "cloud_segmented output is now .las to support tree_id/stem_id extra-bytes." Users who require PLY output can use `rayconvert` on the LAS. |
| cloud_trees_mesh.ply | Old: round-trips trees.txt through ForestStructure::load | Same round-trip, same geometry. The extra per-tree attributes in trees.txt are stored in `treeAttributes()` which `generateSmoothMesh` does not read. Output is geometrically identical. ✓ |

### 7.3 Summary of intentional divergences (for CHANGELOG)

1. **`cloud_trees.txt` header**: new header prepends `tree_id,stem_id,` before `x,y,z,radius`. Parsers that ignore per-tree attributes are unaffected. Parsers that do not tolerate the new header format must be updated.

2. **`cloud_segmented` format**: forced to `.las` when input is `.ply` (required for `stem_id` storage). Old `.ply` segmented output is no longer produced by `rayextract trees`.

3. **`cloud_seeds.txt`**: new file, not produced by old `rayextract trees`. Produced by new `rayextract segment`.

All other outputs (`cloud_trees_mesh.ply`, `cloud_shortest_paths.ply`) are unchanged.
