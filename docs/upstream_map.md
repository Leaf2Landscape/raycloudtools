# Upstream Map — raycloudtools segmentation/reconstruction refactor

Reference document for the segment/reconstruct split and multistem support work
described in CLAUDE.md. Cites file:line throughout. No code is proposed here.

---

## A. CLI dispatch map

### How subcommands are registered today

**File: `raycloudtools/rayextract/rayextract.cpp`**

| What | Line(s) | Detail |
|------|---------|--------|
| Subcommand name objects | 115 | `ray::TextArgument forest("forest"), trees("trees"), trunks("trunks"), terrain("terrain"), leaves("leaves"), grid("grid");` |
| parseCommandLine registrations | 163–179 | Six sequential calls, one per subcommand. Each returns `bool`; all six run unconditionally on every invocation. `parseCommandLine` (in `raylib/rayparse.cpp`) matches argv[1] against the first positional TextArgument and returns `true` only when all mandatory positional args match. |
| Guard (nothing matched) | 181–184 | `if (!extract_trunks && !extract_forest && ... && !extract_grid) { usage(); }` |
| Dispatch if-else-if chain | 187–365 | `if (extract_trunks)` → `else if (extract_trees)` → `else if (extract_forest)` → `else if (extract_terrain)` → `else if (extract_leaves)` → `else if (extract_grid)` → `else { usage(true); }` |
| `main()` | 368–371 | `return ray::runWithMemoryCheck(rayExtract, argc, argv);` — wraps `rayExtract()` for leak detection |

Subcommand handler line ranges:
- `trunks`: 187–199
- `trees`: 201–295
- `forest`: 298–334
- `terrain`: 339–349
- `leaves`: 351–354
- `grid`: 356–359

### Minimum changes to add `segment` and `reconstruct`

Four mechanical edits to `rayextract.cpp`, no other file touched:

1. **Line 115** — append two TextArgument declarations:
   ```cpp
   ray::TextArgument segment("segment"), reconstruct("reconstruct");
   ```

2. **After line 179** — add two `parseCommandLine` calls following the same pattern as the others. `segment` takes `cloud_file` as a mandatory positional arg (mesh is optional via `--ground`). `reconstruct` takes `cloud_file` and `mesh_file` as mandatory positionals.

3. **Line 181** — extend the guard condition with `&& !extract_segment && !extract_reconstruct`.

4. **After line 295** — insert two `else if (extract_segment)` and `else if (extract_reconstruct)` blocks before the existing `else if (extract_forest)` block (or immediately before `else { usage(true); }` — ordering only matters for readability).

5. **Lines 27–104 (usage function)** — add `if (extract_type == "segment" || none)` and `if (extract_type == "reconstruct" || none)` blocks with usage strings.

No changes to `raylib/CMakeLists.txt`, `rayparse.cpp`, or the build system are needed for the dispatch machinery itself. New library functions (`ray::segment`, `ray::reconstruct`) will be compiled into `raylib` and linked automatically.

---

## B. Phase list for `rayextract trees`

Numbered from input load to last output write.

| # | Phase | Entry point (file:line) | One-sentence description |
|---|-------|------------------------|--------------------------|
| 1 | Load cloud | `rayextract.cpp:204–209` | `Cloud::load()` into memory; `removeStartPos()` subtracts `ends[0]` from all starts and ends in-place and returns the offset vector |
| 2 | Load mesh | `rayextract.cpp:211–216` | `readPlyMesh()` into `ray::Mesh`; mesh translated by `-offset` |
| 3 | Params setup | `rayextract.cpp:218–271` | Populate `TreesParams` from parsed CLI arguments |
| 4 | `Trees` constructor entry | `raytrees.cpp:34` | `Trees::Trees(cloud, offset, mesh, params, verbose)` begins |
| 5 | `getRootsAndSegment` | `raytrees.cpp:39–40` | Dijkstra from ground-mesh vertices through all bounded cloud points; fills `points_[]` (parent, root, distance_to_ground, score, weight) and returns `roots_list` (groups of mesh-vertex indices, one group per detected tree) |
| 6 | `calculatePointDistancesToEnd` | `raytrees.cpp:46` | Bottom-up propagation along parent links: each Vertex accumulates `distance_to_end = max(distance over all descendant paths)` |
| 7 | `generateRootSections` | `raytrees.cpp:49` | Creates one `BranchSection` per element of `roots_list`; sets `sections_[i].root = i`; sets `max_distance_to_end` per root from its mesh-vertex points |
| 8 | Build `children[]` | `raytrees.cpp:52–58` | Inverts `points_[i].parent` links into a forward adjacency list `children[parent] → {children}` |
| 9 | Trunk radius estimation loop | `raytrees.cpp:63–211` | For each root section (`sec_` = 0 … roots\_list.size()-1): uses `girth_height_ratio` to set a measurement window; calls `extractNodesAndEndsFromRoots` + `estimateCylinderRadius`; sets `sections_[sec_].taper`, `total_taper`, `total_weight`, `tree_height`, `forest_taper_` |
| 10 | Branch reconstruction loop | `raytrees.cpp:191–274` | Iterates all sections including those added during this loop; BFS upward from roots; detects bifurcations via `findPointClusters`; calls `estimateCylinderRadius` and `estimateCylinderTaper` per section; appends child sections |
| 11 | `calculateSectionIds` | `raytrees.cpp:277–278` | Assigns each Vertex a `section_id` (index into `sections_[]`) based on which BranchSection owns its path |
| 12 | `filterLargestDiameterTrees` | call `raytrees.cpp:284`; defined inline `raytrees.h:195` | Optional (`--largest_diameter`): clears `children` of all non-max-DBH root sections to suppress them from output |
| 13 | `generateLocalSectionIds` | `raytrees.cpp:287` | Assigns per-tree 0-based `id` to each `BranchSection` in breadth-first order; prints "N trees saved" |
| 14 | `removeOutOfBoundSections` | `raytrees.cpp:291–294` | Optional (`--grid_width`): clears `children` of root sections whose base lies in the grid overlap zone |
| 15 | `segmentCloud` | `raytrees.cpp:296–298` | Writes `cloud.tree_ids[]` and `cloud.colours[]` (RGB via `convertIntToColour`); builds `contiguous_section_ids_[]`; writes `root_segs[]` |
| 16 | `removeOutOfBoundRays` | `raytrees.cpp:300–302` | Optional (`--grid_width`): removes cloud rays whose root section is outside the non-overlapping cell |
| 17 | `trees.save()` | `rayextract.cpp:275` | Writes `cloud_trees.txt` (piecewise cylindrical representation, one line per tree) |
| 18 | `saveShortestPaths()` | `rayextract.cpp:278–280` | Optional (`--save_paths`): writes `cloud_shortest_paths.ply` (ASCII PLY, vertex + edge elements) |
| 19 | `cloud.translate(offset)` + `cloud.save()` | `rayextract.cpp:283–284` | Re-adds offset; writes `cloud_segmented.ply` or `cloud_segmented.las` with per-point tree colours and `tree_ids` |
| 20 | `forest.load` + `generateSmoothMesh` + `writePlyMesh` | `rayextract.cpp:288–295` | Round-trips `_trees.txt` through `ForestStructure::load`; generates `cloud_trees_mesh.ply` via smooth capsule meshing |

---

## C. State table

Columns: phase → reads (state in) → writes (state out) → private intermediates produced and discarded.
Member variable names are exact.

| Phase | Reads | Writes | Private intermediates |
|-------|-------|--------|----------------------|
| 5 `getRootsAndSegment` | `cloud.ends`, `cloud.starts`, `cloud.colours[i].alpha` (if `alpha_weighting`), `mesh.vertices()` | `points_[]` fields: `pos`, `start`, `parent`, `root`, `distance_to_ground`, `score`, `weight`; `roots_list` (return value) | `QueueNode` priority queue; KD-tree (`Nabo::NNSearchD`); `heightfield`, `lowfield` 2-D arrays; `counts`, `sums`, `bests`, `max_heights` — all discarded on return |
| 6 `calculatePointDistancesToEnd` | `points_[].parent`, `points_[].pos` | `points_[].distance_to_end` | — |
| 7 `generateRootSections` | `roots_list`, `points_[].distance_to_end` | `sections_[]`: `roots`, `root`, `max_distance_to_end` (initial root sections only) | — |
| 8 Build `children[]` | `points_[].parent` | `children[]` (`std::vector<std::vector<int>>`, local to constructor) | — |
| 9 Trunk estimation | `sections_[sec_].roots`, `points_[]`, `children[]`, `params_->girth_height_ratio`, `params_->girth_height_ratio * tree_height` | `sections_[sec_]`: `taper`, `total_taper`, `total_weight`, `tree_height`, `tip`, `len`, `accuracy`, `junction_weight`, `radius_scale`, `ends`, `children`; `forest_taper_`, `forest_weight_`, `forest_weight_squared_` | `nodes`, `best_nodes`, `best_ends` (local vectors per iteration) |
| 10 Branch reconstruction | `sections_[sec_].*`, `points_[]`, `children[]` | `sections_[]` extended with new child sections; per-section: `tip`, `taper`, `total_taper`, `total_weight`, `len`, `accuracy`, `split_count`, `ends`, `children` | `clusters` from `findPointClusters`; `nodes` per section; `children` (lambda-local copy for `bifurcate`) |
| 11 `calculateSectionIds` | `sections_[].roots` (mesh-vertex indices), `children[]` | `section_ids[]` (`std::vector<int>`, size = `points_.size()`, local to constructor) | — |
| 13 `generateLocalSectionIds` | `sections_[].children`, `sections_[].parent` | `sections_[].id` (per-tree 0-based integer) | — |
| 15 `segmentCloud` | `section_ids[]`, `points_[].root`, `sections_[].*`, `params_->segment_branches` | `cloud.tree_ids[]`; `cloud.colours[]` (RGB); `contiguous_section_ids_[]` (member); `root_segs[]` (local, passed to phase 16) | — |
| 17 `trees.save()` | `sections_[]`, `contiguous_section_ids_[]`, `offset` | `_trees.txt` on disk | — |
| 19 `cloud.save()` | `cloud.ends`, `cloud.starts`, `cloud.times`, `cloud.colours`, `cloud.tree_ids`, `cloud.passthrough`, `cloud.extra_bytes_vlr` | `_segmented.ply` or `_segmented.las` on disk | — |

---

## D. Segment / reconstruct split point

### Desired outputs

Per CLAUDE.md:

- **Segmentation** → `cloud_segmented.las` (per-point `tree_id`, `stem_id`) + `cloud_seeds.txt` (per-stem trunk metadata)
- **Reconstruction** → `cloud_trees.txt` + `cloud_trees_mesh.ply`

### Where the split falls in the current phase list

The cleanest conceptual split is between phase 5 (`getRootsAndSegment`) and phase 7 (`generateRootSections`): after segmentation every bounded point has a `tree_id` (via `points_[j].root`) but reconstruction has not begun.

In practice the split must be exposed as the boundary between two library functions `ray::segment()` / `ray::reconstruct()`, as specified in CLAUDE.md.

### Dirty intermediates that cross the boundary

If reconstruction is given only `(tree_id, stem_id)` per point and per-stem seed positions, it cannot proceed without re-deriving the following:

| Intermediate | Role in reconstruction | Re-derivable from segmentation outputs? |
|--------------|----------------------|----------------------------------------|
| `points_[j].parent` | BFS direction in `extractNodesAndEndsFromRoots` and `extractNodesFromEnds` (phases 9–10) | **Yes** — re-run `connectPointsShortestPath` within each `(tree_id, stem_id)` point set, seeded from the stem base in `seeds.txt`. Equivalent cost to the original Dijkstra but scoped per stem. |
| `points_[j].distance_to_end` | Crop threshold (`crop_length`), taper scale, section height ordering (phases 6, 9, 10) | **Yes** — recomputed bottom-up from re-derived parent links in O(N). |
| `points_[j].root` | Tree-group index used in `segmentCloud` (phase 15) | **Trivially** — set to any canonical integer per `(tree_id, stem_id)`. |
| `children[]` local array | Forward adjacency (phase 8 onward) | **Yes** — inverted from parent links in O(N). |

**Conclusion**: No intermediate must be persisted. `ray::reconstruct()` can re-derive all path state by running a constrained local Dijkstra per `(tree_id, stem_id)` group, seeded from the stem position in `seeds.txt`. This matches the CLAUDE.md invariant that "segmentation outputs must be sufficient inputs for reconstruction."

There is no pre-existing phase boundary in the current `Trees` constructor that maps to a clean file-I/O checkpoint; the split requires factoring the constructor into the two library functions.

---

## E. girth_height_ratio

### Every site where it is read or set

| Site | File:line | Read or set | Detail |
|------|-----------|-------------|--------|
| Declaration | `raytrees.h:25` | declare | `double girth_height_ratio; // how far up tree to measure girth` |
| Default | `raytrees.cpp:18` | set | `0.12` in `TreesParams::TreesParams()` |
| CLI argument type | `rayextract.cpp:131` | declare | `ray::DoubleArgument girth_height_ratio(0.001, 0.5)` |
| CLI option object | `rayextract.cpp:140` | declare | `ray::OptionalKeyValueArgument girth_height_ratio_option("girth_height_ratio", 'i', &girth_height_ratio)` |
| CLI usage string | `rayextract.cpp:60` | display | `"--girth_height_ratio 0.12 - (-i) the amount up tree's height to estimate trunk girth"` |
| CLI assign to params | `rayextract.cpp:237` | set | `params.girth_height_ratio = girth_height_ratio.value();` |
| **Algorithm use (single site)** | **`raytrees.cpp:81–82`** | **read** | `double girth_height = params_->girth_height_ratio * tree_height;` |

### What the algorithm use gates

`girth_height_ratio` has exactly one algorithmic site: `raytrees.cpp:81–82`, inside the trunk estimation loop (phase 9), executed once per root section.

```cpp
// raytrees.cpp:81–82
double girth_height = params_->girth_height_ratio * tree_height;
```

`girth_height` then drives:

1. **Trunk measurement window** (`raytrees.cpp:89–93`): three attempts with `max_dist = girth_height * j/2.0` for j=1,2,3, spanning `[0.5 × girth_height, 1.5 × girth_height]`. Each attempt calls `extractNodesAndEndsFromRoots(nodes, base, children, max_dist*2/3, max_dist)` to gather cloud points for radius fitting.

2. **Trunk radius estimate** (`raytrees.cpp:~106–148`): `estimateCylinderRadius` is called on each node set; the attempt with lowest estimated radius is selected as `best_accuracy` / `best_nodes`. This sets `sections_[sec_].taper` via `estimateCylinderTaper`.

3. **Global taper accumulation** (`raytrees.cpp:~149–182`): the per-trunk taper feeds `forest_taper_`, `forest_weight_`, `forest_weight_squared_`, which in turn govern `meanTaper()` — the function used to scale ALL branch radii throughout phase 10 (`raytrees.cpp:347–370`).

### Code paths that must be bypassed or made conditional when (tree_id, stem_id) is supplied

When stem separation is pre-defined, two problems arise:

1. **Radius estimation gathers points from adjacent stems.** `extractNodesAndEndsFromRoots` at `raytrees.cpp:89–93` follows `children[]` upward from the root mesh vertices without any stem-boundary check. If two stems share overlapping BFS traversal ranges, points from the wrong stem enter the radius fit. **Bypass required**: constrain the BFS to points whose `stem_id` matches the current section's stem, or run reconstruction within per-`(tree_id, stem_id)` point sets from the start.

2. **A pre-supplied seed radius should replace or prime the estimate.** If `seeds.txt` supplies a trusted radius for a stem, the three-attempt estimation loop (`raytrees.cpp:89–148`) should either be skipped (use seed radius directly) or treat the seed radius as a prior. The relevant block is `raytrees.cpp:89–148` (the `for (int j = 1; j<=3; j++)` loop and subsequent best-selection logic).

Concretely, the code regions to make conditional on "seed radius available":
- `raytrees.cpp:89–148` — radius estimation loop (replace with seed value)
- `raytrees.cpp:149–182` — taper estimation from `best_nodes` (still needed but must accept an externally-supplied `rad`)
- No bypass needed for `raytrees.cpp:183–211` (tip placement and child section creation)

---

## F. trees.txt schema

**Writer**: `raytrees.cpp:1218–1269` (`Trees::save`)
**Parser**: `rayforeststructure.cpp:93–284` (`ForestStructure::load`)

### File layout

```
# Tree file. Optional per-tree attributes (e.g. 'height,crown_radius, ') followed by 'x,y,z,radius' and any additional per-segment attributes:
x,y,z,radius,parent_id,section_id[,weight,len,accuracy,junction_weight]
<tree-1 root seg>, <tree-1 child seg>, <tree-1 child seg>, ...
<tree-2 root seg>, ...
```

- **Comment line** (`raytrees.cpp:1227`): begins with `#`; stored in `ForestStructure::comments` by the parser.
- **Header line** (`raytrees.cpp:1228–1232`): always `x,y,z,radius,parent_id,section_id`; if verbose also `,weight,len,accuracy,junction_weight`.
- **One data line per tree** (`raytrees.cpp:1234–1266`): all branch segments for one tree concatenated, segments separated by `, ` (comma-space). Root segment is first; remaining segments are in breadth-first order.

### Per-segment fields

| Column | Writer line | Type written | Value |
|--------|-------------|-------------|-------|
| x | `raytrees.cpp:1241, 1255` | float (offset re-added) | `section.tip[0] + offset[0]` |
| y | same | float | `section.tip[1] + offset[1]` |
| z | same | float | `section.tip[2] + offset[2]` |
| radius | same | double | `Trees::radius(section)` — taper-estimated cylinder radius |
| parent_id | `raytrees.cpp:1241` (root: literal `-1`); `raytrees.cpp:1255` (child: `sections_[node.parent].id`) | int | `-1` for root; 0-based per-tree id of parent section otherwise |
| section_id | `raytrees.cpp:1241, 1255` | int | `contiguous_section_ids_[sec]` — 0-based per-tree, set in `generateLocalSectionIds` |
| weight | `raytrees.cpp:1232` | double | verbose only; `section.weight` |
| len | `raytrees.cpp:1232` | double | verbose only; `section.len` |
| accuracy | `raytrees.cpp:1232` | double | verbose only; `section.accuracy` |
| junction_weight | `raytrees.cpp:1232` | double | verbose only; `section.junction_weight` |

### Parser column handling (`rayforeststructure.cpp:93–284`)

- **Mandatory substring**: header must contain `"x,y,z,radius"` (`rayforeststructure.cpp:120–126`); failure is fatal.
- **Pre-mandatory text** = per-tree attributes (`rayforeststructure.cpp:129–151`): parsed as doubles into `tree.treeAttributes()`.
- **Post-mandatory attributes** (`rayforeststructure.cpp:153–176`): the string `"parent_id"` is recognised by name and sets `has_parent_id = true` (`rayforeststructure.cpp:167–170`); all other named attributes are stored in `attributes[1]` (branch attribute names) and parsed into `segment.attributes[]` as doubles (`rayforeststructure.cpp:250`).
- **Comma-count validation** (`rayforeststructure.cpp:206–210`): `(num_commas − commas_per_tree) % commas_per_segment == 0`. This allows multiple segments per line but requires all data lines to have a comma count consistent with the header. A fatal error is returned if the check fails.
- **Unknown columns silently accepted** and stored in `segment.attributes[]` (`rayforeststructure.cpp:250`). With the current header `x,y,z,radius,parent_id,section_id`, `section_id` is stored in `segment.attributes[0]`.

### Encoding stem_id

**Recommended — append as a named per-segment attribute**:

Extend header to `x,y,z,radius,parent_id,section_id,stem_id`. The parser stores `stem_id` in `segment.attributes[1]` as a double (integer value). All existing code that reads trees.txt and ignores `segment.attributes` continues to work unchanged. To promote `stem_id` to a typed field in `TreeStructure::Segment`, only treetools (or any tool that explicitly iterates `segment.attributes`) needs updating; that change is isolated.

This is strictly more compatible than adding a per-tree attribute (prepending before `x,y,z,radius`), which would require the forest extractor and any tool that assumes one line = one inventory tree to be updated.

---

## G. cloud_trunks.txt schema

**Writer**: `raytrunks.cpp:678–698` (`Trunks::save`)

```
# tree trunks file:
x,y,z,radius
<x>, <y>, <z>, <radius>
...
```

- **Comment header** (`raytrunks.cpp:686`): `# tree trunks file:`
- **Column header** (`raytrunks.cpp:687`): `x,y,z,radius`
- **Data format** (`raytrunks.cpp:695`): comma-space delimited (`", "`); 4 fields; only `trunk.active == true` entries written
- **z** = trunk base height: `base = trunk.centre - trunk.dir * trunk.length * 0.5` (`raytrunks.cpp:694`)

### Does `Trunks::load()` tolerate extra columns?

**No. `raytrunks.cpp:717–739` is strict: extra columns cause a fatal return.**

```cpp
// raytrunks.cpp:717–718
const int num_commas = std::count(line.begin(), line.end(), ',');
if (num_commas == 3)   // exactly 4 fields
{ ... parse ... }
else
{
    std::cerr << "bad input, there should be 4 fields per line: x, y, z, radius." << std::endl;
    return std::vector<...>();   // raytrunks.cpp:739 — empty return, caller treats as error
}
```

### Does `rayextract forest --trunks` call `Trunks::load()`?

**No.** `rayextract.cpp:320–331` — the `--trunks` path calls `ForestStructure::load()` (`rayextract.cpp:323`), not `Trunks::load()`. The loaded `ForestStructure::trees` is then iterated to extract `segment.tip` and `segment.radius` only (`rayextract.cpp:329`).

```cpp
// rayextract.cpp:320–331
if (trunks_option.isSet())
{
    ray::ForestStructure forest;              // local ForestStructure, not ray::Forest
    if (!forest.load(trunks_file.name()))    // uses ForestStructure::load
        usage(true);
    for (auto &tree : forest.trees)
        trunks.push_back({tree.segments()[0].tip, tree.segments()[0].radius});
}
```

### Can the seeds file be a strict superset of cloud_trunks.txt?

**Yes, for the `forest --trunks` use case.** `ForestStructure::load()` stores unknown post-mandatory columns in `segment.attributes[]` (`rayforeststructure.cpp:250`) and the forest extractor discards all attributes, using only `tip` and `radius` (`rayextract.cpp:329`).

A seeds file with header `x,y,z,radius,tree_id,stem_id` and data lines with 5 commas each:
- Passes `ForestStructure::load()` validation (`commas_per_segment = 6`, one segment per line, consistent across all lines) ✓
- `tree_id` and `stem_id` parsed as doubles into `segment.attributes[0]` and `segment.attributes[1]` ✓
- `rayextract forest --trunks` ignores them and uses only `tip`/`radius` ✓

**The strict `Trunks::load()` would still break** on a seeds file with extra columns, but `Trunks::load()` is not in the forest pipeline and does not need to be updated.

---

## H. PLY / LAS reader / writer for labelled clouds

### PLY — no extra-field support today

| Fact | File:line |
|------|-----------|
| `RayPlyEntry` is a compile-time-sized struct: 9 floats (`RAYLIB_DOUBLE_RAYS=OFF`) or 12 floats (`RAYLIB_DOUBLE_RAYS=ON`) | `rayply.h:17` (DOUBLE_RAYS branch), `rayply.h:20` (default branch) |
| `writePlyRayCloud()` signature has no `tree_ids` parameter | `rayply.cpp:165–188` |
| `Cloud::save()` for PLY passes no `tree_ids` | `raycloud.cpp:40–41` |
| PLY header is hard-coded (9 properties) | `rayply.cpp:37–80` |
| PLY reader matches only hard-coded property names; unknown properties silently skipped in row-size accumulation but not extracted | `rayply.cpp:330–770` |

**tree_id and stem_id cannot be stored in PLY without modifying `RayPlyEntry`, the header writer, the body writer, and the reader.** This is a non-trivial change. A simpler path is to mandate `.las` output for the labelled segmented cloud from `rayextract segment`.

### LAS/LAZ — full extra-bytes support today for tree_id

| Fact | File:line |
|------|-----------|
| `tree_id` registered as int32 EXTRA_BYTES attribute named `"tree_id"`, type 5 | `raylaz.cpp:673` |
| Byte layout in extra_bytes: sx(0–3), sy(4–7), sz(8–11), tree_id(12–15), alpha(16) when tree_id present; alpha at 12 when absent | `raylaz.cpp:850–854`, `raylaz.cpp:268` |
| Written as `int32_t` at `point_->extra_bytes + 12` | `raylaz.cpp:850–854` |
| VLR detection is **by name** (scans EXTRA_BYTES VLR for `"tree_id"`) | `raylaz.cpp:92–134`, specifically `raylaz.cpp:122–123` |
| Read back as `int32_t` from byte offset 12; pushed to `tree_ids_out` | `raylaz.cpp:183–188` |
| `Cloud::save()` for LAS passes `tree_ids` through | `raycloud.cpp:38–39` |
| If `cloud.tree_ids` is empty, LAS is written without the `tree_id` attribute and alpha stays at byte 12 | `raylaz.cpp:668–674` (conditional registration) |

### Adding stem_id to LAS

The same pattern used for `tree_id` applies directly:

| Step | Where to add |
|------|-------------|
| Register int32 attribute `"stem_id"` | `raylaz.cpp:673` (after tree_id registration) |
| Write at `point_->extra_bytes + 16` | `raylaz.cpp:854` area — alpha moves to byte 20 |
| Detect `"stem_id"` in VLR scan | `raylaz.cpp:122–123` area — add `has_stem_id_attr` flag |
| Read back from byte 16 into `stem_ids_out` | `raylaz.cpp:183–188` area |
| Add `std::vector<int32_t> stem_ids` to `Cloud` | `raycloud.h:41` (after `tree_ids`) |
| Thread through `Cloud::save()` → `writeLasRayCloud()` | `raycloud.cpp:39` |

### Field naming convention

Attribute names live in the EXTRA_BYTES VLR records. Current raycloud-owned attributes (in order): `"sx"`, `"sy"`, `"sz"`, `"alpha"`, `"tree_id"` (`raylaz.cpp:95`). New names `"stem_id"` can follow the same convention; they are distinguished from original-sensor extra-bytes by membership in this known-names list (`raylaz.cpp:95–101`).

### Summary: what needs a writer change

| Format | tree_id today | stem_id addition | Change scope |
|--------|--------------|-----------------|--------------|
| PLY | Not stored | Not stored | Significant: struct, header, body writer, reader |
| LAS/LAZ | Stored, round-tripped | ~20 lines following existing pattern | Minimal: extend existing extra-bytes machinery |
