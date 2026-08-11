# Experimental Branch Analysis — `ModuleRunOff.F90`
**Baseline:** `MohidLand_Bentley`  
**Date analysed:** 2026-08-11  
**Method:** `git diff --stat MohidLand_Bentley...<branch>` (three-dot, merge-base diff)

---

## 1. Summary Table

| Branch | Files changed | +lines / −lines | Theme | Correctness risk | Verdict |
|---|---|---|---|---|---|
| `DynamicWaveBreakAppart` | 1 | +633 / −14 | Loop fusion: merge velocity + advection into one pass | Medium-HIGH (advection logic changed; no validation) | **incomplete** |
| `MergevelocityAdvection_No_Simd` | 1 | +632 / −14 | Same as above minus `!$OMP SIMD` | None (SIMD removal) | **marginal** |
| `MergeVelocityAdvection_2` | 1 | +1240 / −175 | Same concept + dual-run A/B harness (water-column variant) | HIGH (harness bugs, per-step alloc) | **incomplete** |
| `MergeVelocityAdvection_NoWaterColumn` | 1 | +1235 / −175 | Same + removes water-column guard from advection | HIGH (comparison disabled, outer update commented) | **incomplete** |
| `MergeVelocityAdvection_NoWaterColumn_PreLoad` | 1 | +1240 / −173 | Activates `_SemWaterColumn` + `stop`-on-mismatch | HIGH (same production bugs + runtime `stop`) | **incomplete** |
| `MergeVelocityAdvection_NoWC_PreLoadVel` | 1 | +1540 / −173 | Adds `_SemWaterColumn_Vel`: preload 20+ neighbor velocities | HIGH (same harness) | **promising** concept, **incomplete** |
| `MergeVelocityAdvection_localgravity` | 1 | +1544 / −176 | Same + `localDT_x_Gravity` precompute in DynamicWaveXX | Low (bit-identical hoist) | **marginal** standalone |
| `SetWorkSizeTests` | 1 | +1592 / −184 | `BasinPointsWorkSize` pre-bound + `HasRainFall` timing fix | Low-Medium | **promising** |
| `ActivePointsTests` | 1 | +1999 / −364 | 1D active-point list replaces 2D masked `SetMatrixValue` | Medium | **promising** concept, **incomplete** |
| `BasinPointsTests` | 1 | +2081 / −488 | Extends 1D list to basin-point loops + output routines | Low-Medium | **promising** |
| `SplitVolumeFromDynamicWave` | 1 | +1629 / −178 | Move `myWaterVolume` accumulation out of DynamicWave | HIGH (mask change, debug `stop`) | **incomplete** |

---

## 2. Commit lineage

All eleven branches share a common ancestor in `MohidLand_Bentley`. The commit chain is:

```
MohidLand_Bentley
├─ DynamicWaveBreakAppart  (d39f3038 + 9225e0cc)
│    ├─ MergevelocityAdvection_No_Simd  (77fc5085 – removes !$OMP SIMD)
│    └─ MergeVelocityAdvection_2  (+ db6bd9d9, a5cb3f5b, df162601, 54cd5d07, 13dc077b)
│
└─ MergeVelocityAdvection_NoWaterColumn  (shares DWBreak base + d20809f0, 6ff0f185, …13dc077b)
     └─ MergeVelocityAdvection_NoWaterColumn_PreLoad  (+ 7a91db22 – switch to _SemWaterColumn)
          └─ MergeVelocityAdvection_NoWC_PreLoadVel  (+ bfeadcca, b07a1a12, 3c2a8c7e)
               ├─ MergeVelocityAdvection_localgravity  (+ f9c5d319 – localDT_x_Gravity)
               │    └─ SetWorkSizeTests  (+ 0022d1d4, 0202cb80, 5b151ed0)
               │         ├─ ActivePointsTests  (+ c8d5030a, 3be1751a, fa2542fa)
               │         │    └─ BasinPointsTests  (+ 2c5ed7c2, 665c6176)
               └─ SplitVolumeFromDynamicWave  (+ f9c5d319, 5e29b8ef, a27b257a, 769978be)
```

Note: `BasinPointsTests` (local) does not have a `remotes/origin/` counterpart; `SetWorkSizeTests` also has no remote. Both are local-only.

---

## 3. Detailed Subsections

---

### 3.1 `DynamicWaveBreakAppart`

**Diffstat:** `ModuleRunOff.F90` | +633 / −14 (647 total)

**Commits:**
- `d39f3038 MergeComputefacevelocity with advection`
- `9225e0cc Fix compilation`

**Intent:** The current production code computes the advection terms (`XLeftAdv`, `XRightAdv`, `YTopAdv`, `YBottomAdv`) inside both `DynamicWaveXX_default_CG` and `DynamicWaveYY_default_CG` independently, meaning the same stencil reads happen twice per timestep. This branch merges velocity magnitude computation AND advection term computation into a single `ComputeFaceVelocityModulus` pass, caching results in new persistent arrays `Me%AdvectionTermU(i,j)` / `Me%AdvectionTermV(i,j)`. The `DynamicWaveXX/YY_default_CG` routines then simply read `Me%AdvectionTermU(i,j) * LocalDT` instead of recomputing the advection stencil.

**What actually changed:**

| Location | Transformation |
|---|---|
| `T_RunOff` struct (L798) | Added `AdvectionTermU`, `AdvectionTermV` as `real(8), pointer` arrays |
| `AllocateVariables` | Allocates both arrays over full `Me%Size` |
| `ReadDataFile` | New `OPTIMIZATION` keyword → `Me%Optimization` (default `.false.`) |
| `ModifyRunOff` `DynamicWave_` case | Dispatches: `Optimization=.true.` → new `ComputeFaceVelocityModulus`; else → `ComputeFaceVelocityModulus_original` |
| `ComputeFaceVelocityModulus` (renamed to `_original`) | Unchanged logic, only renamed |
| New `ComputeFaceVelocityModulus` | Merged loop: computes `VelModFaceU/V` AND the full 4-neighbor advection terms into `AdvectionTermU/V`. Adds `!$OMP SIMD` on the inner `do i` loop |
| `DynamicWaveXX_default_CG` (old → `_original`) | Unchanged; renamed |
| New `DynamicWaveXX_default_CG` | Reads `Me%AdvectionTermU(i,j) * LocalDT` for advection instead of recomputing; also dispatches between `DynamicWaveXX_default_CG` / `_original` via `Me%Optimization` |
| `DynamicWaveYY_default_CG` | Same treatment (renamed, new version reads `AdvectionTermV`) |

**Correctness risk:**
- **Medium-HIGH.** The new `ComputeFaceVelocityModulus` uses slightly different advection conditions in places — e.g., the U-face Y-contribution checks use `Me%ComputeFaceV(i+1, j-1) + Me%ComputeFaceV(i+1, j) > 0` (OR condition) whereas the original DynamicWaveXX used `== 2` (AND condition) for that pair. Any logic mismatch will alter `lFlowX/Y`, perturbing `myWaterVolume` and potentially flipping `OpenPoints` cells.
- `!$OMP SIMD` on a loop body full of `if (Me%ComputeFaceU(i,j) == Compute)` branches is ineffective at best; no correctness issue but no gain either.
- `Me%Optimization` defaults to `.false.`, so **the new code is never exercised in production without the keyword** — the harness comparison is absent, so correctness is untested.

**Overlap with done work:** The advection precompute idea is independent of Phase 1–12. However, any implementation must clear the Phase-1 landmine check: the noRain `OpenPoints` cell-flip test.

**Verdict: `incomplete`.** The concept (merge velocity + advection into one pass) is architecturally sound and could meaningfully reduce the total cell iteration count. However, the new advection stencil logic is unverified, there is no A/B harness to validate it, and the `OPTIMIZATION` keyword defaults to off. Needs proper A/B integration before it is safe.

---

### 3.2 `MergevelocityAdvection_No_Simd`

**Diffstat:** +632 / −14 (646 total)

**Commits:** `DynamicWaveBreakAppart` + `77fc5085 Removes SIMD clause`

**Intent:** Test whether the `!$OMP SIMD` directive on the inner loop of the merged `ComputeFaceVelocityModulus` helps or hurts, by removing it.

**What actually changed:** Identical to `DynamicWaveBreakAppart` except the single line `!$OMP SIMD` is removed from the inner `do i` loop of the new `ComputeFaceVelocityModulus`.

**Correctness risk:** None — removing a SIMD hint never changes results.

**Overlap:** Same unvalidated advection logic as `DynamicWaveBreakAppart`.

**Verdict: `marginal`.** Answers one micro-question (SIMD yes/no) about `DynamicWaveBreakAppart`, but the underlying branch is itself unvalidated. If `DynamicWaveBreakAppart`'s concept is pursued, this variant confirms the `SIMD` directive should likely be dropped (branch prediction in the masked inner loop defeats vectorization anyway).

---

### 3.3 `MergeVelocityAdvection_2`

**Diffstat:** +1240 / −175

**Commits:** `DynamicWaveBreakAppart` base + `db6bd9d9 Get all matrix values needed beforehand`, `a5cb3f5b Prepare code for comparison`, `df162601 Allocates comparison matrixes`, `54cd5d07 added deallocate for local allocatables`, `13dc077b Fixes issues in comparison between routines`

**Intent:** Build a dual-run A/B validation harness on top of `DynamicWaveBreakAppart`. Each DynamicWave_ timestep: (1) run original `ComputeFaceVelocityModulus_original` + `DynamicWaveXX/YY_default_CG_original` + `UpdateWaterLevels` as baseline, saving state; (2) restore state; (3) run the new merged `ComputeFaceVelocityModulus` + new `DynamicWaveXX/YY_default_CG`; (4) compare results and `stop` on mismatch.

**What actually changed:**

| Location | Transformation |
|---|---|
| `ModifyRunOff` locals | Adds `Restart_entry`, 8 full-grid allocatables (`myWaterVolume/Column/lFlowX/Y _OriginalMethod/_Original`, `ActivePoints_*`) |
| `ModifyRunOff` DynamicWave_ case | Per-timestep: allocates 10 full-grid arrays, snapshots state, runs baseline path (`_original`), saves result, restores state, runs new path, comparison loop (commented-out `stop` in this branch) |
| `ModifyWaterDischarges` + `UpdateWaterLevels` outer calls | **Commented out** — these are moved inside the DynamicWave case only |
| `ComputeFaceVelocityModulus` | Uses the merged velocity+advection version from `DynamicWaveBreakAppart` (with water column check) |

**Critical bugs in production:**
1. **Per-timestep heap alloc/dealloc** of 8 full-grid `real(8)` arrays (~90 MB for a 1.54M-cell grid, every timestep) — extreme GC pressure and cache thrashing.
2. **`ModifyWaterDischarges` and `UpdateWaterLevels` commented out** in the outer `ModifyRunOff` loop — for non-DynamicWave flow types these routines are never called, breaking volume conservation.
3. **Double-execution of physics** (baseline + new path both update `myWaterVolume/Pred` before restoring): if state restore is incomplete, the second (performant) call corrupts flux state — exactly the ERR010 instability pattern documented in the repo memory.
4. Comparison block is commented out (no `stop`) — harness is running double physics for zero validation benefit.

**Correctness risk: HIGH.** Not runnable safely in production state.

**Verdict: `incomplete`.** The A/B harness design is sound in principle (mirroring the Phase 3–11 methodology) but this specific implementation has all four production-breaking bugs. The comparison block being commented out means it yields no correctness data even as a test.

---

### 3.4 `MergeVelocityAdvection_NoWaterColumn`

**Diffstat:** +1235 / −175

**Commits:** `MergeVelocityAdvection_2` base + `d20809f0 Removes waterColumn usage and only computes velU VelV`, `6ff0f185 Removed redundant code`

**Intent:** Introduce a variant of `ComputeFaceVelocityModulus` that drops the `Me%myWaterColumn` guard from the advection computation — instead of gating advection on `waterColumn_left > AlmostZero .and. waterColumn_right > AlmostZero`, the new variant (`ComputeFaceVelocityModulus` without "_Sem") computes advection for all active faces unconditionally. The motivation: avoid a potential stale-water-column read that exists when the velocity and advection computation precedes `UpdateWaterLevels` for that timestep.

**What actually changed:** Identical harness structure to `_2`. The new path calls the unguarded `ComputeFaceVelocityModulus` (no water-column gate). Comparison code is still commented out.

**Correctness risk: HIGH.** Same production bugs as `_2`. Additionally, removing the water-column guard from advection could activate advection for just-dried cells (those with `myWaterColumn > 0` from the last step but draining this step), which could increase instabilities at wet-dry fronts.

**Verdict: `incomplete`.** Explores one algorithmic variant of the concept but carries all parent harness bugs.

---

### 3.5 `MergeVelocityAdvection_NoWaterColumn_PreLoad`

**Diffstat:** +1240 / −173

**Commits:** `NoWaterColumn` + `7a91db22 Switch to new routine`, `bfeadcca Fixes missing variables in omp private`, `b07a1a12 Fixes velocity using closed faces`

**Intent:** Switch the new (performant) path from the unguarded `ComputeFaceVelocityModulus` to `ComputeFaceVelocityModulus_SemWaterColumn` — a new variant that explicitly handles closed faces and removes the per-cell water-column read for velocity magnitude (it reads `VelModFaceU` from already-computed values). Also enables the comparison block (`stop` on `|diff| > 1E-10`) — making this the first branch where the harness would actually abort on mismatch.

**What actually changed:**
- `ModifyRunOff` new path: `call ComputeFaceVelocityModulus_SemWaterColumn` (instead of plain `ComputeFaceVelocityModulus`)
- Comparison loop is **live** (not commented out), with `stop` on any cell exceeding `1E-10` absolute threshold
- Fixes for missing OMP `PRIVATE` variables and closed-face velocity handling

**Correctness risk:** HIGH. The active `stop` at `1E-10` is extremely tight — even a floating-point reordering from thread scheduling would trip it. A run that doesn't immediately abort would be strong evidence the results are numerically identical. But the other harness bugs (per-timestep alloc, double physics, commented-out outer `UpdateWaterLevels`) still apply.

**Verdict: `incomplete`.** The `stop` makes this actually usable as a correctness probe — but only if the harness bugs are fixed. In current state it is still not runnable in production.

---

### 3.6 `MergeVelocityAdvection_NoWC_PreLoadVel`

**Diffstat:** +1540 / −173

**Commits:** `NoWaterColumn_PreLoad` + `3c2a8c7e Adds pre-load of velocities to avoid duplicate multiplications`

**Intent:** Add a further-optimized variant `ComputeFaceVelocityModulus_SemWaterColumn_Vel` that **preloads** all 20+ neighboring face flows and pre-computes their velocities (`flow / area`) into local scalar registers before entering the advection stencil. For each active cell, this eliminates repeated `Me%FlowXOld(i±k,j±l)/Me%AreaU(i±k,j±l)` divisions — each neighbor's velocity is computed once and reused across multiple advection conditions.

**What actually changed:**
- New subroutine `ComputeFaceVelocityModulus_SemWaterColumn_Vel` (~300 lines): declares 20+ `real(8)` local scalars (`VelU`, `VelU_Left`, `VelU_Right`, etc.) + corresponding flow scalars; loads all of them from `Me%FlowXOld/YOld` and `Me%AreaU/V` once per cell, then uses the pre-loaded values throughout the advection computation
- Large OMP `PRIVATE` clause to cover all new locals
- `ModifyRunOff` new path: calls `ComputeFaceVelocityModulus_SemWaterColumn_Vel`
- Comparison loop still active with `stop`

**Correctness risk:** Same harness issues as parent. The preload technique itself is mathematically equivalent (same arithmetic, different load order) and is very unlikely to change results.

**Potential performance benefit:** Each `FlowXOld(i,j±n)/AreaU(i,j±n)` currently requires a pointer dereference + division. With preloading, 4–8 divisions are replaced by register reads within the advection stencil. This eliminates redundant loads visible both to the compiler and hardware prefetcher. This is the highest-quality idea in the MergeVelocityAdvection family.

**Verdict: `promising` concept, `incomplete` branch.** The preloading of neighbors is a genuine low-risk micro-optimization that could be extracted and applied cleanly. It does not change physics, only the order of memory loads. If separated from the harness overhead and validated with `compare_mohid.py`, it could deliver a 5–15% improvement in `ComputeFaceVelocityModulus`.

---

### 3.7 `MergeVelocityAdvection_localgravity`

**Diffstat:** +1544 / −176 (+4 vs `NoWC_PreLoadVel`)

**Commits:** `NoWC_PreLoadVel` + `f9c5d319 pre-compute some constants`

**Intent:** Hoist `LocalDT * Gravity` out of the per-cell loops in `DynamicWaveXX_default_CG` and `DynamicWaveXX_default_VG` (the `_VG` variant is also updated). Currently `LocalDT * Gravity` appears twice per cell (in `Pressure` and `Friction` calculations). Precomputing `localDT_x_Gravity = LocalDT * Gravity` once before the OMP region replaces two multiplications per cell with a single scalar read.

**What actually changed:**

| Location | Change |
|---|---|
| `DynamicWaveXX_default_CG` | `real :: localDT_x_Gravity`; computed before `!$OMP PARALLEL`; replaces `LocalDT * Gravity` in `Pressure` (L12583) and `Friction` (L12586) |
| `DynamicWaveXX_default_VG` | Same pattern |

Note: `DynamicWaveYY_default_CG/VG` are NOT updated in this branch (they receive the same treatment only in `SplitVolumeFromDynamicWave` and `SetWorkSizeTests` downstream).

**Correctness risk: Low.** `LocalDT * Gravity` is not modified inside the loop, so this is a pure LICM (Loop-Invariant Code Motion) transformation. Bit-identical on any conforming compiler.

**Overlap:** The Intel `ifx` compiler at `optimizeFull`/`/O3` almost certainly already performs this hoisting automatically. The gain is probably zero in a Release build.

**Verdict: `marginal` standalone.** Worth committing as a code-quality improvement (makes the loop intent clearer) but it won't show up in a VTune profile as a meaningful win given the compiler already does it.

---

### 3.8 `SetWorkSizeTests`

**Diffstat:** +1592 / −184

**Commits:** `localgravity` + `0022d1d4 Fixes bug(hasrainfall made false too early)`, `0202cb80 skip when all domain is already covered`, `5b151ed0 Update setworksize to use basinpoitsworksize as well`

**Intent:** Two independent improvements to the `CurrentWorkSize` bounding box mechanism:

1. **`BasinPointsWorkSize` pre-scan:** Compute a static bounding box around all `BasinPoints` cells once at construction. `SetWorkSize` then scans `ActivePoints` only within `BasinPointsWorkSize` (not the full `WorkSize`), and returns early if the current box is already within 5 cells of `BasinPointsWorkSize` on all sides.

2. **`HasRainFall` timing bug fix:** `Me%HasRainFall = .false.` was previously set at the end of `ModifyRainFall_Infiltration` (inside the `if (Me%Compute)` block). If `SetWorkSize` was called before rain processing finished in the same timestep, it could prematurely shrink `CurrentWorkSize` while cells were still being activated by rain. Moved to after `ReadUnLockExternalVar` in `ModifyRunOff`'s outer block.

**What actually changed:**

| Location | Change |
|---|---|
| `T_RunOff` struct | `type(T_Size2D) :: BasinPointsWorkSize` |
| `InitializeVariables` | Initialises `CurrentWorkSize` to an inverted (impossible) box to force the first scan; initialises `BasinPointsWorkSize` to full `WorkSize` |
| Construction path (after `ComputeNextDT`) | `call SetWorkSize(.true.)` — scans `BasinPoints`, stores in `BasinPointsWorkSize` |
| `SetWorkSize` | Optional `Constructor` argument; scans `BasinPoints` when `.true.`, `ActivePoints` otherwise; limits scan to `BasinPointsWorkSize` extent; `SkipILB/JLB/IUB/JUB` early-return if within 5 cells of basin extent |
| `ModifyRainFall_Infiltration` | `HasRainFall = .false.` and `HasInfiltration = .false.` removed |
| `ModifyRunOff` outer block | `HasRainFall = .false.` and `HasInfiltration = .false.` added after `ReadUnLockExternalVar` |

**Correctness risk: Low-Medium.** The `HasRainFall` fix is a genuine bug fix — the old placement could cause `SetWorkSize` to skip the scan on the first no-rain step, leaving `CurrentWorkSize` at the full basin box when it should start shrinking. The "skip if within 5 cells" optimization is safe because `BasinPointsWorkSize` is a tight outer bound; the only way `CurrentWorkSize` could have ILB < BasinPointsWorkSize.ILB is by a programming error elsewhere.

**Overlap with done work:** Phase 12 restricted `ComputeNextDT_CourantScan` to `CurrentWorkSize`. This branch makes `CurrentWorkSize` itself tighter for all NoRain scenarios and reduces the per-timestep cost of `SetWorkSize` itself. These are fully complementary.

**Verdict: `promising`.** Contains a genuine bug fix and a clean, generic work-reduction optimization. The `HasRainFall` timing fix alone is worth merging. The `BasinPointsWorkSize` pre-scan reduces `SetWorkSize`'s own scanning cost and benefits all consumers of `CurrentWorkSize`. This is the most immediately-mergeable branch after cleaning the harness overhead from the parent.

---

### 3.9 `ActivePointsTests`

**Diffstat:** +1999 / −364

**Commits:** `SetWorkSizeTests` + `c8d5030a First version before copilot makes his suggestion for improvements`, `3be1751a introduces flag to try and reduce call to setActivePointsMapArray`, `fa2542fa Apply new method to outputs where applicable`

**Intent:** Replace the masked 2D `SetMatrixValue(Me%myWaterColumnOld, Me%CurrentWorkSize, Me%myWaterColumn, Me%ActivePoints)` calls — which iterate every cell in `CurrentWorkSize` and check `ActivePoints(i,j)` — with direct iteration over a 1D compact list `Me%ActivePointsI(n)` / `Me%ActivePointsJ(n)` of only the active cells. Eliminates the dead-cell branch overhead and allows `SCHEDULE(STATIC)` (contiguous 1D array, uniform work).

**What actually changed:**

| Location | Change |
|---|---|
| `T_RunOff` struct | `ActivePointsI(:)`, `ActivePointsJ(:)`, `NumberOfActivePoints`, `NumberOfBasinPoints`, `ActivePointsNeedUpdate` |
| `AllocateVariables` | Allocates `ActivePointsI/J` to max possible size; builds `BasinPoints`-based initial list at construction |
| `ModifyRunOff` per-timestep | `setActivePointsMapArray` called when `ActivePointsNeedUpdate` and `.not.HasRainFall`; otherwise uses basin count |
| New `setActivePointsMapArray` subroutine | Scans `ActivePoints` over `CurrentWorkSize`, fills the 1D list |
| `ModifyRunOff` `if (Me%Compute)` block | Runs Method 1 (`SetMatrixValue`), saves result, runs Method 2 (1D indexed OMP loop with `STATIC`), compares with `stop` on mismatch > 1E-10 |
| `ActivePointsNeedUpdate = .true.` | Set whenever `ActivePoints` is restored in the harness (correct — list is stale after restore) |

**Correctness risk: Medium.** The 1D list approach is mathematically equivalent to the 2D masked loop. The main risk is timing: `setActivePointsMapArray` is called after `ComputeNextDT` but before the DynamicWave physics in each outer iteration — if `ActivePoints` changes within the same outer step (which it does, in `UpdateWaterLevels`), the list could be stale for part of the step. The `ActivePointsNeedUpdate` flag addresses this but the exact timing relative to `UpdateWaterLevels` needs care. The per-timestep alloc harness from parent still present.

**Overlap with done work:** Phase 10 halved the OMP fork/join count for `SetMatrixValue` by merging X/Y pairs. This branch eliminates the dead-cell scan entirely — a different and additional optimization. Both could coexist.

**Verdict: `promising` concept, `incomplete` branch.** The 1D compact active-point list is a clean idea that eliminates wasted work inside `SetMatrixValue` calls. If separated from the parent harness overhead and combined with Phase 10's X/Y merge, it could deliver meaningful gains for the `SetMatrixValue` calls that remain. The active A/B validation harness proves correctness in the comparison block.

---

### 3.10 `BasinPointsTests`

**Diffstat:** +2081 / −488

**Commits:** `ActivePointsTests` + `2c5ed7c2 Apply the Activepoints1D array methodology to basinpoints as well`, `665c6176 Apply the methodology to other setmatrixvalue calls`

**Intent:** Extend the 1D compact-list technique to **basin-point loops** (which are constant and never change at runtime). Specifically:
1. `BasinPointsI/J` — a permanent 1D list of all basin cells, built once at construction
2. The `firstRestart` path (`FlowXOld = InitialFlowX` over `BasinPoints`) is switched from `SetMatrixValue(... Me%ExtVar%BasinPoints)` to a direct `BasinPointsI/J` indexed loop with `SCHEDULE(STATIC)`
3. `ComputeCenterValues` and related output loops that currently do `do j; do i; if (BasinPoints(i,j)==1)` are rewritten as `do n = 1, NumberOfBasinPoints` with `i = BasinPointsI(n)` and `SCHEDULE(STATIC)`

**What actually changed:**

| Location | Change |
|---|---|
| `T_RunOff` | `BasinPointsI(:)`, `BasinPointsJ(:)` |
| `AllocateVariables` | Fills both `ActivePointsI/J` and `BasinPointsI/J` in the same construction loop |
| `ModifyRunOff` `firstRestart` block | `SetMatrixValue(Me%FlowXOld/YOld, ..., Me%ExtVar%BasinPoints)` → 1D indexed OMP DO SCHEDULE(STATIC) loop |
| `ComputeCenterValues` (double and R4 variants) | The `do j; do i; if (BasinPoints(i,j)==1)` loops replaced by `do n = 1, NumberOfBasinPoints; i = BasinPointsI(n); j = BasinPointsJ(n)`. `SCHEDULE(DYNAMIC, CHUNK)` → `SCHEDULE(STATIC)` |
| Similar to above in `ComputeCenterValues_R4` | Same transformation |

**Correctness risk: Low-Medium.** `BasinPoints` is static — the 1D list is exact. The `STATIC` schedule for uniform-work loops (each basin-point cell does the same computation) is correct and slightly more efficient than `DYNAMIC`. The transformation is logically equivalent to the existing 2D `if (BasinPoints)` loops.

**Overlap with done work:** Phase 7 applied `x**2. → x*x` to `ComputeCenterValues`/`_R4`. Both changes can be applied together — they are orthogonal.

**Verdict: `promising`.** The `BasinPointsI/J` technique is production-ready in principle — since `BasinPoints` is immutable, the list is built once and never stale. The `ComputeCenterValues` refactoring to `STATIC` + 1D indexing eliminates the 2D dead-cell scan and the `DYNAMIC` overhead for those loops. Needs to be separated from the parent harness scaffolding before merging.

---

### 3.11 `SplitVolumeFromDynamicWave`

**Diffstat:** +1629 / −178

**Commits:** `NoWC_PreLoadVel` + `f9c5d319 pre-compute some constants` (localDT_x_Gravity), `5e29b8ef Tests updating watervolume in updatewatervolumes`, `a27b257a removes activecell left check`, `769978be remove unnecessary var`

**Intent:** Restructure the volume bookkeeping: currently `DynamicWaveXX_default_CG` and `DynamicWaveYY_default_CG` each accumulate `Me%myWaterVolume` as they compute `lFlowX/lFlowY`. This branch moves all volume accumulation into a dedicated `UpdateWaterLevels_original` variant that computes the full divergence `(lFlowX(i,j) - lFlowX(i,j+1) + lFlowY(i,j) - lFlowY(i+1,j)) * LocalDT` in one pass.

Additionally: removes all `Me%ActivePoints_Left` assignments from `DynamicWaveXX_default_CG` (the `ActivePoints_Left` flag was set per X-face flow to enable `UpdateWaterLevels`'s non-critical path to track which cells received X flow).

**What actually changed:**

| Location | Change |
|---|---|
| `DynamicWaveXX_default_CG` | `dVol = lFlowX * LocalDT; myWaterVolume += dVol` removed from LimitToCriticalFlow path. All `Me%ActivePoints_Left(i,j) = 0/1` assignments removed. `CHUNK` local var removed (unused without ActivePoints_Left). |
| `DynamicWaveXX_default_VG` | Same `localDT_x_Gravity` precompute as localgravity |
| New inline `UpdateWaterLevels` (GridIsConstant path) | Added as a new subroutine: `myWaterVolume(i,j) += (lFlowX(i,j) - lFlowX(i,j+1)) * LocalDT` with debug `stop` on negative volume |
| A/B harness | Baseline: `UpdateWaterLevels_original`; performant: new `UpdateWaterLevels` |
| `ActivePoints_Left` comparison | Commented out in the A/B comparison loop |

**Correctness risk: HIGH.**
- Removing `ActivePoints_Left` from `DynamicWaveXX` means `UpdateWaterLevels`'s variable-grid path (which uses `if (Me%ActivePoints_Left(i,j+1) == 1)` to determine which cells need X-flow update) would receive all-zero `ActivePoints_Left` — silently producing wrong volumes in the VG path.
- The new `UpdateWaterLevels` GridIsConstant version has two debug `stop` statements: `stop` on negative volume (immediate crash on any dry-cell edge) and `stop` on slightly-negative volume (`< 0.0`) without the `AllmostZeroNegative` threshold. This makes it unusable in production — wet-dry fronts routinely produce small negative volumes that are flushed to zero by the tolerance.
- The Y-contribution (`lFlowY`) is absent from the new `UpdateWaterLevels` GridIsConstant code shown in the diff — only `lFlowX` is accumulated. This is either incomplete or a copy omission.

**Verdict: `incomplete`.** The architectural concept (decouple volume accumulation from flux computation) has merit for enabling future parallelism improvements, but the implementation has multiple showstopper bugs: missing Y-flux accumulation, breakage of `ActivePoints_Left` mask, and production-crashing `stop` statements.

---

## 4. Clusters & Recommendations

### Cluster A — Loop Fusion / Advection Precompute (`DynamicWaveBreakAppart` family)

**Branches:** `DynamicWaveBreakAppart`, `MergevelocityAdvection_No_Simd`, `MergeVelocityAdvection_2`, `MergeVelocityAdvection_NoWaterColumn`, `MergeVelocityAdvection_NoWaterColumn_PreLoad`, `MergeVelocityAdvection_NoWC_PreLoadVel`, `MergeVelocityAdvection_localgravity`

**Common idea:** Merge advection computation into `ComputeFaceVelocityModulus`, cache as `AdvectionTermU/V`, eliminate duplicate stencil traversal in `DynamicWaveXX/YY`.

**Comparison within cluster:**  
`DynamicWaveBreakAppart` → No_Simd (remove `SIMD`) → `_2` (add A/B harness) → `NoWaterColumn` (remove water-column gate) → `NoWaterColumn_PreLoad` (switch to `_SemWaterColumn`, enable `stop`) → `NoWC_PreLoadVel` (add preloaded neighbors — **furthest development**) → `localgravity` (add `localDT_x_Gravity` — trivial add-on)

The **furthest and most valuable variant is `MergeVelocityAdvection_NoWC_PreLoadVel`** with `ComputeFaceVelocityModulus_SemWaterColumn_Vel`: it contains both the advection precompute concept and the neighbor-velocity preloading. The `localDT_x_Gravity` hoist in `localgravity` is a trivial extra.

**Status of the whole cluster:** None of these branches are production-ready — they all carry the broken per-timestep-alloc harness. The harness must be stripped before any physics work can be evaluated.

**Recommendation:** Extract just `ComputeFaceVelocityModulus_SemWaterColumn_Vel` + `Me%AdvectionTermU/V` arrays + the `DynamicWaveXX/YY` new versions as isolated changes. Validate with `compare_mohid.py` on both scenarios (the noRain OpenPoints cell-flip is the critical test). If zero-diff, commit and profile. This is **Opus-led** territory (the advection stencil logic difference is a landmine; see Phase 1 history).

---

### Cluster B — Work-Set Restriction (`SetWorkSizeTests`, `ActivePointsTests`, `BasinPointsTests`)

**Common idea:** Reduce the number of cells visited per routine — either by a tighter bounding box (`SetWorkSizeTests`) or by replacing 2D masked loops with compact 1D active/basin-point lists.

**Comparison within cluster:**  
`SetWorkSizeTests` (bounding box, `HasRainFall` fix) → `ActivePointsTests` (1D active list for `SetMatrixValue`) → `BasinPointsTests` (1D basin list for output/reset routines)

**Ranking by risk/reward:**
1. **`SetWorkSizeTests` `HasRainFall` timing fix** — zero-risk bug fix, merge immediately.
2. **`SetWorkSizeTests` `BasinPointsWorkSize`** — clean, generic, complements Phase 12. Needs validation on both scenarios but expected zero-diff.
3. **`BasinPointsTests` basin 1D list** — `BasinPoints` is static so the list is always correct; `ComputeCenterValues` refactoring is safe. Expected zero-diff. Moderate scope.
4. **`ActivePointsTests` active 1D list** — Timing-sensitive (list must be rebuilt after `UpdateWaterLevels`); the `ActivePointsNeedUpdate` flag mechanism needs careful review. Medium complexity.

**Recommendation:** `SetWorkSizeTests` is the priority — the `HasRainFall` fix is a correctness issue. `BasinPointsTests` basin-list technique for `ComputeCenterValues` is clean and safe. Both are **Sonnet-appropriate** (bit-identical, straightforward transforms). `ActivePointsTests` active-list timing is subtler — review before committing.

---

### Cluster C — Algorithmic Restructuring (`SplitVolumeFromDynamicWave`)

**`SplitVolumeFromDynamicWave`** — moves volume accumulation out of `DynamicWaveXX/YY` and into `UpdateWaterLevels`. The concept enables cleaner X/Y decoupling and opens potential for further parallelism, but the current implementation has showstopper bugs (`ActivePoints_Left` removal breaks VG path; missing Y-flux; debug `stop`s). **Opus-led, not ready.**

---

### Priority ranking (investigate deeper)

| Priority | Branch / concept | Why | Who |
|---|---|---|---|
| 1 | `SetWorkSizeTests` — `HasRainFall` timing fix | Bug fix | Sonnet |
| 2 | `SetWorkSizeTests` — `BasinPointsWorkSize` bounding | Generic, complements Phase 12 | Sonnet |
| 3 | `BasinPointsTests` — 1D basin list for `ComputeCenterValues` | Low-risk, bit-identical | Sonnet |
| 4 | `MergeVelocityAdvection_NoWC_PreLoadVel` — extract `_SemWaterColumn_Vel` neighbor preload | High potential, needs landmine check | Opus |
| 5 | `ActivePointsTests` — 1D active list for `SetMatrixValue` | Medium complexity, timing-sensitive | Opus |
| 6 | `SplitVolumeFromDynamicWave` — volume restructuring | Interesting architecture, too many bugs | Opus (after fixes) |
| — | `DynamicWaveBreakAppart` / `No_Simd` / `_2` / `NoWaterColumn*` / `localgravity` | Superseded by / subsets of `NoWC_PreLoadVel` or trivial | Skip |

The `localDT_x_Gravity` optimization from `localgravity` (`f9c5d319`) is worth cherry-picking as a 2-line change — it is in `SetWorkSizeTests` and all downstream branches already.
