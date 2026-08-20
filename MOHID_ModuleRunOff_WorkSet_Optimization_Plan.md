# MOHID Land – ModuleRunOff Work-Set & Volume Optimization Plan (Phase 13+)

> Continuation of `MOHID_Performance_Optimization_Plan.md` (Phases 0–12, all DONE/committed on
> the `perf/Phase12` chain). This plan covers a **new set** of optimizations sourced from 11
> experimental branches. Update this file incrementally as phases land; other chat sessions
> continue from here.

## Context

- **Target file:** `Software/MOHIDLand/ModuleRunOff.F90` (~19,200 lines)
- **Hot path:** `HydrodynamicApproximation = DynamicWave_`, `GridIsConstant = .true.`, `LimitToCriticalFlow = .true.` → `DynamicWaveXX_default_CG` / `DynamicWaveYY_default_CG`
- **Build:** VS2017 / Intel Fortran (`ifxCompiler`), `Release Double OpenMP|x64`; profiling via `Profile Double OpenMP|x64` (see the Phase 0 section of the original plan for VTune setup)
- **Baseline for THIS plan:** confirmed **`perf/Phase12`** chain (Phases 1–12). New phases are cut from the confirmed+cleaned prior phase.
- **Test models:**
  - WithRain: `D:\Work\Projects\Performance2D\Benchmarks\UK_EPA_Benchmark_Test8A\Performance_LargeModel_quickversion\exe\`
  - NoRain: `...\Performance_LargeModel_quickversion_noRain\exe\`
- **Validation tool:** `compare_mohid.py` (repo root). `python compare_mohid.py <res dir>`. Known accepted FAIL: `FloodPeriod.dat`. HDF5/text outputs live one level up from `exe\` at `...\res\`.

---

## Source: Experimental Branch Analysis (Phase 0 — DONE)

Full report: `BranchAnalysis_ExperimentalBranches.md` (diffs of 11 branches vs `MohidLand_Bentley`, three-dot/merge-base). Branches forked from `MohidLand_Bentley`, **not** from the perf chain → porting requires reconciling with Phases 1–12 (Phase 4 unrolled face dispatch, Phase 10 merged X/Y `SetMatrixValue` pairs, Phase 11 cached `coeff²`, Phase 12 restricted the Courant scan).

### Branch verdicts (condensed)

| Branch | Theme | Verdict |
|---|---|---|
| `SetWorkSizeTests` | `BasinPointsWorkSize` pre-bound + `HasRainFall` reset move | promising (split into P13 + P14) |
| `ActivePointsTests` | 1D active-point list replaces masked `SetMatrixValue` | promising concept, incomplete (P15) |
| `BasinPointsTests` | 1D basin list for basin/output loops | safe but **low payoff here** (P16) |
| `MergeVelocityAdvection_NoWC_PreLoadVel` | fuse velocity+advection, preload neighbor velocities (`_SemWaterColumn_Vel`) | promising, incomplete — **furthest MergeVelocityAdvection variant** (P17) |
| `MergeVelocityAdvection_*` (others) | subsets/variants of the above | superseded / skip |
| `MergeVelocityAdvection_localgravity` | `localDT_x_Gravity` LICM hoist | marginal (compiler already does it) — skip |
| `SplitVolumeFromDynamicWave` | move volume accrual out of DynamicWave | incomplete (missing Y-flux, breaks `ActivePoints_Left`, debug `stop`s) → **P18 with fixes** |

### Verification of the analysis (done this session — corrections)

1. **`HasRainFall` is NOT a clean bug fix.** Today `ModifyRainFall_Infiltration` resets `Me%HasRainFall=.false.` at its end (`~L10791`) *before* `SetWorkSize` runs, so on rain steps `SetWorkSize` runs the shrink-scan (correctly picking up rained `ActivePoints=1` cells). Moving the reset makes rain steps use the **full domain + skip the scan** — a behavior change that can *increase* box size on partial-rain steps and shifts every `.not.HasRainFall`-gated routine. Treat as its own validated phase (P14), not a merge-immediately fix.
2. **Cluster A premise corrected.** XX and YY do **not** compute the same advection twice: `DynamicWaveXX_default_CG` builds X-face advection from `FlowXOld`/`AreaU` (`~L11661`); `DynamicWaveYY_default_CG` builds Y-face advection from `FlowYOld`/`AreaV` (`~L12874`). The real duplication is that `ComputeFaceVelocityModulus` **and** the DynamicWave routines each sweep `CurrentWorkSize` re-reading `FlowXOld/YOld`+`AreaU/V` and recomputing face velocities — that's what P17 fuses.
3. **`BasinPointsTests` low payoff for this model.** Its prime target, the double `ComputeCenterValues` (`~L16789`, full-`WorkSize` + `BasinPoints` guard), is **cold** (single-precision output uses `ComputeCenterValues_R4`, which already iterates `CurrentWorkSize`). Swapping `ComputeCenterValues_R4` to a static basin list would **regress** NoRain (basin ≫ active box). So P16 is safe but mostly for portability.
4. **Volume-accrual asymmetry explained (drives P18).** YY updates both `myWaterVolume(i,j)` and `(i-1,j)` directly (same-column neighbor → same OMP thread, `~L12909`). XX updates only `(i,j)` and defers the `(i,j-1)` side via `ActivePoints_Left` (`~L11692`) because the X neighbor is a *different column* → different thread → data race.

---

## Phased Plan

Branch strategy: continuous chain, each `perf/PhaseN` cut from the confirmed+cleaned `perf/Phase(N-1)`, starting from `perf/Phase12`. Validation methodology = the in-run A/B harness from the original plan (`CheckProfileScalarDiff` / `CheckProfileMatrixDiff`, tolerances `ProfileTolAbs_=1e-5`, `ProfileTolRel_=1e-4`).

### Phase 13 — `SetWorkSize` `BasinPointsWorkSize` pre-scan ✅ DONE  *(Sonnet, Low risk)*
- **From:** `SetWorkSizeTests` (the `BasinPointsWorkSize` parts only, NOT the `HasRainFall` move).
- **What:** Add `type(T_Size2D) :: BasinPointsWorkSize` to `T_RunOff`; compute a static basin bounding box once at construction; restrict `SetWorkSize`'s per-timestep `ActivePoints` scan (`~L11052`) to that box + early-return when already within ~5 cells of the basin extent.
- **Why safe:** `ActivePoints ⊆ BasinPoints ⊆ BasinPointsWorkSize`, so the resulting `CurrentWorkSize` is a result-safe superset; downstream loops (incl. Phase-12 Courant scan halo) stay bit-identical. Only `SetWorkSize`'s own cost drops (~7.9s NoRain wall).
- **Validate:** `compare_mohid` zero-diff both scenarios; `SetWorkSize` `ModuleStopWatch` CPU/Wall delta.
- **Result (2026-08-17):** ✅ `compare_mohid` no new failures both scenarios. `SetWorkSize` wall **NoRain 7.9s → 0.106s**, **WithRain → 0.012s** (CPU 0.000 both — routine now early-returns / scans only the tiny basin box). Cost essentially eliminated. Committed on `perf/Phase13`; cut `perf/Phase14` from here.

### Phase 14 — `HasRainFall` reset relocation ⏳ TODO  *(Opus, Medium — behavior-changing)*
- **From:** `SetWorkSizeTests` commit `0022d1d4`.
- **What:** Move `Me%HasRainFall=.false.` (and `HasInfiltration`) out of `ModifyRainFall_Infiltration` (`~L10791`) to after `ReadUnLockExternalVar` in `ModifyRunOff`'s outer block.
- **Risk:** rain steps switch to full-domain box + skip scan; interacts with `ModifyGeometryAndMapping`/output routines gated on `.not.HasRainFall`. Perf sign ambiguous (could grow the box on partial-rain steps).
- **Validate:** A/B + `compare_mohid` both scenarios; explicitly confirm **no WithRain regression**. May be dropped if net-negative.

### Phase 15 — `ActivePoints` 1D compact list ❌ REJECTED (dead end, reverted)  *(Opus, Medium — timing-sensitive)*
- **From:** `ActivePointsTests`.
- **What:** Add `ActivePointsI/J(:)`, `NumberOfActivePoints`, `ActivePointsNeedUpdate`; build via a `setActivePointsMapArray` scan; replace masked `CurrentWorkSize`+`ActivePoints` copy/update loops (e.g. Phase-10 `SetFlowOldXY` ~10.9s NoRain; `SetMatrixValue` sites `~L8030-L8114`) with 1D-indexed `SCHEDULE(STATIC)` loops.
- **Risk:** list goes stale after `UpdateWaterLevels`; the `ActivePointsNeedUpdate` rebuild discipline must be exact. Reconcile with Phase 10's merged X/Y passes.
- **Validate:** A/B (`CheckProfileMatrixDiff`) on affected arrays; `compare_mohid` both scenarios.

#### Phase 15 — RE-SCOPED (this session, Opus)  *(after code reading)*
The original "convert per-iteration `SetFlowOldXY` to a rebuild-on-mutation 1D list" premise **does not hold** for this model:
- `UpdateWaterLevels` mutates `ActivePoints` **every sub-iteration** — it *grows* it (`if (ActivePoints_Left(i,j+1)==1) ActivePoints(i,j)=1`, the X-flux the XX pass defers) and *shrinks* it (dry cells → `=0`) — at [ModuleRunOff.F90 L14005/L14029/L14041/L14063](Software/MOHIDLand/ModuleRunOff.F90#L14005).
- The 45.4 s consumer `SetFlowOldXY(...ActivePoints)` runs **once per sub-iteration** at [L8135](Software/MOHIDLand/ModuleRunOff.F90#L8135). Because `ActivePoints` changed in the prior iteration, the compact list would need `setActivePointsMapArray` (a serial box scan) **every iteration** — the same masked scan it replaces, but serial vs the current parallel (9.3×) loop → **no-op / regression**. The built-once/consumed-many precondition fails for the per-iteration consumer.

**Re-scoped to (user decision) → then measured → REJECTED:**
- **Part A (Option 2, implemented + measured, then reverted):** built the compact list **once per outer `ModifyRunOff` step**, while `ActivePoints` is still the stable state inherited from the previous step (right after `ModifyGeometryAndMapping`, *before* the `doIter` loop mutates it), used **only** for the once-per-step `ActivePoints`-masked consumers at [L8055–L8061](Software/MOHIDLand/ModuleRunOff.F90#L8055): `myWaterColumnOld` copy + `SetInitialFlowXY`, fused into one STATIC OMP region `SetInitialWorkSet`. Correct (**bit-identical by construction**; A/B harness `CheckProfileMatrixDiff8` **never tripped** over the full NoRain run) — but a **net performance loss**:

  | label (NoRain, 10T) | CPU s | Wall s |
  |---|---|---|
  | `SetInitialWorkSetBaseline` (old masked loops) | 115.9 | **15.64** |
  | `setActivePointsMapArray` (serial list rebuild) | 11.9 | **24.35** |
  | `SetInitialWorkSet` (1D-list consume) | 67.7 | **9.30** |

  New = rebuild + consume = 24.35 + 9.30 = **33.6 s** vs baseline **15.64 s** → **+18 s regression**. **Root cause (now measured):** the list can only be built by a serial box scan; `setActivePointsMapArray` runs single-threaded (11.9 CPU / 24.35 wall ≈ one memory-bound core) and *alone* costs more than the entire **parallel** (7.4×) masked baseline it feeds. The cheaper 1D consume (9.3 vs 15.6) can't recover the serial-scan tax. **This is the same wall as the per-iteration `SetFlowOldXY` case, and it also sinks the once-per-step consumers.** ⇒ code reverted to the Phase 13 baseline; `perf/Phase15` carries only this documented negative result. **The only way the 1D list becomes a win is P19 (incremental maintenance, no box scan) — see Candidate future phases.**
- **Part B (Option 3, pointer-swap for `SetFlowOldXY`) — INVESTIGATED → REJECTED (user decision, keep `SetFlowOldXY` as-is):** `lFlowX/Y` and `FlowXOld/YOld` are both `real(8) pointer` at identical `Me%Size` bounds → a module-level swap is mechanically trivial. **But it is not equivalent to the masked copy:** `DynamicWaveXX_default_CG` writes `lFlowX(i,j)` (value *or* `0.0`) for **every** `CurrentWorkSize` cell ([L11858/L11862](Software/MOHIDLand/ModuleRunOff.F90#L11858)), so after it a full swap sets `FlowOld` at recently-*dried* (now-inactive) cells to the last `lFlow` (≈0), whereas the masked copy leaves a **stale nonzero** `FlowOld` there — and advection reads `FlowXOld(i,j+1)` **unguarded** at [L11724](Software/MOHIDLand/ModuleRunOff.F90#L11724). ⇒ the swap perturbs the trajectory (Phase-1 `OpenPoints` cell-flip landmine, NoRain-sensitive). No bit-identical "avoid the copy" reformulation exists (the copy is inherent double-buffering). Parked like the cube-root landmine — revisit only under explicit accept-small-perturbation + long-run NoRain drift validation.

### Phase 16 — `BasinPoints` 1D list (optional, low payoff here) ⏳ TODO  *(Sonnet, Low risk)*
- **From:** `BasinPointsTests`.
- **What:** `BasinPointsI/J(:)` static list for full-`WorkSize` basin-guarded loops (double `ComputeCenterValues`, `firstRestart` resets). Bit-identical.
- **Note:** double `ComputeCenterValues` is cold in this model → do only if cheap or for other models/grids. **Do NOT** convert `ComputeCenterValues_R4` (would regress NoRain).

### Phase 17 — Fuse velocity+advection + neighbor-velocity preload ⏳ TODO  *(Opus, HIGH risk)*
- **From:** `MergeVelocityAdvection_NoWC_PreLoadVel` — extract `ComputeFaceVelocityModulus_SemWaterColumn_Vel` + `AdvectionTermU/V` cache; `DynamicWaveXX/YY_default_CG` read cached advection instead of re-traversing the stencil; preload neighbor face velocities (`flow/area`) into registers.
- **Mandatory pre-work:** prove the merged advection stencil is logically identical to the originals (`~L11661` / `~L12874`) — the report flagged an OR-vs-AND condition mismatch (`ComputeFaceV(i+1,j-1)+ComputeFaceV(i+1,j) > 0` vs `== 2`).
- **Drop** the branches' broken per-timestep-alloc harness entirely.
- **Validate:** full in-run A/B (`CheckProfileMatrixDiff` on `lFlowX`/`lFlowY`), both scenarios, watch NoRain `OpenPoints` cell-flip (Phase-1 landmine).

### Phase 18 — Decouple volume accumulation from DynamicWave (XX & YY) ⏳ TODO  *(Opus, HIGH risk — FINAL)*
- **From:** `SplitVolumeFromDynamicWave`, with bug fixes; extend to YY.
- **Concept:** move `myWaterVolume` accrual out of the flux loops into a per-cell **divergence pass** `Δvol(i,j) = (lFlowX(i,j) − lFlowX(i,j+1) + lFlowY(i,j) − lFlowY(i+1,j))·LocalDT`. Reads computed fluxes, writes only own cell → fully parallel; removes `ActivePoints_Left` and the XX cross-column race.
- **Main risk:** the in-loop `LimitToCriticalFlow` limiter reads `myWaterVolumePred` (predictor updated as flows are computed, `~L11700-L11712`). Preserve that — keep predictor updates in-loop and move only the *final* volume accrual, or prove divergence-form equivalence.
- **Fix branch showstoppers:** (a) add missing **Y-flux** term; (b) reconcile `ActivePoints_Left` removal with `UpdateWaterLevels` VG path; (c) remove debug `stop`s, restore `AlmostZeroNegative` tolerance.
- **XX first, then YY** (YY already race-free → smaller benefit; evaluate after XX).
- **Validate:** A/B snapshot `myWaterVolume`/`myWaterColumn`/`myWaterVolumePred`/`lFlowX`/`lFlowY`/`ActivePoints`; `CheckProfileMatrixDiff` on `myWaterVolume` AND `lFlowX`/`lFlowY` every step; both scenarios; `OpenPoints` cell-flip watch; then performant-only vs `_original`.

---

## Explicitly excluded
- `localDT_x_Gravity` hoist — compiler already does this LICM at `/O3`.
- `SplitVolumeFromDynamicWave` as-is — reconstructed in P18 (branch is incomplete).
- Converting `ComputeCenterValues_R4` to a static basin list — regresses NoRain.
- `MergeVelocityAdvection_{No_Simd,_2,NoWaterColumn,NoWaterColumn_PreLoad}` — subsets of P17's source.

## Ordering rationale
Cluster B (13–16) first: safest, genuinely hot (SetWorkSize), builds the 1D-list infrastructure. Then the two `DynamicWave` landmines (17, 18) last and **never in flight simultaneously**.

---

## Candidate future phases (not yet scheduled)

### P19 — Incrementally-maintained active/open 1D list  *(the real lever behind the P15 wall)*
Both the per-iteration `SetFlowOldXY` copy and the per-step output loops keep hitting the same wall: the compact list is **rebuilt by scanning the box** (`O(box)`), so a per-iteration/per-step rebuild costs the same as the masked scan it replaces (and serial vs the current parallel loop → net loss). **Measured (Phase 15, NoRain):** the serial `setActivePointsMapArray` rebuild = 24.35s wall, worse than the 15.64s parallel masked baseline it replaced → confirmed dead end for any scan-based rebuild. The breakthrough is to **maintain the list incrementally** — update it in place only when a cell toggles `0↔1`, which happens in exactly one per-iteration spot, `UpdateWaterLevels` (grow via `ActivePoints_Left`, shrink on dry-out at [ModuleRunOff.F90 L14005/L14029/L14041/L14063](Software/MOHIDLand/ModuleRunOff.F90#L14005)). Maintenance then becomes `O(changes)` instead of `O(box)`, which would turn the whole family (per-iteration `SetFlowOldXY` **and** the every-step output loops) into genuine NoRain wins.
- **Cost / risk:** `UpdateWaterLevels` is an OpenMP parallel loop → incremental append needs atomic or thread-local compaction, plus a periodic re-sort to keep the list monotone for cache-friendly `SCHEDULE(STATIC)`. Sizeable, HIGH-risk architectural change; validate with the in-run A/B harness (`CheckProfileMatrixDiff8` on the consumers) + NoRain cell-flip watch. Do standalone, never in flight with 17/18.

### Output-routine reuse of the 1D list — investigated, LOW priority (marginal here)
Rebuilding the list for the output routines is only *correct* if rebuilt **right before the output block** (after channel interaction / `RouteDFourPoints` / `Modify_Boundary_Condition` / discharges have finished mutating `ActivePoints`/`OpenPoints`), and scoped to the two loops with **no `else`** branch: `ComputeCenterVelocities_R4` non-distortion path and `OutputFloodingAll_R4`. **Excluded:** `ComputeCenterValues_R4` and the distortion branch of `ComputeCenterVelocities_R4` — their `else` zeroes just-dried cells, so a compact list would leave stale nonzero HDF/timeserie output. **Payoff is marginal:** the only every-step `OpenPoints`-guarded consumer (`OutputFloodingAll_R4`) is *single*, so a per-step box-scan rebuild does not beat its current masked scan; the multi-consumer benefit only appears on infrequent output steps (`ComputeCenterVelocities_R4` is output-step-only). `OutputFloodingAll_R4` also has a `Sum` reduction → list iteration changes the grouping (new FP nondeterminism; `FloodPeriod.dat` already the accepted-fail). ⇒ only worthwhile *after* P19 makes the list free to maintain.

---

## Running Notes / Findings (append as we go)

- 2026-08-11: Phase 0 (branch analysis) done + verified; corrections logged above. Working tree ≈ `MohidLand_Bentley` (Phase 11 `OverLandCoefficientXSquare` absent) — start new phases by checking out `perf/Phase12`.
- 2026-08-17: **Phase 13 DONE, no new failures.** `SetWorkSize` wall NoRain 7.9s → 0.106s, WithRain → 0.012s (see logs `LOG_WithRain_phase13.dat` / `LOG_NoRain_phase13.dat`). For reference, other NoRain wall costs at this point (10T): `ModifyGeometryAndMapping` 30.7s, `ComputeNextDT` 27.8s, `DynamicWaveXX_default_CG` 36.8s, `DynamicWaveYY_default_CG` 32.5s, `ComputeFaceVelocityModulus` 22.2s, `SetFlowOldXY` 45.4s, `UpdateWaterLevels` 17.7s, `OutputFloodingAll_R4` 13.5s, `ComputeCenterValues_R4` 15.5s. (`SetFlowOldXY` NoRain 45.4s wall is now a notable target — relevant to Phase 15's active-list.)
- 2026-08-20: **Phase 15 REJECTED — dead end, measured & reverted.** The compact 1D `ActivePoints` list is a net loss even for the safest once-per-step consumers: NoRain the serial `setActivePointsMapArray` rebuild is **24.35s wall alone** vs the **15.64s** parallel masked baseline it replaces (new total rebuild+consume = 33.6s → **+18s**). Confirmed root cause: list compaction is inherently serial (one memory-bound core) and loses to the 7.4×-parallel masked copy; no bit-identical parallel way to build it. A/B harness (`CheckProfileMatrixDiff8`) never tripped → the conversion was correct, just slower. `SetFlowOldXY` left at 45.6s (Part B pointer-swap rejected earlier: not bit-identical). Code reverted to Phase 13; `ModuleRunOff.F90` unchanged. **Next viable idea = P19 (incremental list maintenance).**
- _(add per-phase results, VTune deltas, gotchas here)_

## Files Modified (fill in as phases land)

| Phase | File | Change | Branch/commit |
|---|---|---|---|
| 13 | `Software/MOHIDLand/ModuleRunOff.F90` | Added `BasinPointsWorkSize` to `T_RunOff`; construction-time basin bounding box; `SetWorkSize` restricted-scan + early-return | `perf/Phase13` |
| 15 | — (reverted) | `ActivePoints` 1D compact list — implemented + measured, **net regression (+18s NoRain)**, code reverted to Phase 13 baseline | `perf/Phase15` (doc only) |
