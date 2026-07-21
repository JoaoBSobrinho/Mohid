# MOHID Land – ModuleRunOff Performance Optimization Plan

## Context

- **Target file:** `Software/MOHIDLand/ModuleRunOff.F90` (~17,662 lines)
- **Hot path:** `Me%HydrodynamicApproximation = DynamicWave_`, `Me%GridIsConstant = .true.`, `Me%LimitToCriticalFlow = .true.`
  - Routes to `DynamicWaveXX_default_CG` → `DynamicWaveYY_default_CG`
- **Build:** Visual Studio 2017 / Intel Fortran (`ifxCompiler`), `Release Double OpenMP|x64`
  - `optimizeFull`, `ipoSingleFile`, `OpenMPParallelCode`, `RealKIND="realKIND8"`, `FloatingPointModel="precise"`, `FloatingPointExceptionHandling="fpe0"`
- **VTune:** `C:\Program Files (x86)\Intel\oneAPI\vtune\latest\bin64\vtune.exe`
- **Profile build exe:** `Solutions\VisualStudio2017_IntelFortran18\MOHIDNumerics\MOHIDLand\x64\Profile Double OpenMP\MOHIDLand.exe`
- **Test models:**
  - WithRain: `D:\Work\Projects\Performance2D\Benchmarks\UK_EPA_Benchmark_Test8A\Performance_LargeModel_quickversion\exe\`
  - NoRain: `D:\Work\Projects\Performance2D\Benchmarks\UK_EPA_Benchmark_Test8A\Performance_LargeModel_quickversion_noRain\exe\`

---

## VTune Profiling Infrastructure (DONE)

- `Profile Double OpenMP|x64` build config added to:
  - `Solutions/VisualStudio2017_IntelFortran18/MOHIDNumerics/MOHIDLand/MOHIDLand.vfproj` (line ~184)
  - `Solutions/VisualStudio2017_IntelFortran18/MOHIDNumerics/MOHIDBase1/MOHIDBase1.vfproj` (line ~195)
  - `Solutions/VisualStudio2017_IntelFortran18/MOHIDNumerics/MOHIDBase2/MOHIDBase2.vfproj` (line ~195)
  - `Solutions/VisualStudio2017_IntelFortran18/MOHIDNumerics/MOHIDNumerics.sln`
- Profile config settings: `ipoOff`, `WholeProgramOptimization="false"`, `debugEnabled`, `GenerateDebugInformation="linkDebuggingInformation"` (needed for VTune symbol resolution)
- VTune batch scripts (sw hotspot mode):
  - `...\Performance_LargeModel_quickversion\exe\run_vtune_hotspot.bat`
  - `...\Performance_LargeModel_quickversion\exe\run_vtune_threading.bat`
  - `...\Performance_LargeModel_quickversion_noRain\exe\run_vtune_hotspot.bat`
  - `...\Performance_LargeModel_quickversion_noRain\exe\run_vtune_threading.bat`

---

## VTune Profiling Results (Baseline)

### Top hotspots – WithRain 10T vs 1T

| Function | 10T CPU (s) | 1T CPU (s) | Speedup | Root Cause |
|---|---|---|---|---|
| `_for_ieee_signaling_gt/lt/ge_k8_` | ~250 combined | — | — | real(8) comparisons under fpe0 |
| `_libm_pow_l9` | 164 | — | — | `**2.` and `**(4./3.)` in friction; `abs(cmplx())` |
| `for_is_nan_t_` | 150 | — | — | `abs(cmplx(U,V))` complex modulus + sqrt under fpe0 |
| `COMPUTENEXTDT` | ~90 | ~185 | ~2x | Courant scan, sqrt under fpe0 |
| `_kmpc_barrier` | — | — | low | OMP load imbalance (NoRain case) |
| `MODIFYGEOMETRYANDMAPPING` | — | ~107 | 1.2x | Memory bandwidth / access pattern |
| `SETMATRIXVALUES2D_R8_FROMMATRIX` | — | ~119 | 1.5x | Memory bandwidth |

---

## Phase 0 – VTune Profiling Infrastructure ✅ DONE

## Phase 1 – Code-Level Optimizations in ModuleRunOff.F90 ✅ DONE (see results below)

All changes applied to `Software/MOHIDLand/ModuleRunOff.F90`:

### 1a. `ComputeFaceVelocityModulus` – eliminate `abs(cmplx(...))` ✅
```fortran
! Before:
Me%VelModFaceU(i, j) = abs(cmplx(U, Vaverage))
Me%VelModFaceV(i, j) = abs(cmplx(Uaverage, V))

! After:
Me%VelModFaceU(i, j) = sqrt(U*U + Vaverage*Vaverage)
Me%VelModFaceV(i, j) = sqrt(Uaverage*Uaverage + V*V)
```
**Eliminates:** `_libm_pow_l9` from complex modulus path. **Does NOT eliminate `for_is_nan_t_`** — explicit `sqrt()` under `fpe0` calls the same NaN-check wrapper. `for_is_nan_t_` moved source but did not decrease (see Phase 1 Results). Root fix is Phase 2 FPE relaxation.

### 1b. `DynamicWaveXX_default_CG` + `DynamicWaveYY_default_CG` – friction formula ✅
```fortran
! Before (2 pow calls):
Me%OverlandCoefficientX(i,j)** 2. / (HydraulicRadius ** (4./3.))

! After (1 pow call, n*n replaces n**2):
Me%OverlandCoefficientX(i,j) * Me%OverlandCoefficientX(i,j) / &
    (HydraulicRadius * HydraulicRadius ** (1.0/3.0))
```
**Eliminates:** one `_libm_pow_l9` call per cell (the `**2.` call). Reduces from 2 to 1 pow call.  
**Note:** `cbrt()` is NOT a Fortran intrinsic in Intel ifx — use `HydraulicRadius ** (1.0/3.0)` instead.

### 1c. `DynamicWaveXX_default_CG` – CriticalFlow real(8) comparisons ✅
```fortran
! Declarations changed:
!   real(8) :: CriticalFlow  →  real :: CriticalFlow, lFlowX_local
! OMP PRIVATE clause: added lFlowX_local

! Before (real(8) comparisons → _for_ieee_signaling_gt_k8_):
if (abs(Me%lFlowX(i, j)) > CriticalFlow) then
    if (Me%lFlowX(i, j) > 0) then
        Me%lFlowX(i, j) = CriticalFlow
    else
        Me%lFlowX(i, j) = -1.0 * CriticalFlow

! After (real(4) comparisons):
lFlowX_local = real(Me%lFlowX(i, j), 4)
if (abs(lFlowX_local) > CriticalFlow) then
    if (lFlowX_local > 0.0) then
        Me%lFlowX(i, j) = CriticalFlow
    else
        Me%lFlowX(i, j) = -CriticalFlow
```
**Eliminates:** `_for_ieee_signaling_gt/lt_k8_` calls in the LimitToCriticalFlow path.

### 1d. `DynamicWaveYY_default_CG` – same as 1c with `lFlowY_local` / `Me%lFlowY` ✅

> **Scope:** Only `_default_CG` variants (GridIsConstant=true path). The `_default_VG`, `_CG`, `_VG` variants are unchanged.

---

## Phase 1 – Results

### WithRain – 10 threads

| Function | Original (s) | After Phase 1 (s) | Δ |
|---|---|---|---|
| `DynamicWaveYY_default_CG` | 292 | 216 | **-26%** |
| `DynamicWaveXX_default_CG` | 276 | 215 | **-22%** |
| `_for_ieee_signaling_gt_k8_` | 186 | 151 | **-19%** |
| `_libm_pow_l9` | 164 | 134 | **-18%** |
| `for_is_nan_t_` | 150 | 160 | **+7% (regression!)** |
| `ComputeFaceVelocityModulus` | 104 | 69 | **-34%** |
| `_kmpc_barrier` (spin) | 40 | 29 | **-28%** |
| `_for_ieee_signaling_lt_k8_` | 36 | 29 | **-19%** |
| `ModifyGeometryAndMapping` | 77 | 63 | -18% |
| `ComputeNextDT` | 62 | 53 | -15% |
| `SetMatrixValues2D_R8_FromMatrix` | 62 | 52 | -16% |

> **Estimated wall-time reduction (WithRain 10T): ~15–20%**

### WithRain – 1 thread

| Function | Original (s) | After Phase 1 (s) | Δ |
|---|---|---|---|
| `DynamicWaveXX_default_CG` | 184 | 141 | **-23%** |
| `DynamicWaveYY_default_CG` | 183 | 137 | **-25%** |
| `_for_ieee_signaling_gt_k8_` | 138 | 107 | **-22%** |
| `_libm_pow_l9` | 118 | 97 | **-18%** |
| `for_is_nan_t_` | 114 | 129 | **+13% (regression!)** |
| `ComputeFaceVelocityModulus` | 78 | 48 | **-38%** |
| `for_is_nan_s_` | 0 | 12 | **new (real4 sqrt path)** |

### NoRain – 10 threads

| Function | Original (s) | After Phase 1 (s) | Δ |
|---|---|---|---|
| `DynamicWaveXX_default_CG` | 195 | 204 | +5% (marginal regression) |
| `DynamicWaveYY_default_CG` | 207 | 190 | -8% |
| `ComputeNextDT` | 185 | 185 | **0% — unchanged** |
| `_kmpc_barrier` (spin) | 149 | 150 | ~0% |
| `_for_ieee_signaling_gt_k8_` | 110 | 111 | ~0% |
| `_libm_pow_l9` | 94 | 93 | ~0% |
| `for_is_nan_t_` | 94 | 97 | ~0% |
| `SetMatrixValues2D_R8_FromMatrix` | 119 | 116 | -3% |

> **Phase 1 had minimal effect on NoRain case** — the `LimitToCriticalFlow` path (changes 1c/1d) is not hot in this scenario. ComputeNextDT completely unchanged.

### Root-cause analysis of `for_is_nan_t_` regression

> **NOTE: The arithmetic reformulations described here were REVERTED — see Step 1 above.** The `abs(cmplx(U,V))` → `sqrt(U*U+V*V)` change and the friction reformulations caused a noRain regression (`OpenPoints` cell flip). Phase 1 source changes are now limited to the 4 `max()` guards only, which are no-ops and produce zero diff vs original.

Under the old (reverted) Phase 1a: replacing `abs(cmplx(U, V))` with `sqrt(U*U + V*V)` moved the NaN-check but didn't eliminate it — every explicit `sqrt()` call under `FloatingPointExceptionHandling="fpe0"` goes through the same NaN-checking wrapper. This also explains why `for_is_nan_s_` (real4 variant) appeared as a new hotspot. **Phase 2 addresses this properly at compiler level.**

### Summary table – remaining costs (WithRain 10T)

| Category | Combined CPU time (s) | Root cause |
|---|---|---|
| `for_is_nan_t_` + ieee_signaling variants | ~340 | `fpe0` wrapping of every sqrt/comparison |
| `DynamicWaveXX` + `DynamicWaveYY` | ~431 | Physics loop — further reducible |
| `_libm_pow_l9` | ~134 | `**(1.0/3.0)` cube root in friction term |
| `ComputeNextDT` | ~53 | CFL scan, sqrt per cell |
| `ModifyGeometryAndMapping` | ~63 | Memory bandwidth |
| `SetMatrixValues2D_R8_FromMatrix` | ~52 | Memory bandwidth |

**The single biggest remaining opportunity is eliminating `fpe0` overhead (~340s out of ~1300s remaining = 26% of all hotspot CPU time).**

---

## Phase 2 – FPE Mode Relaxation for Hot Subroutines (STEP 2a IN PROGRESS — awaiting validation)

**Background:** `fpe0` (trap on all floating-point exceptions) forces every `sqrt`, every real(8) comparison, and every `pow` through library wrappers that check for NaN/Inf before returning. This is the root cause of:
- `for_is_nan_t_` — wraps every `sqrt(real8)` call
- `for_is_nan_s_` — wraps every `sqrt(real4)` call
- `_for_ieee_signaling_gt/lt/ge_k8_` — wraps every `>`, `<`, `>=` comparison on `real(8)`
- Partially explains the cost of `_libm_pow_l9` (pow is also wrapped under fpe0)

**Expected gain (WithRain 10T): ~340s = ~26% of remaining hotspot CPU time eliminated.**

> **IMPLEMENTATION NOTE FOR NEXT AGENT:** Steps 0, 1, and 2a are already done (see below). Start from Step 2b. All edits are in `Software/MOHIDLand/ModuleRunOff.F90`. Re-grep line numbers before editing — file is ~17,670 lines and earlier edits shift them.

---

### Step 0 — Validate Phase 1 correctness ✅ DONE

Ran `compare_mohid.py` against both test cases:
- WithRain: 18/19 PASS. `FloodPeriod.dat` FAIL (9.3e-2 s timing shift) is known and accepted.
- NoRain `_phase1` vs current: HDF5 = exact 0 diff.

---

### Step 1 — Add `max(..., 0.0)` sqrt guards ✅ DONE (final, validated)

Applied 4 defensive guards to `Software/MOHIDLand/ModuleRunOff.F90`. These are purely defensive no-ops when physics is valid (`Gravity * depth >= 0` always) — confirmed zero numerical diff vs original in both WithRain and noRain cases.

| Location | Change |
|---|---|
| `DynamicWaveXX_default_CG` (CriticalFlow) | `sqrt(Gravity * WaterDepth)` → `sqrt(max(Gravity * WaterDepth, 0.0))` |
| `DynamicWaveYY_default_CG` (CriticalFlow) | `sqrt(Gravity * WaterDepth)` → `sqrt(max(Gravity * WaterDepth, 0.0))` |
| `ComputeNextDT` (GridIsConstant path) | `sqrt(Gravity * aux)` → `sqrt(max(Gravity * aux, 0.0))` |
| `ComputeNextDT` (variable-grid path) | `sqrt(Gravity * aux)` → `sqrt(max(Gravity * aux, 0.0))` |

> **`ComputeFaceVelocityModulus` NOT changed:** `abs(cmplx(U, Vaverage))` retained. An earlier attempt replaced it with `sqrt(max(U*U + Vaverage*Vaverage, 0.0))` and also reformulated the friction term (`coeff**2.` → `coeff*coeff`, `HydraulicRadius**(4./3.)` → `HydraulicRadius * HydraulicRadius**(1./3.)`). These changes caused the noRain case to fail vs original (`Grid/OpenPoints` cell flip at output timestep 5 — one borderline cell tipped across the active threshold by floating-point accumulation). **All arithmetic reformulations were reverted.** Only the four `max()` guards above survive. The `abs(cmplx())` and `**(4./3.)` friction formula are unchanged from the original.

---

### Step 2a — Eliminate IEEE-comparison overhead via `/assume:noieee_compares` ✅ DONE and VALIDATED

**Validation results:**
- IEEE wrappers (`_for_ieee_signaling_*`, `for_is_nan_t_`, `for_is_nan_s_`) **completely gone** from VTune profile (1-thread WithRain hotspot — see `hotspot_report_WithRain_phase2_step1_Test_1Thread_2.txt`).
- **WithRain:** 36/37 PASS vs `_original`. Only `FloodPeriod.dat` fails (pre-existing from before any optimisation). All files pass vs `_phase1` with zero diff.
- **NoRain:** ALL PASS vs `_original` (zero diff — Phase 1 `max()` guards are no-ops).
- No NaN/Inf anomalies in either case.

> **ROOT CAUSE — CORRECTED.** The `for_is_nan_t_` / `for_is_nan_s_` / `_for_ieee_signaling_gt/lt/ge_*` hotspots (~345s combined) were **NOT** caused by `/fpe:0`. Proof: the `Profile Double OpenMP` build log shows `ModuleRunOff.F90` compiled with **no `/fpe` switch** (i.e. default `/fpe:3`) yet the wrappers were still present. The real cause is **`/standard-semantics`** (emitted by `F2003Semantics="true"`), which implicitly enables **`/assume:ieee_compares`**. That option routes every `real(8)` relational operator through IEEE-correct library routines (`_for_ieee_signaling_*`), each of which calls `for_is_nan_t_` to pre-check operands for NaN. `/fp:precise` is NOT the cause and is kept. `_libm_pow_l9` is genuine cube-root compute (`**(1.0/3.0)`), not FPE/IEEE overhead.

> **Abandoned dead ends (kept for history):**
> - `!DEC$ OPTIONS /fpe=3` per-routine directive — invalid; Intel Fortran ignores it. Removed.
> - `fpe0 → fpe3` alone — necessary but **not sufficient**. It does not remove the wrappers (they come from `ieee_compares`, not fpe). It is still applied because it lets `/assume:noieee_compares` be safe: hardware compares can raise INVALID on a stray NaN, and `fpe3` (no trapping) prevents that from aborting.
> - Per-file `FileConfiguration` overrides on `ModuleRunOff.F90` — **ignored** under `MultiProcessorCompilation="true"` (`/MP` compiles sources in a shared batch). Only **project-level** compiler options take effect. (`ModuleRunOff.F90` is built by `MOHIDLand.vfproj`; `MohidLandEngine.vfproj` only has `OpenMI` configs and is not built by the Double OpenMP solution config.)

> **CRITICAL — VS must reload external edits.** The `.vfproj` edits were made directly in XML. If Visual Studio has the solution open, it builds from its in-memory copy until the project is reloaded. **Close & reopen VS (or Unload → Reload project), then Build → Clean Solution → Rebuild.** A flag-only change may not trigger incremental recompilation of `ModuleRunOff.F90`.

**Applied (both settings, at the configuration level):** `FloatingPointExceptionHandling="fpe3"` **and** `AdditionalOptions="/assume:noieee_compares"`, for both `Release Double OpenMP|x64` and `Profile Double OpenMP|x64`, in all three projects:
- `MOHIDLand.vfproj` (Release L174, Profile L185)
- `MOHIDBase1.vfproj` (Release L185, Profile L196)
- `MOHIDBase2.vfproj` (Release L185, Profile L196)

Debug configs intentionally left at `fpe0` for development error-trapping.

**Branch:** `PerformanceTestingCopilot` (working tree, uncommitted). `perf/phase1-only` has Phase 1 source changes committed for reference.

---

### Step 2b — ~~Roll out to remaining 5 subroutines~~ (OBSOLETE — superseded by config-level `fpe3` + `/assume:noieee_compares` in Step 2a)

<details><summary>Original per-routine plan (no longer needed)</summary>

Apply `!DEC$ OPTIONS /fpe=3` / `!DEC$ END OPTIONS` around each of these (re-grep line numbers first):

| Subroutine | Approximate start | Approximate end |
|---|---|---|
| `DynamicWaveXX_default_CG` | ~L11490 | ~L11738 |
| `ComputeFaceVelocityModulus` | ~L11373 | ~L11459 |
| `ModifyGeometryAndMapping` | ~L10802 | ~L11048 |
| `UpdateWaterLevels` | ~L13855 | ~L13941 |
| `ComputeNextDT` (**last — riskiest**) | ~L17305 | ~L17514 |

> `ComputeNextDT` must be validated individually — a NaN celerity silently corrupts the timestep.

</details>

---

### Optional experiment ✅ this is now the ADOPTED approach (see Step 2a)

Setting `FloatingPointExceptionHandling` from `fpe0` to `fpe3` at config level in all three vfproj files is exactly what was done in Step 2a — program-wide upside in one shot, no per-routine directives needed.

---

### Risk assessment

| Risk | Mitigation |
|---|---|
| Division by zero → silent Inf | Guarded by `Compute`/`BasinPoints`/`OpenPoints` masks in the loops |
| `sqrt` of negative depth → silent NaN | Step 1 `max(..., 0.0)` guards ✅ applied |
| CriticalFlow comparison imprecision | `CriticalFlow` is `real(8)` throughout — no precision issue |
| NaN propagation in `Me%lFlowX/Y` | Surfaces during HDF5 write; caught by `compare_mohid.py` |
| `ifx` directive syntax wrong | Step 2a build test catches this before rollout |

---

## In-Run A/B Profiling Harness (methodology — applies to Phase 3 onward)

**Problem:** Wall-clock and VTune CPU numbers drift between runs because of ambient load (Teams meetings, background processes, thermal throttling). Comparing a *baseline run* against a *separate optimized run* therefore mixes the optimization effect with machine noise.

**Solution:** Run the **baseline** and **performant** variants of the profiled routine(s) **inside the same process**, back-to-back, every timestep. Both variants then experience identical ambient conditions, so VTune's per-function CPU times give a clean, noise-cancelled delta in a *single* run.

### Protocol (per profiled routine)
1. **Snapshot** the persistent state the routine mutates (only needed for state-mutating routines).
2. Run the **baseline** variant.
3. Capture its outputs into `_baseline` variables.
4. **Restore** the snapshot so the performant variant starts from identical inputs.
5. Run the **performant** variant.
6. **Compare** performant vs baseline against the plan tolerances (abs `1e-5` **and** rel `1e-4` must *both* be exceeded to FAIL); **`stop` with a diagnostic** on mismatch.
7. **Overwrite** the outputs with the **baseline** values so the simulation trajectory always follows the trusted baseline (keeps multi-run behaviour bit-identical while validating the optimization live).

> For **read-only** routines (e.g. `ComputeNextDT`, which only returns a scalar `nextDTCourant`) steps 1/4/7 collapse to "just use the baseline scalar" — no snapshot/restore is needed.

### Runtime switch
A `runoff.dat` keyword controls the mode:

| `PROFILE_PERFORMANT_ONLY` | Behaviour |
|---|---|
| `0` (default) | **Dual-run A/B** — baseline + performant both run, compared, trajectory follows baseline. Use for profiling & live validation. |
| `1` | **Performant only** — no baseline, no compare, zero overhead. Use for the **final** HDF5/timeseries validation against the original. |

Implemented via `Me%ProfilePerformantOnly` (default `.false.`) in `ModuleRunOff.F90`. Tolerances are module parameters `ProfileTolAbs_ = 1.0e-5`, `ProfileTolRel_ = 1.0e-4`. Scalar comparison helper: `CheckProfileScalarDiff` (already in the code).

### Per-phase execution checklist (for the implementing agent)

> Follow these steps at the **start of every phase** that uses the A/B harness. This is the authoritative recipe — do not improvise.

1. **Branch:** cut `perf/PhaseN` from the *confirmed & cleaned* `perf/Phase(N-1)` (see Branch strategy).
2. **Classify the routine:**
   - **Read-only** (returns scalars only, mutates no persistent `Me%…` state) → use `CheckProfileScalarDiff`. No snapshot/restore. (Phase 3 pattern.)
   - **State-mutating** (writes persistent matrices, e.g. `Me%lFlowX/Y`) → you MUST add/reuse `CheckProfileMatrixDiff` (spec below) **at the start of this phase**, and implement the full snapshot → baseline → capture → restore → performant → compare → restore-to-baseline protocol.
3. **Duplicate the routine:** keep the current production code as `<Routine>_baseline`; put the optimized code in `<Routine>` (or `<Routine>_perf` extracted from a thin dispatcher). Copy the baseline **verbatim** so it is a true reference. **Then prove it is byte-for-byte identical to the previous phase's routine** — do NOT trust the copy:
   ```powershell
   git show perf/Phase(N-1):Software/MOHIDLand/ModuleRunOff.F90 > $env:TEMP\prev.F90
   # extract the routine body from prev.F90 and diff it against <Routine>_baseline;
   # the only allowed differences are the routine name and dummy-arg plumbing.
   ```
   If the diff shows any logic/expression/loop-bound change, the baseline is corrupt — fix it before proceeding.
4. **Wire the dispatcher** gated on `Me%ProfilePerformantOnly`:
   - `.true.` → call performant only.
   - `.false.` → run baseline + performant, compare via the appropriate `CheckProfile*Diff`, then follow the baseline result/state.
5. **Validate & profile** (see per-phase validation workflow).
6. **Clean up** once confirmed (see Cleanup rule), commit, then start Phase N+1.

### `CheckProfileMatrixDiff` specification (create at the start of the first state-mutating phase, then reuse)

Not yet in the code — add it in `ModuleRunOff.F90` next to `CheckProfileScalarDiff` when the first matrix-mutating phase begins. Required behaviour (mirror `CheckProfileScalarDiff` exactly so the two are interchangeable):

- **Signature:** `CheckProfileMatrixDiff(baselineMat, perfMat, routineName, Niter)` where `baselineMat`/`perfMat` are `real(8), dimension(:,:), pointer` (or matching the matrix type being compared — provide a `_R4` overload if a `real(4)` matrix is needed).
- **Scan** the computed cells only (respect `Me%WorkSize` bounds and the same `BasinPoints/OpenPoints == Compute` masks the routine uses — do NOT compare halo/inactive cells).
- **FAIL rule (identical to the scalar helper and `compare_mohid.py`):** for the worst cell, FAIL only if **both** `absDiff > ProfileTolAbs_` **and** `relDiff > ProfileTolRel_`, where `relDiff = absDiff / max(|baseline|, |perf|)` (guard `max(...) > AlmostZero`).
- **On FAIL:** `write(*,*)` the routine name, `Niter`, worst `(i,j)`, baseline value, perf value, abs & rel diffs and the tolerances, then `stop`.
- Reuse the module parameters `ProfileTolAbs_` / `ProfileTolRel_`.

### Snapshot / restore note (state-mutating phases only)

Allocate `<matrix>_snapshot` and `<matrix>_baseline` module or local buffers once (not every timestep — reuse them). Sequence per call: copy live matrix → `_snapshot`; run baseline; copy result → `_baseline`; copy `_snapshot` back to live; run performant; `CheckProfileMatrixDiff(_baseline, live, …)`; copy `_baseline` back to live so the trajectory follows the baseline.

### Cleanup rule
Once an optimization is **confirmed** (A/B passes tolerances and the performant-only run matches the original outputs), **delete the `_baseline` variant, the snapshot/restore buffers, and the dual-run branch of the dispatcher**, keeping only the performant routine. Then cut the next phase branch from that cleaned state.

### Branch strategy
**Continuous chain.** Each `perf/PhaseN` branch is cut from the *confirmed & cleaned* `perf/Phase(N-1)`. This measures each phase's incremental gain on top of prior confirmed gains, and the A/B baseline for phase N is exactly phase N-1's production routine.

---

## Phase 3 – `ComputeNextDT` Optimizations ✅ CONFIRMED & CLEANED

> **Status:** ✅ COMPLETE. Optimized, validated, baseline code deleted, committed on `perf/Phase3`. `ComputeNextDT` now calls `ComputeNextDT_CourantScan` directly with no dispatch overhead. `perf/Phase4` should be cut from this state.

**Location:** `ComputeNextDT` + `ComputeNextDT_CourantScan[_baseline]` (search these names in `ModuleRunOff.F90`)  
**Hotspot:** ~185s CPU NoRain 10T (completely unaffected by Phase 1); ~53s WithRain 10T. Scales ~2x (good parallelism but absolute cost still high).

### Optimizations applied (the performant variant):
1. **Hoisted `Distance_Courant`** out of the loop in the constant-grid path (computed once).
2. **`sqrt(max(aux, 0.0))`** guard (avoids the `for_is_nan` wrapper; also pairs with Phase 2 FPE relaxation).
3. **`sqrt_gravity = sqrt(Gravity)`** precomputed → `celerity = sqrt_gravity * sqrt(max(aux, 0.0))` (multiply + cheaper depth-only sqrt instead of `sqrt(Gravity*aux)`).
4. **Removed the `strideJ` reshape/transpose + inner `c` loop**; replaced with two explicit unrolled East/North face blocks (fewer integer ops, no per-cell 2×2 array indexing).

### Phase 3 A/B Results (Profile Double OpenMP — noise-cancelled in-run comparison)

| Scenario | Threads | `_baseline` CPU (s) | Performant CPU (s) | Gain |
|---|---|---|---|---|
| WithRain | 1T | 29.4 | 20.1 | **-32%** |
| WithRain | 10T | 32.6 | 20.3 | **-38%** |
| NoRain | 1T | 205.8 | 142.9 | **-31%** |
| NoRain | 10T | 166.6 | 126.1 | **-24%** |

> CPU times above are VTune effective CPU time (sum across all threads). For multi-threaded runs, wall-time gain on `ComputeNextDT` is approximately CPU_gain / 10. The NoRain case is by far the bigger beneficiary (baseline was ~7× higher than WithRain).

**A/B correctness:** `CheckProfileScalarDiff` did not trigger on either scenario — `nextDTCourant` agrees within tolerance every timestep.

### Remaining top hotspots after Phase 3 (for context — WithRain 10T)

| Function | CPU (s) | Notes |
|---|---|---|
| `DynamicWaveYY_default_CG` | 122 | Phase 4+ target |
| `DynamicWaveXX_default_CG` | 108 | Phase 4+ target |
| `_libm_pow_l9` | 85 | `**(1./3.)` friction — unchanged |
| `ComputeFaceVelocityModulus` | 65 | `abs(cmplx())` still in |
| `SetMatrixValues2D_R8_FromMatrix` | 58 | Phase 5 target |
| `ModifyGeometryAndMapping` | 44 | Phase 4 target |
| `OutputFloodingAll_R4` | 41 | Not previously identified — new candidate |
| `ComputeCenterVelocities_R4` | 40 | Not previously identified — new candidate |
| `ComputeNextDT_CourantScan_baseline` | 33 | Baseline only (will disappear after cleanup) |
| `UpdateWaterLevels` | 28 | Not previously identified |
| `_kmpc_barrier` (spin) | 25 | Load imbalance — Phase 6 target |
| `ComputeNextDT_CourantScan` | 20 | **Performant (Phase 3 result)** |

### Remaining top hotspots after Phase 3 (NoRain 10T)

| Function | CPU (s) | Notes |
|---|---|---|
| OMP spin `func@0x180099f40` | 330 | Thread spin/wait — severe load imbalance (Phase 6) |
| `ComputeNextDT_CourantScan_baseline` | 167 | Baseline only |
| `ComputeNextDT_CourantScan` | 126 | **Performant (Phase 3 result)** |
| `_kmpc_barrier` (spin) | 125 | Load imbalance — Phase 6 target |
| `SetMatrixValues2D_R8_FromMatrix` | 123 | Phase 5 target |
| `DynamicWaveXX_default_CG` | 120 | Phase 4+ target |
| `DynamicWaveYY_default_CG` | 118 | Phase 4+ target |
| `ModifyGeometryAndMapping` | 88 | Phase 4 target |
| `ComputeFaceVelocityModulus` | 84 | `abs(cmplx())` |
| `_libm_pow_l9` | 63 | Friction `**(1./3.)` |
| `ComputeCenterVelocities_R4` | 52 | New candidate |
| `OutputFloodingAll_R4` | 49 | New candidate |

> **Key new observations from Phase 3 profiles:**
> - `OutputFloodingAll_R4` and `ComputeCenterVelocities_R4` are now clearly visible hotspots (41s and 40s WithRain 10T; 49s and 52s NoRain 10T) — not previously identified. Worth investigating before Phase 4.
> - The NoRain OMP spin (`func@0x180099f40` = 330s spin) and `_kmpc_barrier` (125s) confirm severe load imbalance is the single largest NoRain cost — Phase 6 (scheduling) may be higher priority than Phase 4 for the NoRain scenario.
> - `SetMatrixValues2D_R8_FromMatrix` has jumped to 123s NoRain 10T (vs ~119s in original baseline) — ranking it higher than `DynamicWave` in NoRain.

### Validation workflow for this phase:
1. ~~Build `Profile Double OpenMP`, run with `PROFILE_PERFORMANT_ONLY : 0` → VTune A/B profiling~~ ✅ DONE (results above)
2. ~~Build `Release Double OpenMP`, run with `PROFILE_PERFORMANT_ONLY : 1` on both WithRain and NoRain, then `compare_mohid.py` against originals.~~ ✅ DONE — WithRain 36/37 PASS (pre-existing `FloodPeriod.dat` FAIL only), NoRain 38/38 PASS + direct HDF5 vs `_original` zero diff.
3. ~~Delete `ComputeNextDT_CourantScan_baseline` + the dual-run branch in `ComputeNextDT`; keep only `ComputeNextDT_CourantScan`. Commit. Cut `perf/Phase4` from here.~~ ✅ DONE

---

## Phase 4 – `ModifyGeometryAndMapping` Optimizations ✅ CONFIRMED & CLEANED

> **Status:** ✅ COMPLETE. Optimized, validated, baseline code deleted, committed on `perf/Phase4`. `perf/Phase5` should be cut from this state.

**Location:** `ModifyGeometryAndMapping` in `ModuleRunOff.F90`  
**Hotspot:** ~88s CPU NoRain 10T (Phase 3 baseline), ~44s WithRain 10T.

### Optimization applied:
The `strideJ` 2×2 identity array and the `do c = 1, size(strideJ,1)` dispatch loop unrolled into explicit **East face** (U, `j-1` neighbour) and **North face** (V, `i-1` neighbour) blocks in all three branches (`KinematicWave_`, `GridIsConstant`, variable-grid). Eliminates per-cell `strideJ` array indexing, the `c` loop counter, and the `if (dj==1)` branch from the hot inner loop. Zero arithmetic change — pure structural transformation.

### Phase 4 A/B Results (noise-cancelled in-run comparison, OMP parallel section times)

| Scenario | Threads | `_baseline` CPU (s) | Performant CPU (s) | Gain |
|---|---|---|---|---|
| WithRain | 1T | 32.9 | 27.7 | **−16%** |
| WithRain | 10T | 27.3 | 13.6 | **−50%** |
| NoRain | 1T | 98.1 | 79.1 | **−19%** |
| NoRain | 10T | 57.2 | 21.0 | **−63%** |

The 10T gains are exceptionally strong (-50%/-63%). The `strideJ` dispatch overhead consumed roughly half the routine's per-thread budget, suppressed at 10T by the compressed compute-to-overhead ratio.

**A/B correctness:** `CheckProfileMatrixDiff` never triggered on either scenario.

**Note on harness overhead (A/B profile runs only):** `CheckProfileMatrixDiff` cost ~100s in NoRain 10T (serial scan of all active cells every timestep — harness-only). The NoRain OMP spin also increased from ~330s to ~970s during the A/B profile run due to the serial 12-array snapshot/restore between OMP parallel calls; this overhead disappears entirely in the cleaned production code.

### Validation (Release Double OpenMP, PROFILE_PERFORMANT_ONLY:1):
- **WithRain:** 18/18 PASS vs `_original` (including `FloodPeriod.dat` — zero diff)
- **NoRain:** 18/18 PASS vs `_original` (zero diff on all HDF5 + text files)

### Remaining top hotspots after Phase 4 (WithRain 10T Profile run, `_perf` entries are production-representative):

| Function | CPU (s) | Notes |
|---|---|---|
| `DynamicWaveYY_default_CG` | 114 | Phase 5+ target |
| `DynamicWaveXX_default_CG` | 86 | Phase 5+ target |
| `_libm_pow_l9` | 72 | `**(1./3.)` friction — unchanged |
| `OutputFloodingAll_R4` | 41 | New candidate |
| `ComputeFaceVelocityModulus` | 38 | `abs(cmplx())` |
| `ComputeCenterVelocities_R4` | 33 | New candidate |
| `ModifyGeometryAndMapping` | ~14 | **Phase 4 result** |
| `SetMatrixValues2D_R8_FromMatrix` | 25 | |
| `ComputeNextDT_CourantScan` | 23 | Phase 3 result |

### Remaining top hotspots after Phase 4 (NoRain 1T, cleanest — no OMP spin distortion):

| Function | CPU (s) | Notes |
|---|---|---|
| `ComputeNextDT_CourantScan` | 163 | Phase 3 result |
| `DynamicWaveXX_default_CG` | 123 | |
| `DynamicWaveYY_default_CG` | 116 | |
| `SetMatrixValues2D_R8_FromMatrix` | 113 | |
| `ComputeFaceVelocityModulus` | 89 | |
| `_libm_pow_l9` | 79 | Friction `**(1./3.)` |
| `ModifyGeometryAndMapping` | **79** | **Phase 4 result** |
| `OutputFloodingAll_R4` | 44 | New candidate |
| `ComputeCenterVelocities_R4` | 37 | New candidate |

---

## Phase Reprioritization (post-Phase 4)

After Phase 4, the profile data confirms the two scenarios pull in different directions and the "easy, safe" wins are nearly exhausted. Code inspection of the remaining candidates:

| Candidate | Cost | Code reality | Risk | Owner |
|---|---|---|---|---|
| OMP load imbalance | NoRain ~455s (spin+barrier) | 95 `SCHEDULE(DYNAMIC,…)` sites; scheduling-only | Low correctness, broad + timing-measured | **Opus-led** (strategic) |
| DynamicWave XX/YY | ~230s WithRain / ~238s NoRain | friction `**(1./3.)` genuine cube-root; `real(8)` compares | **HIGH** — documented friction cell-flip regression | **Opus-led** (landmine) |
| `SetMatrixValues2D_R8_FromMatrix` | 113s NoRain | trivial masked copy; already `DYNAMIC`/`STATIC` split | Medium; generic util, many callers, other module | Opus or careful Sonnet |
| `ComputeFaceVelocityModulus` | 38–65s | `abs(cmplx())` still present | **Medium** — Phase 1 tried, REVERTED (noRain cell-flip) | Opus or skip |
| `ComputeCenterVelocities_R4` | 37–52s | 6× `**2.0` pow/cell in `sqrt(x**2.0+y**2.0)` | **Low** — output-only R4 arrays, Phase-1 pattern | ✅ **Sonnet** |
| `OutputFloodingAll_R4` | 41–49s | max/accumulate stats; no pow; memory-bound | Low but unclear gain | Investigate-only |

**Revised order** (the old Phase 5=SetMatrixValues / 6=load-balancing / 7=`_CG`-variants ordering is superseded; the `_CG` non-default variants aren't even on the hot path — the hot ones are `_default_CG`, already done in Phase 1). **Phase 6 (OMP load balancing) was investigated, implemented, measured, and concluded a DEAD END** — see the Phase 6 VERDICT; phases renumbered accordingly:

1. **Phase 5 → `ComputeCenterVelocities_R4`** — clean `**2.0`→`x*x`, output-only, low-risk. **Sonnet.** ✅ DONE
2. **Phase 7 → Output-only `**2.0`→`x*x` cleanup (`ComputeCenterValues`/`_R4`)** — same pattern as Phase 5, output-only, low-risk. **Sonnet.** ✅ DONE
3. **Phase 8 → `SetMatrixValues` / call-frequency reduction** — memory-bound; needs call-graph analysis. **Opus or careful Sonnet.**
4. **Phase 9 → DynamicWave XX/YY** — biggest WithRain compute, HIGH regression risk. **Opus, strict A/B, watch the friction cell-flip.**
5. **Phase 6 → OMP load balancing** — investigated, implemented & measured → ❌ **DEAD END** (near-no-op on the large domain; see Phase 6 VERDICT). **Opus-led.**

**Opus-only (do NOT hand to Sonnet):** DynamicWave friction and `ComputeFaceVelocityModulus abs(cmplx)` (both documented cell-flip landmines); OMP load-balancing strategy (95 sites, timing-only measurement, chunk/GUIDED judgment).

---

## Phase 5 – `ComputeCenterVelocities_R4` Optimizations ✅ CONFIRMED & COMMITTED

> **Status:** ✅ COMPLETE. Optimized, validated with zero diff on both scenarios, committed
> on `perf/Phase5`. `perf/Phase6` should be cut from this state.

**Location:** `subroutine ComputeCenterVelocities_R4` in `Software/MOHIDLand/ModuleRunOff.F90`  
**Hotspot:** ~37s WithRain 10T / ~52s NoRain 10T (newly visible after Phase 3).

### Optimization applied:
Six `**2.0` pow calls per cell — three occurrences of
`sqrt(Me%CenterVelocityX_R4(i,j)**2.0 + Me%CenterVelocityY_R4(i,j)**2.0)` (one each in
the `GridIsRotated`, non-rotated, and `Distortion` branches). These compile to
`_libm_pow_l9`, the exact hotspot class eliminated in Phase 1.

Replaced every `X**2.0` with `X*X`: introduced local `real(4) :: cx, cy` per branch, computed
once, then `sqrt(cx*cx + cy*cy)`; added `cx,cy` to the corresponding `!$OMP PARALLEL PRIVATE(...)`
clauses. Pure Phase-1 transformation — no changes to divisions, masks, or loop bounds. Did NOT
touch the identical `**2.0` pattern present in `UpdateFVSOutputVariables_CG_R4/_VG_R4` or
`OutputFloodingAll_R4`/`ComputeCenterValues_R4` (out of scope for this phase — candidate for a
future phase).

### Validation (compare_mohid.py, zero diff — not just within tolerance):
- **WithRain:** `Run55` dir vs `Run55_phase5_baseline` dir — 18/18 PASS, max abs/rel = 0.
  `RunOff_55.hdf5` vs `RunOff_55_phase5_baseline.hdf5` — 27 datasets, PASS, max abs/rel = 0.
- **NoRain:** `Run55` dir vs `Run55_phase5_baseline` dir — 18/18 PASS, max abs/rel = 0.
  `RunOff_55.hdf5` vs `RunOff_55_phase5_baseline.hdf5` — 57 datasets, PASS, max abs/rel = 0.

### Methodology note (output-only routine → no A/B harness):
The outputs (`CenterVelocityX_R4`, `CenterVelocityY_R4`, `VelocityModulus_R4`,
`CenterFlowX/Y_R4`) flow ONLY into output/statistics arrays (`Me%Output%…`, HDF5), never
back into the simulation trajectory (`lFlowX/lFlowY/iFlowX/iFlowY/myWaterColumn/AreaU/AreaV/DT`).
Because they cannot perturb the trajectory, the in-run A/B harness (dispatcher/`_baseline`/
`_perf`/`CheckProfileMatrixDiff`) is **not required** — validate directly via `compare_mohid`
vs `_original` on both scenarios. **Confirm output-only via a usage grep first**; if any
feedback path into the trajectory is found, STOP and use the harness instead (note: these are
`real(4)` arrays → would need a `CheckProfileMatrixDiff_R4` overload per that helper's spec).

### Scope guardrails — do NOT touch:
- `DynamicWaveXX/YY_default_CG` friction terms (`**(1./3.)`, real(8) compares) — regression landmine.
- `ComputeFaceVelocityModulus` `abs(cmplx())` — Phase 1 tried & REVERTED (noRain cell-flip).
- Any OMP `SCHEDULE` clauses (separate future phase).

---

## Phase 6 – OpenMP Load Balancing (❌ IMPLEMENTED & MEASURED → DEAD END — Opus-led)

> ### ❌ VERDICT — thread-count / barrier tuning is a DEAD END for this large domain
> Phase 6 was fully implemented (self-tuning `NUM_THREADS(RunOffBoxThreads(box))`, module param
> `MinCellsPerThread_ = 512`, on 9 `CurrentWorkSize` OMP regions) and measured with an in-run
> diagnostic counter + a temporary A/B harness + a `KMP_BLOCKTIME=0` run. **Conclusion: the gate
> is a near-no-op on the LargeModel and blocktime tuning yields nothing.** Two independent measurements
> agree:
>
> **1. Gate diagnostic (`RunOffBoxThreads` counter, printed at `KillRunOff`).**
> - **NoRain:** 204,478 calls, only **12.1%** (24,726) returned `< full team`; box `nCells`
>   min/avg/max = **20 / 445,965 / 1,540,472**. The gated calls all carry `< 5,120` cells, so
>   work-weighted the gate touches **≈ 0.14 % of the cell-work** (max `24,726 × 5,119 ≈ 1.27e8`
>   vs total `204,478 × 445,965 ≈ 9.12e10`). ~99.9 % of compute still runs the full 10-thread team.
> - **WithRain:** **0.00%** gated, box constant `1,540,472`, always 10 threads → **provable no-op,
>   zero regression** (full domain every step).
> - Corroborates the earlier VTune read (barrier wait counts down only ~9–18 %).
>
> **2. `KMP_BLOCKTIME=0` run (NoRain, `LOGKMPSet.dat` vs `LOG_noKMPSet.dat`) — CONFOUNDED / no benefit.**
> - The KMP=0 run is uniformly ~7–33 % slower **including serial startup routines `KMP_BLOCKTIME`
>   cannot physically affect** (`ConstructBasin` +23 %, `ConstructRunOff` +27 %, `RunOffOutput` +33 %)
>   ⇒ the delta is run-to-run ambient/thermal noise, **not** blocktime (two separate whole-app runs
>   do **not** cancel noise — that is exactly why the in-process A/B harness existed).
> - The heaviest, region-dense parallel routines were inflated **least** (`DynamicWaveYY` +4.9 %,
>   `DynamicWaveXX` +7.2 %), not most ⇒ **no wake-up-cost signal** either.
> - `ModifyRunOff` wall **rose** (147.9 → 164.2 s) instead of dropping ⇒ the ~692 s busy-spin is
>   **not reclaimable critical-path slack**; it is threads kept hot between closely-spaced regions.
>
> **Decisions:** (a) Keep the `NUM_THREADS` gate only as a correctness-safe, zero-cost generic
> guardrail — it is **not** the Phase 6 win. (b) **Reject** raising `MinCellsPerThread_`: to bite
> the bulk work it would need ≈ 44,000 (= 445,965 / 10), which under-threads medium/large boxes that
> *do* benefit from parallelism → net loss + brittle (WithRain would not regress until > 154,047,
> but NoRain loses first). (c) **STOP thread tuning; pivot to reducing WORK** — see "Pivot" below.
> All temporary harness code (`_AB` wrappers, `RunOffForceThreads`, diagnostic vars/print, snapshot
> buffers) is to be removed before any commit; the ERR010 blow-up it caused is fixed (see Resolved Bug).
>
> ### Pivot — the real NoRain levers (clean hot list from the less-noisy default run)
> Per-routine **Wall** inside `ModifyRunOff` (147.9 s), with `CPU/Wall` ≈ threads used:
>
> | Routine | Wall (s) | ≈threads/10 |
> |---|---|---|
> | `ComputeNextDT` | 22.2 | 8.8× |
> | `DynamicWaveXX_default_CG` | 20.4 | 9.0× |
> | `DynamicWaveYY_default_CG` | 19.2 | 8.8× |
> | `ModifyGeometryAndMapping` | 13.5 | 7.4× |
> | `ComputeFaceVelocityModulus` | 12.3 | 8.5× |
> | `SetFlowOldXY` | 10.9 | 9.2× |
> | `OutputFloodingAll_R4` | 8.5 | 7.1× |
> | `SetWorkSize` | 7.9 | 7.3× |
> | `UpdateWaterLevels` | 7.5 | 7.1× |
>
> These already run at **7–9× of 10 threads** → parallelism is *already good*; further gains must
> **reduce work, not add threading**. Top levers: **`ComputeNextDT` (22.2 s) + `ModifyGeometryAndMapping`
> (13.5 s)** sweep the full 1.54M-cell domain while the wet box is tiny → restrict to the active
> bounding box (A/B-guarded; Phase-1 landmine: `ComputeNextDT` is a min-reduction, `ModifyGeometryAndMapping`
> mutates geometry). Then `DynamicWaveXX/YY` (~40 s combined) = per-cell arithmetic (Phase 5/7/11
> bit-identical style).

**Observation:** NoRain OMP spin (`func@0x180099f40`) ~330s + `_kmpc_barrier` ~125s = ~455s
of load imbalance in production 10T (unaffected by Phase 1–4 code changes since it is a
scheduling issue). This is the single largest NoRain cost and dwarfs every individual compute
hotspot. Minor for WithRain (~25s barrier).

**Cause:** Uneven distribution of "wet" cells across rows.

### Current scheduling (95 `!$OMP DO SCHEDULE(DYNAMIC, ChunkJ/CHUNK)` sites):
```fortran
!$OMP DO SCHEDULE(DYNAMIC, CHUNKJ)
```

### Options:
1. ~~Tune `CHUNKJ` chunk size~~ — investigated & REJECTED as the NoRain lever (Experiment 0: dispatch
   is ~0.5s, coarsening made the barrier worse). May still help large WithRain domains — revisit.
2. ~~Use `SCHEDULE(GUIDED)`~~ — REJECTED (user: countless wet/dry variations → GUIDED unreliable)
3. **Scale threads-per-region to the active-box size (self-tuning, no keyword)** ← the direction for
   the future Phase 6 (see "Direction" below).

> **Why Opus-led:** broad (95 sites), correctness-safe (scheduling doesn't change math →
> `compare_mohid` trivially passes) but the *effect* is timing-only, so the in-run A/B harness
> does not apply — needs careful wall-time measurement methodology and chunk/GUIDED judgment.

### Survey findings (branch `perf/Phase6`, cut from `perf/Phase5` a620bad5)

- **All 95 `SCHEDULE(DYNAMIC, …)` sites in ModuleRunOff collapse to one value.** 21 sites use
  `CHUNKJ` (the global directly), 74 use a local `CHUNK` set via `CHUNK = ChunkJ`. Fortran is
  case-insensitive so `CHUNKJ` == `ChunkJ`.
- **`ChunkJ` is a GLOBAL** (`ModuleGlobalData.F90` L2181, default `1`), computed once at
  construction in `ModuleBasin.F90` L897: `ChunkJ = max((JUB-JLB)/ChunkJFactor, 1)` with
  `ChunkJFactor` default `99999` ⇒ **`ChunkJ = 1`** for any normal grid ⇒ every loop runs
  `SCHEDULE(DYNAMIC, 1)` (finest granularity: best balance, worst dispatch-atomic overhead).
- **`ChunkJ` is shared** with PorousMedia, RunoffProperties, MOHIDBase, and all of MOHIDWater →
  must NOT change its computation. Solution: ModuleRunOff-local chunk (decoupled).
- **Loops parallelize over columns `j`** (outer), rows `i` (inner). Physics loops iterate
  `Me%CurrentWorkSize`; **`ComputeNextDT_CourantScan` and `ModifyGeometryAndMapping` iterate the
  full `Me%WorkSize`.**
- **`Me%CurrentWorkSize` is a dynamic bounding box grown around active cells ONLY when
  `.not. Me%HasRainFall`** (L10616-10653, `SetWorkSize` L11097). WithRain uses the full domain.
  **This is exactly why NoRain is imbalanced (small irregular box over dynamic,1) and WithRain
  is not.**
- `openmp_num_threads` (global, `ModuleGlobalData`) holds the real thread count → usable in a
  chunk formula.

### Root cause — REVISED after Experiment 0 (the chunk hypothesis was WRONG)

**Experiment 0 (committed Phase 5 baseline + `CHUNK_I_FACTOR : 40`, NoRain 10T, my Phase 6 code
NOT built)** disproved the dispatch-contention hypothesis:

| libiomp function | CPU | role |
|---|---|---|
| `_kmpc_dispatch_next_4` | **0.541s** | DYNAMIC chunk hand-out |
| `_kmpc_dispatch_init_4` | 0.125s | loop dispatch setup |
| `func@0x180099f40` | 326s (**315.7s spin**, wait cnt **4,730,938**) | idle worker spin at fork/join |
| `_kmpc_barrier` | 243s (**243.0s spin**, wait cnt **4,850,561**) | region-exit barrier |
| `func@0x180047260` | 107s (**102.5s spin**) | join/reduction wait |

- **Dynamic scheduling costs ~0.67s total** → chunk size is NOT the bottleneck; there is no
  dispatch-atomic contention. The premise behind `RunOffChunk` was wrong for NoRain.
- The real cost is **~692s of threads idle-spinning at fork/join barriers** vs ~656s of useful
  compute — **~half the CPU wasted**. It is structural: **4.85M barrier calls** = a huge number of
  *tiny* parallel regions over the small NoRain active box.
- **Coarsening made it worse:** `CHUNK_I_FACTOR : 40` sets `ChunkJ = fullJ/40` (a large fixed chunk
  from the *full* domain); on the small NoRain box that starves threads → `_kmpc_barrier` ~doubled
  (243s vs ~125s baseline). Confirms the active box is small and the lever is *parallelism*, not
  *chunk*.

**Implication:** for a small box, `RunOffChunk` collapses to 1 anyway (= baseline), and no chunk
value can reduce the per-region barrier count. Phase 6 pivots to **not going parallel when the box
is too small**.

### Status — HISTORY (superseded by the VERDICT block at the top of this section)

> The notes below are the *earlier* deferral rationale, kept for context. Phase 6 has since been
> implemented and measured to a dead end — read the VERDICT block above for the current status.

Phase 6 work is **reverted** (`ModuleRunOff.F90` back to clean Phase 5 `a620bad5`) and moved to be
the **last** phase. Do Phase 7 (and any others) first. Reasons (user):
1. A tuned cell-count threshold is **model/grid-specific** → brittle across other grids/setups.
2. Adding a `runoff.dat` keyword (`OMP_PARALLEL_MIN_CELLS` / `OMP_CHUNK_PER_THREAD`) forces a
   **downstream UI change** (the GUI that writes the data files) — not worth it unless gains are large.
3. Needs careful review; revisit when there is time.

**What we keep from this investigation (valuable, don't lose):** the Survey findings and the
**Experiment 0 root cause** above — dispatch is ~0.5s, so **chunk size is not the lever**; the cost
is ~692s of fork/join **barrier idle-spin** over ~4.85M tiny parallel regions on the small NoRain
active box.

### Direction for the (future) real Phase 6 — self-tuning, NO keyword

Instead of an `IF()` on/off gate with a magic threshold, **scale the thread count per region to the
work size** using a grid-independent internal constant (a module PARAMETER, *not* a data-file
keyword — so no UI change):

```fortran
! module parameter, grid-independent (tune once in code):
integer, parameter :: MinCellsPerThread_ = <e.g. 512>

! before each CurrentWorkSize-based parallel region:
nCells  = (Me%CurrentWorkSize%IUB-Me%CurrentWorkSize%ILB+1) &
        * (Me%CurrentWorkSize%JUB-Me%CurrentWorkSize%JLB+1)
nThreadsBox = max(1, min(openmp_num_threads, nCells / MinCellsPerThread_))
!$OMP PARALLEL NUM_THREADS(nThreadsBox) ...
```

- Tiny box → `NUM_THREADS(1)` = serial, no fork/join/barrier (kills the NoRain overhead).
- Large box (WithRain full domain, or NoRain once water spreads) → full threads, as today.
- Self-scaling to *any* grid and thread count; the only constant is "min useful cells per thread",
  which is physically about amortizing fork/join overhead, not grid-specific.
- Still parallelism-only → `compare_mohid` must be zero-diff (same iterations; REDUCTION valid at
  any thread count). Validate on BOTH scenarios; verify no WithRain regression.
- Apply to the `CurrentWorkSize` box loops (`ComputeFaceVelocityModulus`, `DynamicWaveXX/YY_default_CG`,
  `OutputFloodingAll_R4`, and candidates `UpdateWaterLevels`, `ComputeCenterVelocities_R4`,
  `CalculateTotalStoredVolume`). Leave the full-`WorkSize` loops (`ComputeNextDT_CourantScan`,
  `ModifyGeometryAndMapping`) always-parallel.
- Open question to weigh at review time: `omp_get_max_threads`/nested-thread-count interactions,
  and whether `NUM_THREADS` per region has measurable call overhead vs the barrier it removes.

> The chunk idea (Lever B, `RunOffChunk`) *may* still help large WithRain domains, but per Experiment 0
> it is not the NoRain lever and is **not** worth a keyword on its own. Fold it in only if the
> NUM_THREADS work shows a clear, keyword-free way to set it (or drop it).

### Measurement plan (EXECUTED — see the VERDICT block at the top of Phase 6)

- ~~No A/B harness (parallelism-only → math unchanged)~~ — a **temporary** in-run A/B harness +
  diagnostic counter *was* added to measure the gate noise-cancelled; it has since been removed
  (it caused the `CheckStability ERR010` blow-up, now fixed — see Resolved Bug). Result: the gate
  touches ≈ 0.14 % of cell-work (NoRain) / 0 % (WithRain).
- `compare_mohid.py` was zero-diff on BOTH scenarios (parallelism-only, as expected).
- `KMP_BLOCKTIME=0` NoRain run done → **inconclusive/no benefit** (confounded by run-to-run noise;
  serial startup inflated more than the parallel regions). Busy-spin is not reclaimable critical-path
  slack. **Net: thread/blocktime tuning is a dead end; pivot to reducing work.**

---

## Phase 7 – Output-only `**2.0`→`x*x` cleanup (`ComputeCenterValues` / `ComputeCenterValues_R4`) ✅ CONFIRMED

> **Status:** ✅ COMPLETE. Extends the Phase 5 `ComputeCenterVelocities_R4` cleanup to the two
> remaining output-only routines that shared the identical `sqrt(X**2.0 + Y**2.0)` pattern.
> Committed on `perf/Phase7_OutputPow` (cut from `perf/Phase6` ef643b84, which itself is
> Phase 5 code + the Phase 6 deferral doc only — no source change).

**Location:** `subroutine ComputeCenterValues` and `subroutine ComputeCenterValues_R4` in
`Software/MOHIDLand/ModuleRunOff.F90`.

### Scope decision — FVS output routines explicitly excluded
The same `**2.`/`**2.0` pattern also exists in `UpdateFVSOutputVariables_CG`,
`UpdateFVSOutputVariables_CG_R4`, `UpdateFVSOutputVariables_VG`, `UpdateFVSOutputVariables_VG_R4`
(the FVS-solver output-statistics routines). **User explicitly excluded these from Phase 7** — no
 changes were made there. They remain a candidate for a possible future phase if ever revisited.

### Full survey of `**2.`/`**2.0` sites in `ModuleRunOff.F90` (45 sites found)

**SKIP (21 sites) — trajectory feedback confirmed, explicit guardrail, or dead code:**

| Lines | Routine | Reason |
|---|---|---|
| 9175, 9201 | `ComputeStateFVS_CG` | `velMod` → friction (`tau_u/tau_v`) → `velocityU/V` → `Me%VelModFaceU/V`, `Me%myWaterColumn`. Trajectory-critical. |
| 9353, 9377 | `ComputeStateFVS_VG` | Same pattern. Trajectory-critical. |
| 9370 | `ComputeStateFVS_VG` | Commented-out line — dead code. |
| 9594–9616, 9763–9785, 9933–9955, 10104–10126 | `UpdateFVSOutputVariables_CG/_CG_R4/_VG/_VG_R4` | Output-only (same proof as below) but **excluded per user instruction — no FVS routine changes this phase.** |
| 11584, 11839, 12102, 12436 | `DynamicWaveXX_*` | `Me%OverlandCoefficientX**2.` friction term — guardrail (`OverlandCoefficient*`/`Friction`). |
| 12789, 13049, 13330, 13648 | `DynamicWaveYY_*` | `Me%OverlandCoefficientY**2.` friction term — guardrail. |
| 13323, 13327, 13641, 13645 | `DynamicWaveYY_*` | Commented-out lines — dead code. |
| 15337, 15338 | `FlowFromChannels` | Quadratic solve for `NewLevel` → `dVol` → `Flow` → channel/runoff volume exchange (`Me%myWaterVolume`/`myWaterColumn`). Trajectory-critical — found via full survey, not on original candidate list. |
| 17450, 17497 | `ComputeNextDT_CourantScan` | `Distance_Courant` → `nextDTCourant` → `DT`. Guardrail (`DT`). |

**OUTPUT-ONLY, FIXED (8 sites, 2 routines):**

| Routine | Lines |
|---|---|
| `ComputeCenterValues` | 16933–16934, 16960–16961 |
| `ComputeCenterValues_R4` | 17151–17152, 17173–17174 |

**Usage-grep evidence (output-only proof):** `Me%FlowModulus`/`Me%FlowModulus_R4` (42 refs) and
`Me%VelocityModulus`/`Me%VelocityModulus_R4` (49 refs) flow only into `Me%Output%*`
(`MaxFlowModulus`, `VelocityAtMaxWaterColumn`, `FloodRisk`, `WeightedVelocity`) and HDF5/`Array2D`/
`Data2D` writers — never into `lFlowX/Y`, `iFlowX/Y`, `myWaterColumn/Volume`, `AreaU/V`, `DT`,
`VelModFaceU/V`, `OverlandCoefficient*`, or `Friction`. `Me%CenterFlowX/Y[...]` and
`Me%CenterVelocityX/Y[...]` (the sqrt operands) are themselves written only from read-only
consumption of `Me%VelModFaceU/V` or `iFlowX/Y` earlier in the same routine — never written back
into trajectory arrays. Same conclusion as Phase 5's already-validated `ComputeCenterVelocities_R4`.

`OutputFloodingAll`/`OutputFloodingAll_R4` (originally flagged as candidates) contain **no**
`**2.` patterns at all — nothing to change there.

### Optimization applied
For each of the 8 sites: introduced local scalars — `real :: cfx, cfy, cvx, cvy` in
`ComputeCenterValues` (double, matches `Me%CenterFlowX`/`CenterVelocityX`), `real(4) :: cfx, cfy,
cvx, cvy` in `ComputeCenterValues_R4` (matches `Me%CenterFlowX_R4`/`CenterVelocityX_R4`) — assigned
from `Me%CenterFlowX/Y[...]`/`Me%CenterVelocityX/Y[...]` once per cell, then
`sqrt(cfx*cfx + cfy*cfy)` / `sqrt(cvx*cvx + cvy*cvy)` replacing `sqrt(X**2. + Y**2.)`. Added
`cfx,cfy,cvx,cvy` to the enclosing `!$OMP PARALLEL PRIVATE(...)` clauses (both
`WriteMaxFlowModulus` true/false branches in each routine). Pure Phase-1/Phase-5 transformation —
bit-identical, zero diff expected.

### Validation (correctness)
No A/B *correctness* harness needed (output-only, same methodology as Phase 5) — validate directly
with `compare_mohid.py` on WithRain + NoRain (HDF5 at `res/` root as `RunOff_55*.hdf5`); zero-diff
expected on both scenarios. Known accepted FAIL: `FloodPeriod.dat`.

### Performance-measurement attempt (in-run A/B timing) — TRIED, then REMOVED
The `x**2.`→`x*x` change eliminates `pow` calls, but comparing VTune CPU times from two **separate**
builds/runs reintroduces the ambient-noise problem the Phase 3+ in-run A/B harness solves — and the
real delta on these small routines is smaller than run-to-run noise. So a timing-only in-run dual
computation was added (baseline `x**2.` block under its own `StartWatch` label + performant `x*x`
block, gated by `PROFILE_PERFORMANT_ONLY`, with `REDUCTION` sums fed through `CheckProfileScalarDiff`
to block dead-code elimination and as a bonus correctness check).

**Two findings killed it:**

1. **The routines are cold paths — measurement is below the noise floor.** VTune 1T (Profile build):
   `ComputeCenterValues_R4` ≈ **0.14s** in WithRain (three sub-regions 0.095 + 0.030 + 0.015s),
   **absent** in NoRain; `ComputeCenterValues` (double) **never runs** (model uses single-precision
   output). Reason: `ComputeCenterValues_R4` runs its Modulus block **only on output timesteps**
   (`ComputeEverything = WriteHdf .or. WriteTimeSerie .or. WriteMaxFlowModulus`); every other step it
   delegates to `ComputeCenterVelocities_R4` — which is the real per-timestep hotspot (**12.4s
   WithRain 1T**) and was **already optimized in Phase 5**. So the Phase 7 pow-elimination is correct
   and harmless but its performance benefit is ≈0.

2. **The harness had a false-positive `stop` bug (`_R4` routine only).** The baseline sum-loop reads
   the *existing* `Me%FlowModulus_R4`/`VelocityModulus_R4` for inactive cells
   (`OpenPoints /= BasinPoint`), while the performant loop runs the original conditional zeroing
   (`if (… /= 0.0) … = 0.0`) on those same cells before summing. Any cell that **dried out** (was
   active, now inactive) still holding a stale nonzero modulus → baseline sum includes the stale
   value `V`, performant sum includes `0` → the sums diverge → `CheckProfileScalarDiff` exceeds
   tolerance and `stop`s the run at the first output step with a newly-dry cell. This is a
   **harness measurement artifact, not a defect in the `x*x` change** (per-cell active values are
   bit-identical; `compare_mohid.py` and the performant-only run pass cleanly). The double-precision
   `ComputeCenterValues` did not have this bug (it unconditionally writes `0.0` to inactive cells, so
   baseline and performant agree). A correct in-place fix would be to accumulate the sums **only
   inside the active branch** in both loops — but given finding (1) it wasn't worth keeping.

**Resolution:** the entire timing harness was removed (source restored to the clean `x*x`-only
Phase 7 state, commit `c52de3f8`). The committed Phase 7 code is purely the bit-identical
`x**2.`→`x*x` transformation on the 8 `ComputeCenterValues`/`_R4` sites — validated by
`compare_mohid.py` zero-diff. Lesson recorded: for **cold** routines, skip in-run A/B timing
entirely — the change is justified as a bit-identical Phase-1-class cleanup, not by a measured delta.

---

## Phase 8 – `OutputFloodingAll_R4` cache/hoist cleanup ✅ CONFIRMED & COMMITTED

> **Status:** ✅ COMPLETE. Committed on `perf/Phase8_FloodStats` (commit `323430af`, cut from
> `perf/Phase7_OutputPow` `8d324ae7`). Chosen instead of Phase 9/`SetMatrixValues` as a lighter,
> lower-risk output-only candidate identified during the Phase 7 hotspot survey
> (`OutputFloodingAll_R4` ≈ 22.5s WithRain 1T, no `pow`, memory/stats-bound).

**Location:** `subroutine OutputFloodingAll_R4` in `Software/MOHIDLand/ModuleRunOff.F90`
(~L18694), called every timestep from `Outputs` → `ModifyRunOff`.

### Output-only proof
Usage-grepped all 9 variables the routine writes (`Me%Output%MaxWaterColumn_R4`,
`VelocityAtMaxWaterColumn_R4`, `TimeOfMaxWaterColumn`, `MaxFloodRisk_R4`, `FloodPeriods`,
`FloodArrivalTime`, `TotalFloodedArea`, `MaxTotalFloodedArea`, `TimeOfMaxTotalFloodedArea`)
across all of `Software/**/*.F90`. Every reference is confined to `ModuleRunOff.F90`, and every
sink is `WriteGridData`/`GridData2D_Real`, `maxval`/`maxloc` (log only), `deallocate`, or
`Me%Output%FloodPeriod(i,j)` (itself only feeding `WriteGridData`). Zero references anywhere near
`lFlowX/Y`, `iFlowX/Y`, `myWaterColumn/Volume`, `AreaU/V`, `DT`, `VelModFaceU/V`,
`OverlandCoefficient*`, `Friction` — confirmed output-only, no A/B harness needed (same
methodology as Phase 5/7).

### Optimization applied
Most loop-invariant hoisting (`FloodWaterColumnLimit`, `NFloodPeriodLimits`, the `ComputePoints`
pointer selection, caching `Me%myWaterColumn(i,j)`) **already existed** in this routine before
Phase 8 — not a fresh win. Two remaining redundant-read eliminations applied:

1. **Cache `Me%VelocityModulus_R4(i,j)`** into a new local `VelMod` (`real(4)`, added to the
   `!$OMP PRIVATE` clause) — was read twice per active cell (once in the `MaxWaterColumn` branch,
   once in the `FloodRisk` calc); now read once.
2. **Hoisted `Me%Output%FloodRiskVelCoef` and `Me%Output%FloodArrivalWaterColumnLimit`** (both
   read-only scalars, never mutated in this routine) into locals `FloodRiskVelCoefLocal` /
   `FloodArrivalWaterColumnLimitLocal`, computed once before the parallel region.

Both new locals preserve the exact original declared types, so arithmetic promotion order is
unchanged — the change is a pure redundant-load elimination, not an arithmetic reformulation.
Diff size: +10/-5 lines. Left untouched: the inner `do n = 1, NFloodPeriodLimits` `FloodPeriods`
indexing (a local-array copy would add per-timestep allocation overhead — net negative), the
`Me%GridIsConstant` branch (cheap/predictable, not worth duplicating the loop), the
drainage-network `WorkSize` loop below (outside the hot per-timestep `CurrentWorkSize` region),
and the non-`_R4` `OutputFloodingAll` (guardrail — not mirrored).

### Validation
`compare_mohid.py`:
- **WithRain:** 36/37 PASS. Only fail is the pre-existing, accepted `FloodPeriod.dat` timing-shift
  vs the Phase 1 baseline (unrelated to Phase 8). HDF5 (27 datasets) PASS.
- **NoRain:** **38/38 PASS**, zero fails (57 HDF5 datasets).
- Explicit check vs `RunOff_55_phase5_baseline.hdf5` (pre-Phase 7/8, both of which are pure
  redundant-read eliminations expected to be bit-identical): PASS on both scenarios, max diff
  ≈1.19e-07/2.38e-07 — the float32 machine-epsilon noise floor ($2^{-23}$), identical in
  magnitude to the pre-existing diff vs `_original`, i.e. accumulated single-ULP rounding from
  the already-confirmed Phases 1–5, not a new regression from Phase 8.

### Performance
VTune single-run comparisons (separate builds, different days) showed `OutputFloodingAll_R4`
at ~30s WithRain 1T / ~41s NoRain 1T vs. the ~22.5s figure noted from the Phase 7 survey — **not
a reliable signal** (same ambient-noise caveat as the Phase 7 lesson; this routine is too cheap
for cross-run VTune deltas to mean anything, and no in-run A/B timing harness was built for the
same cold/hot-path-cost-benefit reasons as Phase 7). The change is justified as a correct,
bit-identical micro-optimization (redundant-load elimination), not by a measured delta.

---

## Phase 9 / Phase 11 – `DynamicWaveXX_default_CG` / `DynamicWaveYY_default_CG` deeper optimization (PARTIAL — coeff² precompute done, cube-root deferred)

> **Status:** ✅ Safe half done on `perf/Phase11_DynamicWave` (cut from `perf/Phase10_SetMatrixCallFreq`
> `16a37dcf`). Cube-root reformulation intentionally deferred (Phase-1 landmine).

The biggest single compute hotspot overall (~230s WithRain, ~238s NoRain combined). The
dominant remaining cost is `_libm_pow_l9` from the friction term `**(1.0/3.0)` (a genuine
cube-root) plus the physics arithmetic.

### Key insight (Phase 11)
The friction term has **two** pow calls per cell per timestep:
`VelModFace * OverlandCoefficient(i,j)**2. / (HydraulicRadius**(4./3.))`.
`ConstructOverLandCoefficient` runs **once** at init (`ConstructRunOff` L1084) → `OverlandCoefficientX/Y`
are **constant for the whole run**, so `coeff**2.` is a per-cell constant recomputed every timestep
for nothing.

### Applied (bit-identical, zero-risk — no A/B harness needed)
Precompute `OverlandCoefficientX/Y**2.` **once** into new fields `OverLandCoefficientXSquare`/`YSquare`
(filled inside `ConstructOverLandCoefficient`'s guarded loops using the **identical** `** 2.` expression,
so the cached bits equal exactly what the hot loop produced), then read the cached array in the friction
line of both `_default_CG` routines (`DynamicWaveXX` L11644, `DynamicWaveYY` L12849). Removes one
`_libm_pow_l9` call per cell per timestep with **zero** result change → `compare_mohid.py` zero-diff
expected (same class as Phase 5/7/8 bit-identical cleanups). Cold guardrail variants (`_default_VG`,
`_CG`, `_VG`) intentionally left untouched (only `_default_CG` runs in this model).

### Deferred (HIGH REGRESSION RISK — Phase-1 landmine)
The genuine cube-root `HydraulicRadius**(4./3.)` (varies per timestep) has **no bit-identical**
reformulation. Phase 1 attempted friction reformulations here (`coeff**2.`→`coeff*coeff`,
`HydraulicRadius**(4./3.)`→`HydraulicRadius*HydraulicRadius**(1./3.)`) and had to REVERT them —
they tipped a borderline `OpenPoints` cell across the active threshold in the noRain case
(see Phase 1 / Phase 2 Step 1 notes). Any future attempt must go through the full in-run A/B
harness (`CheckProfileMatrixDiff` on `lFlowX`/`lFlowY`) and be validated on BOTH scenarios with
extreme care. Opus-led. **Not attempted in Phase 11 by user decision.**

---

## Phase 10 – `SetMatrixValues2D_R8_FromMatrix` / call-frequency reduction ✅ CONFIRMED & COMMITTED

> **Status:** ✅ COMPLETE. Committed on `perf/Phase10_SetMatrixCallFreq`, cut from
> `perf/Phase8_FloodStats` (`323430af`).

**Location:** `Software/MOHIDLand/ModuleRunOff.F90`, `ModifyRunOff`'s per-sub-iteration hot path.
**Hotspot:** ~113s NoRain 1T / ~123s NoRain 10T for the shared
`SetMatrixValues2D_R8_FromMatrix` copy routine (`Software/MOHIDBase1/ModuleFunctions.F90`
~L1300), called from ~80+ sites across the codebase.

### Call-graph analysis
`SetMatrixValues2D_R8_FromMatrix` itself is already a trivial masked/unmasked copy with a
`DYNAMIC`/`STATIC` schedule split — memory-bandwidth bound, so no code-level upside in the copy
itself. Because it has ~80 unrelated callers, the routine can't be changed in isolation; the
lever is **reducing how often `ModifyRunOff` calls it**. Traced every `SetMatrixValue(...)` call
site inside `ModifyRunOff`: no site could be safely deleted (each either feeds
`DynamicWaveXX/YY` friction/advection — trajectory-critical — or is stability-critical), but 3
sites are **paired X/Y calls** (2 separate calls, each its own OMP fork/join, copying twinned
arrays under the same mask) that can be merged into 1 call with a single OMP region:
`FlowXOld`/`FlowYOld` (every sub-iteration — highest frequency), `InitialFlowX`/`InitialFlowY`
(once per `ModifyRunOff` call), and `lFlowX`/`lFlowY` restore-on-restart (once per restart
retry).

### Optimization applied
Merged each X/Y pair into a single subroutine (`SetFlowOldXY`, `SetInitialFlowXY`) with one
`!$OMP PARALLEL`/`DO` region copying both arrays under the shared mask check in one pass, halving
the fork/join count for that call site. Measured via an in-run dual-timing block (both the
original 2-call form and the merged form run back-to-back under separate `StartWatch`/`StopWatch`
labels, gated by a temporary keyword) so VTune/`ModuleStopWatch` could isolate the cost delta
from the ~80 unrelated callers sharing the same compiled OMP-region symbol — the harness and its
keyword were removed once the gain was confirmed (see cleanup below).

**`lFlowX`/`lFlowY` restore-on-restart (`SetLFlowRestoreXY`) was implemented, then REVERTED**:
VTune showed it never appeared in either scenario's hotspot report (the `Niter > 1` restart-retry
branch is rare/negligible in this benchmark), so the added complexity wasn't worth it. Final
kept scope = `SetFlowOldXY` + `SetInitialFlowXY` only.

### Performance (`ModuleStopWatch`, clean per-label CPU/Wall, no cross-caller contamination)
| Routine | NoRain baseline → merged | WithRain baseline → merged |
|---|---|---|
| `SetInitialFlowXY` | 27.9/29.5s → 19.3/19.4s (~31% CPU cut) | 9.8/10.1s → 8.6/8.8s (~12% cut) |
| `SetFlowOldXY` | 52.2/52.5s → 39.8/42.0s (~24% CPU cut) | 10.3/10.5s → 9.3/9.2s (~10% cut) |

Both merges show consistent, real wins (baseline always strictly slower, all 4 comparisons).
NoRain benefits more in absolute seconds (more `doIter`/restart-loop churn than WithRain in this
benchmark).

### Validation
`compare_mohid.py`: **WithRain 36/37 PASS** (only the pre-existing, accepted `FloodPeriod.dat`
timing-shift fail vs the Phase 1 baseline, unrelated to Phase 10; `RunOff_55.hdf5` PASS, max
abs/rel 1.192e-07). **NoRain 38/38 PASS**, zero fails (all `.srr`/`.dat` + both `RunOff_55.hdf5`
comparisons vs `_original` and `_phase1`, max diffs ~1e-7 float noise).

### Cleanup
Dual-timing harness (`_Baseline`/`_Performant` variants, dispatcher subroutines, the
`ProfilePerformantOnly` field and its `PROFILE_PERFORMANT_ONLY` keyword) removed once the gain
was confirmed — only the merged `SetFlowOldXY`/`SetInitialFlowXY` subroutines remain, called
directly from `ModifyRunOff` under their original names.

---

## Superseded (kept for history) – `DynamicWaveXX_CG` / `DynamicWaveYY_CG` non-default variants

The `_CG` (non-`default`) variants require `Me%ComputeAdvectionU` to be set and are **not on the
hot path** for the benchmark models — they did not appear in any Phase 3/4 hotspot list. Only fix
if a future model configuration surfaces them. (The hot variants are `_default_CG`, addressed in
Phase 1.)

---

## Open Issues / To Investigate (not scheduled)

### ✅ RESOLVED BUG — simulations blew up: `CheckStability - ModuleRunoff - ERR010`
**Status:** FIXED — the temporary Phase 6 **A/B measurement harness** was removed from `ModuleRunOff.F90`. Full details in repo memory `/memories/repo/mohid-fortran-perf.md` (section "RESOLVED BUG").

- **Symptom:** runs stopped with `CheckStability - ModuleRunoff - ERR010` (instability not recovered after restarts).
- **Root cause (confirmed):** the temporary `_AB` dual-call wrappers each ran their routine **twice** (baseline full-team + perf gated). The read-modify-write routines (`DynamicWaveXX/YY`, `UpdateWaterLevels`, …) are **not** pure-overwrite, so the second call corrupted `myWaterVolume`/flux/limiter state → instability. This was the regression, **not** the `NUM_THREADS` gating or the `SetMatrixValues` STATIC change (both `compare_mohid.py` zero-drift validated on the pre-harness build).
- **Fix applied:** reverted the 9 `ModifyRunOff` call sites to plain routine names; deleted the 7 `_AB` wrapper subroutines, the `RunOffForceThreads` module var + its baseline branch in `RunOffBoxThreads`, and the snapshot buffers. **Kept** the validated `NUM_THREADS(RunOffBoxThreads(...))` gating on all 9 OMP regions plus the harmless diagnostic counter + `KillRunOff` print (still answers "does the gate engage?"). Compiles clean.

### `forrtl: warning (526): IEEE_INVALID is signaling` — WithRain run
**Status:** OPEN — do not fix yet, investigate later.

- **Symptom:** WithRain model emits `forrtl: warning (526): IEEE_INVALID is signaling` (a warning at run/termination, not a trap — the run continues).
- **Likely cause:** Phase 2 relaxed FPE from `fpe0` to `fpe3` and added `/assume:noieee_compares`. Under `fpe3` the INVALID exception is *not trapped* but the IEEE flag can still be *raised* (e.g. a hardware compare against a stray NaN, or `0.0/0.0`, `sqrt(neg)`, `Inf-Inf` somewhere). The runtime reports the raised flag at exit as warning 526. Under the old `fpe0` this would have aborted at the point of occurrence instead.
- **Why it may be benign:** the `sqrt(max(..., 0.0))` guards already prevent `sqrt` of negatives; masked division may still produce a transient NaN/Inf in an inactive cell that never affects results (A/B compares + `compare_mohid.py` are passing).
- **To investigate:**
  1. Confirm it is new since Phase 2 (rebuild an `fpe0` debug build and see where it traps → pinpoints the offending operation).
  2. Check for unguarded divisions (`/ aux`, `/ velFace`, `/ Distance_Courant`) where the denominator can be exactly 0 for a computed cell.
  3. Decide: add a targeted guard, or accept as a benign warning and document it.
- **Do NOT** silence it globally by reverting FPE settings — that would reintroduce the IEEE-wrapper hotspots eliminated in Phase 2.

---

## Key Variable Types (Reference)

| Variable | Type | Notes |
|---|---|---|
| `HydraulicRadius` | `real` (real(4)) | Face hydraulic radius |
| `OverlandCoefficientX/Y` | `real` pointer array (real(4)) | Manning n |
| `VelModFaceU/V` | `real(8)` pointer array | Velocity modulus at face |
| `lFlowX/Y` | `real(8)` pointer array | Flow at face (m³/s) |
| `CriticalFlow` | was `real(8)`, now `real` in _default_CG | Critical flow limit |
| `Me%DX`, `Me%DY` | `real` (real(4)) scalar | Grid spacing (GridIsConstant=true) |
| `Gravity` | `real, parameter = 9.81` | In `ModuleGlobalData.F90` line 1758 |
| `U`, `V`, `Uaverage`, `Vaverage` | `real` (real(4)) | In `ComputeFaceVelocityModulus` |

---

## Build & Test Workflow

```
1. Build "Release Double OpenMP|x64" in VS2017
2. Run both test simulations (WithRain + NoRain)
3. Validate outputs against baseline (see Validation Tool below)
4. Build "Profile Double OpenMP|x64"
5. Run: D:\...\Performance_LargeModel_quickversion\exe\run_vtune_hotspot.bat
6. Open VTune result: vtune-gui <result_dir>
7. Compare hotspot timings to tables in Phase 1 Results section
```

---

## Validation Tool — `compare_mohid.py`

**Location:** `compare_mohid.py` (repo root)  
**Requirements:** `pip install h5py numpy` (Python 3.7+, h5py and numpy already installed)

Compares all output files between a baseline and a new run. Supports HDF5 (`.hdf5`) and MOHID text formats (`.dat`, `.srr`). Auto-detects `*_original` / `*_phase1` baseline pairs inside a `res/` folder.

### Usage

**Auto mode** — pass the `res/` folder; tool finds all `*_original` and `*_phase1` baseline pairs automatically:
```
python compare_mohid.py D:\...\Performance_LargeModel_quickversion\res
python compare_mohid.py D:\...\Performance_LargeModel_quickversion_noRain\res
```

**Explicit mode** — compare any two files or directories directly:
```
python compare_mohid.py D:\...\res\RunOff_55_original.hdf5  D:\...\res\RunOff_55.hdf5
python compare_mohid.py D:\...\res\Run55_original           D:\...\res\Run55
python compare_mohid.py D:\...\res\Run55_phase1             D:\...\res\Run55
```

**Options:**
```
--tol-abs 1e-5    Absolute tolerance (default). Both abs AND rel must exceed to FAIL.
--tol-rel 1e-4    Relative tolerance (default).
--verbose, -v     Show per-dataset detail and top-5 worst diffs in passing files.
--no-color        Disable ANSI colour.
```

**Exit codes:** `0` = all pass, `1` = one or more FAIL/ERROR, `2` = bad arguments.

### What it compares

| File type | Method |
|---|---|
| `.hdf5` | Every numeric dataset compared individually via h5py; max abs/rel diff reported |
| `.dat` | All floating-point numbers extracted and compared element-by-element |
| `.srr` | Same as `.dat` |
| `.fin`, `.map`, `.ver`, etc. | Skipped |

### Known baseline differences (Phase 1 + Phase 2 vs original)

Phase 1 `max()` guards are no-ops when physics is valid — confirmed zero diff vs original in both cases.

| Case | File | vs original | Notes |
|---|---|---|---|
| WithRain | All except `FloodPeriod.dat` | **PASS** | — |
| WithRain | `FloodPeriod.dat` | **FAIL** (9.3e-2 s, rel 2.8e-3) | Pre-existing before any optimisation; accepted |
| NoRain | All files | **PASS** (zero diff) | Arithmetic reformulations were reverted to achieve this |

The `FloodPeriod.dat` FAIL is a known and accepted Phase 1 side-effect (0.09 second shift in flood timing due to the CriticalFlow real4 → real8 precision change). For Phase 2 validation, the acceptance criteria are:
- All HDF5 datasets: max abs < 1e-5 (float32 noise threshold)
- `FloodPeriod.dat`: compare against Phase 1 baseline (`Run55_original` = Phase 1 results), not the original pre-Phase-1 baseline
- No new NaN/Inf values in any output

---

## Files Modified

| File | Change |
|---|---|
| `Software/MOHIDLand/ModuleRunOff.F90` | Phase 1 optimizations (9 changes) |
| `Software/MOHIDLand/ModuleRunOff.F90` | Phase 2 Step 1: 6 `sqrt(max(..., 0.0))` guards (L11421, L11447, L11681, L12900, L17383, L17429) |
| `Software/MOHIDLand/ModuleRunOff.F90` | Phase 2 Step 2a: `!DEC$ OPTIONS /fpe=3` wrapping `DynamicWaveYY_default_CG` (L12697 + L12956) |
| `Solutions/.../MOHIDLand/MOHIDLand.vfproj` | Added `Profile Double OpenMP|x64` config |
| `Solutions/.../MOHIDBase1/MOHIDBase1.vfproj` | Added `Profile Double OpenMP|x64` config |
| `Solutions/.../MOHIDBase2/MOHIDBase2.vfproj` | Added `Profile Double OpenMP|x64` config |
| `Solutions/.../MOHIDNumerics.sln` | Added `Profile Double OpenMP|x64` to solution |
| `D:\...\run_vtune_hotspot.bat` (×2) | VTune hotspot scripts (WithRain + NoRain) |
| `D:\...\run_vtune_threading.bat` (×2) | VTune threading scripts |
| `compare_mohid.py` | Output validation tool (HDF5 + dat/srr comparison) |
