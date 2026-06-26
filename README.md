# Brain–Spine MEG Coherence Pipeline

MATLAB analysis pipeline for computing coherence between cortical sources,
spinal cord sources, and rectified EMG using MEG data. The pipeline uses DICS
beamforming (FieldTrip) and LCMV virtual electrodes to characterise
corticospinal connectivity in the beta band (10–35 Hz).

---

## Dependencies

| Toolbox | Notes |
|---|---|
| [FieldTrip](https://www.fieldtriptoolbox.org/) | Source analysis, beamforming, visualisation |
| [SPM](https://www.fil.ion.ucl.ac.uk/spm/) | EEG/MEG data loading (`spm_eeg_load`, `spm2fieldtrip`) |
| [NeuroSpec 2.11](http://www.neurospec.org/) | Coherence spectra, directionality (`sp2a2_R2_mt`) |
| brainspineconnectivity (`bsc_path`) | Lab-internal utilities |
| FieldTrip `external/matplotlib` | `magma` colourmap |

All toolbox paths are set at the top of each script under **USER CONFIG**.

---

## Data

Each subject's MEG data is a preprocessed SPM `.mat` file following the naming
convention:

```
<data_root>/sub-<SubjectID>/ses-001/meg/
    pmergedoe1000mspddfflo45hi45hfcstatic_<run>_array1.mat
```

- Run is `001` for all subjects except `OP00224` (run `002`).
- Channel `EXG1` carries the EMG signal; it is rectified in-place before
  analysis.
- Trial labels: `1` = static contraction, `2` = rest.

### Subject lists

| Group | IDs | n |
|---|---|---|
| Brain (MEG array covers cortex) | OP00212, OP00213, OP00215, OP00219, OP00225, OP00221, OP00224 | 7 |
| Spine (MEG array covers torso/cord) | OP00212, OP00213, OP00215, OP00219, OP00220, OP00221, OP00224, OP00225, OP00226 | 9 |
| Both (coherence spectra) | Intersection of above | 7 |

---

## Geometry files

Both files are required by multiple steps. Update their paths under USER CONFIG
in each script.

| File | Contents |
|---|---|
| `geometries_experimental_withbrain.mat` | `sources_brain`, `sources_cent`, `mesh_brain`, `mesh_wm`, `mesh_torso`, `mesh_bone`, `mesh_lungs`, `mesh_heart` |
| `T.mat` | 4×4 affine matrix mapping native MEG space → MNI |
| `leadfield_experimental_bslaw_experimental.mat` | `leadfield_bs` — Biot-Savart leadfield for spine sources |

---

## Pipeline overview

```
step1  →  step2  →  step3  →  step4  →  step5
  │          │          │          │
Brain-     Spine-    Spinal    Brain       Coherence
EMG        EMG       VE        VE (M1)     spectra &
DICS       DICS      (LCMV)    (LCMV)      statistics
```

Run steps in order. Steps 3 and 4 depend on the saved results of steps 1 and 2
respectively. Step 5 requires both VEs from steps 3 and 4.

---

## Script reference

### `step1_brainEMG_coherence_DICS.m`

**Purpose:** Computes DICS beamformer coherence between cortical sources and
rectified EMG (10–35 Hz) for the 7 brain subjects.

**Method:**
- Single-shell head model on `mesh_brain`.
- Common spatial filter computed on the full (combined) frequency data; applied
  separately to static and rest conditions via label permutation test
  (500 permutations, `rng(1)`).
- Gaussian spatial smoothing (FWHM = 8 mm).
- Family-wise error correction: 95th percentile of the permutation max
  distribution.
- Coherence difference (`coh_diff = static − rest`) thresholded; peaks
  transformed to MNI space via `T.mat`.

**Key settings:** `lambda = '10%'`, `fband = [10 35]`, `numpermutation = 500`,
`fwhm_mm = 8`.

**Flags:**
- `run_subjects = 1` — run the full subject loop.
- `run_group_only = 1` — skip subject loop, load saved results, re-run group
  analysis only.

**Outputs (`save_dir`):**

| File | Description |
|---|---|
| `brainResult_sub<ID>_brain_pct10.mat` | `coh_diff`, `thr95`, `mask`, `invp`, `x_mni` per subject |
| `groupRes_brain_DICS_brain_pct10.mat` | `subjResults` struct for all 7 subjects |

**Figures (`save_dir/figures`):**

| File | Description |
|---|---|
| `step1_sub<ID>_brainEMG_coherence_brain_pct10.fig` | Per-subject cortical coherence map |
| `step1_group_brain_max_brain_pct10.fig` | Group peak location on brain mesh |
| `step1_group_brain_prevalence_brain_pct10.fig` | Group prevalence map (threshold = 0.3) |

---

### `step2_spineEMG_coherence_DICS_v2.m`

**Purpose:** Computes DICS coherence between spinal cord sources and rectified
EMG for the 9 spine subjects, using the Biot-Savart (BS) leadfield.

**Method:**
- Infinite-medium volume conductor (appropriate for BS leadfield).
- Label-based channel matching between MEG data and leadfield (handles
  per-subject bad channel removal cleanly).
- Gaussian smoothing (FWHM = 20 mm); common spatial filter + 500-permutation
  label permutation test.
- Dominant source orientation extracted via SVD of the DICS filter at the peak
  source; stored as `dom_ori` [x; y; z] unit-vector components.
- Group prevalence map; connected-component clustering on prevalence ≥ 0.2 to
  define the VE ROI for step 3.

**Additional figures (vs v1):**
- Per-subject null maxima location histogram.
- Participant-1-only null distribution stratified by dominant orientation axis
  (x = right–left, y = cranial–caudal, z = dorsal–ventral).
- Group dominant orientation bar chart ("Figure C").

**Key settings:** `lambda = '10%'`, `fband = [10 35]`, `numpermutation = 500`,
`fwhm_mm = 20`, `out_suffix = '_BS'`.

**Outputs (`save_dir`):**

| File | Description |
|---|---|
| `subResult_sub<ID>_BS.mat` | `coh_diff`, `thr95`, `mask`, `invp_smooth`, `pvals` per subject |
| `groupRes_spine_DICS__BS.mat` | `subjResults` struct including `dom_ori` for all 9 subjects |
| `cluster_spineEMG_pos_BS.mat` | `ROIpos` — largest connected prevalence cluster; input for step 3 |

---

### `step3_VE_BS_prevalenceROI.m`

**Purpose:** Builds a spinal cord virtual electrode (VE) for each subject using
LCMV beamforming, centred on the data-driven prevalence-cluster ROI from step 2.

**Method:**
- ROI loaded from `cluster_spineEMG_pos_BS.mat`; centre and bounding radius
  computed.
- LCMV spatial filter fitted to the full data covariance; `ft_virtualchannel`
  with SVD extracts the dominant component within the ROI sphere.

> **Note:** The output suffix `_BS1` is used for initial reproducibility
> testing. Rename to `_BS` (drop the `1`) once confirmed, so that step 5 can
> find the files.

**Outputs (`save_dir`):**

| File | Description |
|---|---|
| `VE_spine_prevalence_sub<ID>_forspectra_BS1.mat` | `VE`, `roi_center`, `R`, `ROIpos` per subject |

---

### `step4_brain_VE_M1.m`

**Purpose:** Builds a cortical (M1) virtual electrode for each of the 7 brain
subjects using LCMV beamforming at the intersection of the group brain-EMG
significant region and an anatomical M1 sphere.

**Method:**
- Group brain-EMG prevalence map (from step 1, threshold = 0.3) identifies
  significantly active cortical voxels.
- M1 sphere (radius 20 mm) centred on MNI `[−38 −26 56]` (left M1 hand knob)
  is transformed to native MEG space via `T.mat`.
- Intersection of sphere and active voxels defines the ROI; falls back to the
  full sphere if the intersection is empty.
- Single-shell head model; LCMV filter fitted to full data covariance;
  `ft_virtualchannel` with SVD.

**Key settings:** `M1_mni_centre = [-38 -26 56]`, `M1_radius_mm = 20`.

**Outputs (`save_dir`):**

| File | Description |
|---|---|
| `M1_ROI_brain_pct10.mat` | `intersection_pos`, `roi_native`, M1 sphere parameters |
| `sub<ID>_VE_brain_M1_brain_pct10.mat` | `VE_brain`, `roi_center`, `R`, `idx_roi` per subject |

---

### `step5_coherence_spectra.m`

**Purpose:** Downstream coherence spectra, directionality, PSI, and statistical
analysis using the brain (M1) and spinal cord VEs from steps 3–4.

**Method:**
- Three signal pairs: Brain–EMG, Brain–Spine, Spine–EMG.
- NeuroSpec `sp2a2_R2_mt` for multitaper coherence spectra and cross-correlation.
- Multitaper partial coherence (custom `mt_partial_coherence`) controlling for
  the third channel.
- Phase slope index (`ft_connectivityanalysis`, method `psi`) for directionality.
- Halliday R2 forward/reverse areas for directionality comparison.
- FOOOF periodic power extraction.
- Paired statistics: Wilcoxon signed-rank and paired t-tests.
- Saves `brainspine_boxplot_data_bslaw_prevalence.mat` for figure regeneration.

**Key settings:** `fband = [10 35]`, `seg_pwr = 11` (segment length 2^11 = 2048
samples), `lat_min/max = ±50 ms`.

**Outputs (`save_dir/figures/spectra_bslaw_prevalence`):**

| File | Description |
|---|---|
| `subOP00212_coherence_spectra_*.fig/.png` | Participant 1 three-pair coherence spectra |
| `group_coherence_spectra_global_*.fig/.png` | Group spectra aligned to each pair's own peak |
| `brainspine_peak_vs_threshold_boxplot_*.fig/.png` | Peak coherence / threshold boxplot |
| `partial_coherence_pct_reduction_*.fig/.png` | Brain-EMG reduction after partialising cord |
| `group_directionality_comparison_*.fig` | Halliday R2 vs PSI panel |
| `group_peak_latencies_*.fig` | Cross-correlation peak latency plot |
| `SuppTable1b_peak_coherence_BS.csv` | Peak coherence table (multitaper + NeuroSpec) |
| `SuppTable_partial_coherence_pct_reduction.csv` | Partial coherence % reduction per subject |
| `brainspine_boxplot_data_bslaw_prevalence.mat` | Saved data for boxplot regeneration |

---

### `GET_FIGURES.m`

**Purpose:** Regenerates all pipeline figures from pre-saved `.mat` results,
without re-running the computationally expensive beamforming steps.

**Usage:** Edit the USER CONFIG section (paths, smoothing FWHMs, `plot_step*`
flags, `saveFigs`), then run. Each step flag can be toggled independently.

**Required input files (all in `save_dir`):**

| File | Used by |
|---|---|
| `groupRes_brain_DICS_brain_pct10.mat` | Step 1 figures |
| `groupRes_spine_DICS__BS.mat` | Step 2 figures |
| `cluster_spineEMG_pos_BS.mat` | Step 3 figure |
| `brainspine_boxplot_data_bslaw_prevalence.mat` | Step 4 (boxplot) figure |
| `geometries_experimental_withbrain.mat` | All steps |
| `T.mat` | Steps 1, 4 |

**Smoothing options:** `doSmooth = 1` enables post-hoc re-smoothing of saved
results at configurable FWHMs (`spine_smooth_fwhm_mm`, `brain_smooth_fwhm_mm`).
Set to `0` to plot results as originally computed.

---

## Reproducibility

`rng(1)` is called at the start of steps 1 and 2 to fix the permutation test
random seed. Results should be bit-for-bit reproducible given the same input
data and toolbox versions.

---

## Typical run order

```matlab
% 1. Brain-EMG coherence maps
run step1_brainEMG_coherence_DICS.m

% 2. Spine-EMG coherence maps + ROI cluster
run step2_spineEMG_coherence_DICS_v2.m

% 3. Spinal virtual electrode (rename _BS1 → _BS outputs when confirmed)
run step3_VE_BS_prevalenceROI.m

% 4. M1 virtual electrode
run step4_brain_VE_M1.m

% 5. Coherence spectra and statistics
run step5_coherence_spectra.m

% Regenerate figures only (after any of the above)
run GET_FIGURES.m
```
