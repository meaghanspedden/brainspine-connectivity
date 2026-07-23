# Brain–Spinal cord OPM Coherence Pipeline

MATLAB analysis pipeline for computing coherence between cortical sources,
spinal cord sources, and rectified EMG using OPM data. The pipeline uses DICS
beamforming (FieldTrip) and LCMV virtual electrodes to characterise
corticospinal connectivity in the beta band (10–35 Hz).

---Data is available here:
https://zenodo.org/records/20824565
and
https://zenodo.org/records/20824827


## Dependencies

| Toolbox | Notes |
|---|---|
| [FieldTrip](https://www.fieldtriptoolbox.org/) | Source analysis, beamforming, visualisation |
| [SPM](https://www.fil.ion.ucl.ac.uk/spm/) | EEG/MEG data loading (`spm_eeg_load`, `spm2fieldtrip`) |
| [NeuroSpec 2.11](http://www.neurospec.org/) | Coherence spectra, directionality (`sp2a2_R2_mt`) |
| FieldTrip `external/matplotlib` | `magma` colourmap |


## Setup

Before running anything, edit `source/brainspine_config.m` with your local
paths — this is the only file that needs changing to run the pipeline on a
new machine. It ships with the paths used on the original analysis machine
as a template; update each field to match your own setup:
- `fieldtrip_path`, `spm_path`, `neurospec_path` — toolbox locations
- `data_root` — raw/preprocessed data directory
- `leadfields_dir` — leadfield/geometry files (`geomfile`, `T_mat_path`, `lf_path` are derived from this)
- `save_dir` — pipeline output directory
- `figures_out_dir` — output directory used by `GET_FIGURES.m`

All other repo-internal paths (`source/`, `plotting/`, `preprocessing/`) are
resolved automatically at runtime, so the repo can be cloned to any location.


## Repository structure

| Folder | Purpose |
|---|---|
| `preprocessing/` | Raw data preprocessing (heartbeat removal, HFO rejection, trial selection) |
| `source/` | Current 5-step analysis pipeline (Biot-Savart leadfields) |
| `plotting/` | Visualisation utilities (colormaps, topoplot, shaded error bars) |
| `stats/` | Statistical utilities (FDR correction, Bayesian prevalence) |


## Usage

### 1. Preprocessing

Run `preprocessing/PREPROC_SPINECOH_V1.m` before the source analysis pipeline.
Key steps:
- Heartbeat estimation and removal (`est_heartbeat.m`, `remove_heartbeat.m`)
- High-frequency oscillation component removal (`hfo_project_out_comps.m`)
- Outlier trial rejection (`spm_opm_removeOutlierTrials.m`)

Per-subject bad channel lists are stored in `preprocessing/channels_removed/` (subjects OP00212–OP00226).

### 2. Source analysis

Scripts live in `source/`. Run steps in order.

```
step1  →  step2  →  step3  →  step4  →  step5
  │          │          │          │
Brain-     Spine-    Spinal    Brain       Coherence
EMG        EMG       VE        VE (M1)     spectra &
DICS       DICS      (LCMV)    (LCMV)      statistics
```

Steps 3 and 4 depend on the saved results of steps 1 and 2 respectively. Step 5 requires both VEs from steps 3 and 4.

### 3. Figures

`source/GET_FIGURES.m` generates all manuscript figures from the saved pipeline outputs.


## Utilities

**Plotting** (`plotting/`): colormaps (`brewermap`, `viridis`), `shadedErrorBar`, `topobrain`, `plotBEM`.

**Stats** (`stats/`): `fdr_bh.m` (Benjamini-Hochberg FDR correction), Bayesian prevalence posterior plots.
