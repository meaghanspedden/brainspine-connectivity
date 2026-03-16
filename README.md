# brainspine-connectivity

Analysis code for:

> Spedden ME, Schmidt M, O'Neill GC, West TO, Mellor S, Tierney TM, Alexander NA, Puvvada S, Callaghan M, Farmer SF, Bestmann S, Barnes GR. *Coherent oscillations across the human brain, spinal cord, and muscle network.*

Data are available on Zenodo: [doi: 10.5281/zenodo.19029856](https://doi.org/10.5281/zenodo.19029856) and [doi: 10.5281/zenodo.19004628](https://doi.org/10.5281/zenodo.19004628)

---

## Overview

This repository contains MATLAB code for concurrent OPM brain and spinal cord and EMG data preprocessing and source analysis. The pipeline detects coherent oscillatory activity (10–35 Hz) across the brain–spinal cord–muscle axis during voluntary hand contraction, using DICS beamforming with a five-compartment spinal cord boundary element model (BEM).

---

## Requirements

- **MATLAB** R2021b or later
- **SPM** (development version): https://github.com/spm/spm
- **FieldTrip**: https://github.com/fieldtrip/fieldtrip
- **NeuroSpec 2.11**: https://github.com/dmhalliday/NeuroSpec
- **Helsinki BEM Framework**: https://github.com/MattiStenroos/hbf_lc_p
- **msg_coreg toolbox**: https://github.com/maikeschmidt/msg_coreg

## Repository structure
 
```
brainspineconnectivity/
├── preprocessing/
│   ├── preproc_spinecoh.m   # Step 0: OPM and EMG preprocessing
│   └── ...                  # Helper functions for preprocessing
└── source/
    ├── RUN_PIPELINE.m       # Steps 1–4: source-level coherence analysis
    ├── RUN_SPECTRA.m        # Steps A–B: virtual electrode spectra and directed coherence
    └── ...                  # Helper functions for source analysis and plotting
```
 
---

## Usage

Scripts should be run in order: **preproc_spinecoh.m** → **RUN_PIPELINE.m** → **RUN_SPECTRA.m**

### Step 0 — Preprocessing (`preproc_spinecoh.m`)

Preprocesses raw OPM (`.lvm`) and EMG (`.vhdr`/`.eeg`/`.vmrk`) data for one participant at a time.

Set the subject ID at the top of the script:
```matlab
sub = 'OP00212';
```

**What it does:**
- Loads raw OPM data and sensor positions (`.tsv`)
- Computes and inspects power spectra; identifies and removes bad channels
- Applies homogeneous field correction (HFC)
- Bandpass filters OPM and EMG data (5–45 Hz) with additional notch filters at 50 and 100 Hz
- Estimates and removes cardiac artefacts using signal space projection (SSP)
- Downsamples OPM data to match EMG sampling rate (1000 Hz)
- Applies synthetic gradiometry using a reference sensor
- Synchronises OPM and EMG data via TTL pulse; crops to matched time windows
- Epochs data into 1-second non-overlapping trials and removes outlier trials
- Merges all runs and conditions into a single SPM file per participant

**Output:** A merged, preprocessed SPM `.mat`/`.dat` file per participant, saved to `save_dir`.

---

### Steps 1–4 — Source coherence pipeline (`RUN_PIPELINE.m`)

Master script for DICS-based source localisation of coherent activity. Edit the **USER CONFIG** section at the top before running.

**Key paths to set:**
```matlab
cfg_pipeline.fieldtrip_path  = '...\fieldtrip';
cfg_pipeline.spm_path        = '...\spm';
cfg_pipeline.bsc_source_path = '...\brainspineconnectivity\source';
cfg_pipeline.data_root       = '...';          % root folder containing sub-XX directories
cfg_pipeline.geomfile        = '...\geometries_cervical_realistic.mat';
cfg_pipeline.lf_v2_path      = '...\bem_v2_leadfield_cervical_realistic_bem_.mat';
cfg_pipeline.T_mat_path      = '...\T.mat';    % MNI-to-native transformation matrix
cfg_pipeline.save_dir        = '...';          % output directory
```

**Steps (each can be toggled on/off):**

| Flag | Step | Description |
|------|------|-------------|
| `run_step1` | Brain–EMG DICS | Localises brain sources coherent with EMG; permutation testing (n=500); group prevalence map |
| `run_step2` | Spine–EMG DICS | Localises spinal cord sources coherent with EMG using BEM v2 lead fields; null distribution diagnostics for OP00212 |
| `run_step3` | Spinal virtual electrode | Defines spinal ROI from Step 2 prevalence cluster; builds virtual electrode time series per participant using LCMV |
| `run_step4` | SpineVE–Brain DICS | Localises brain sources coherent with the spinal virtual electrode; M1 ROI overlap analysis |

**Step dependencies:**
- Step 3 requires `cluster_spineEMG_pos_bemv2*.mat` from Step 2
- Step 4 requires `VE_spine_sub*_forspectra_bemv2*.mat` from Step 3

**Key parameters:**
```matlab
cfg_pipeline.fband          = [10 35];   % frequency band (Hz)
cfg_pipeline.numpermutation = 500;       % permutations for significance testing
cfg_pipeline.doSmooth       = 1;         % Gaussian spatial smoothing
cfg_pipeline.spine_smooth_fwhm_mm = 20; % smoothing kernel for spinal maps
cfg_pipeline.brain_smooth_fwhm_mm = 8;  % smoothing kernel for brain maps
```

---

### Steps A–B — Spectra and directed coherence (`RUN_SPECTRA.m`)

Downstream analysis of virtual electrode time series. Requires outputs from `RUN_PIPELINE.m`. Edit the **USER CONFIG** section before running.

**Steps:**

| Flag | Step | Description |
|------|------|-------------|
| `run_stepA` | Brain virtual electrode | Builds brain virtual electrode from M1 ROI using LCMV beamformer |
| `run_stepB` | Coherence spectra | Pairwise NeuroSpec coherence (Brain–EMG, Brain–Spine, Spine–EMG): directed coherence, normalised cumulant density (cross-correlation), and peak latencies; FieldTrip coherence spectra; 10–35 Hz band power comparison (contraction vs rest) |

**Key parameters to match with `RUN_PIPELINE.m`:**
```matlab
cfg.doSmooth             = 1;
cfg.spine_smooth_fwhm_mm = 20;   % must match RUN_PIPELINE
cfg.M1_roi_suffix = 'functionalVE_spineSmooth_20mm_brainSmooth_8mm';
```

---

## Data format

Raw and preprocessed data are expected in the following structure (see Zenodo for the dataset):

```
sub-XX/
└── ses-001/
    ├── emg/
    │   ├── *.vhdr
    │   ├── *.eeg
    │   └── *.vmrk
    └── meg/
        ├── *_positions.tsv
        ├── *_raw.lvm
        ├── *_preprocessed.dat
        └── *_preprocessed.mat
```

Shared forward model files (BEM meshes, lead fields, MNI-to-native transformation) are in the `Leadfields_meshes/` folder at the top level of the Zenodo repository.

---

## Notes

- Participant OP00212 is the reference participant with full anatomical co-registration. Results for this participant are presented alongside group-level analyses in the paper.
- Subject lists for brain and spinal analyses differ: brain analyses use n=7, spinal analyses use n=9 (see paper for details).
- Analysis was performed in a magnetically shielded room at UCL. The preprocessing pipeline is tailored to QuSpin Neuro-1 triaxial OPM data acquired with the custom scanner-cast described in the paper.

---


