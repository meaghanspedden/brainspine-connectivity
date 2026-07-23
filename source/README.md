# Brain–Spinal cord OPM coherence pipeline

MATLAB analysis pipeline for computing coherence between cortical sources,
spinal cord sources, and rectified EMG using OPM data. The pipeline uses DICS
beamforming (FieldTrip) and LCMV virtual electrodes to characterise
corticospinal connectivity in the beta band (10–35 Hz).

Data is available here:
https://zenodo.org/records/20824565
and
https://zenodo.org/records/20824827

---

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

