# brainspine-connectivity
Preproc and source recon for concurrent brain spinal cord OPM recordings.  


Script overview:
1. preproc_spinecoh
2. prep for source recon:setup_source; create_leadfields; (run once,not pr subject) beamforming (spinal cord only)
3. source_freq: analyses data pr. subject. t test and permutation test for difference between rest and contraction for each spinal cord source point.
4. group_analysis: visualises prevalence of significant results over subjects (based on permutation testing primarily)

Group analysis
source_freq_group calculates plots and saves power differences pr spinal cord source point between rest and contraction per TRIAL for further analysis 
in group_analysis_classical (t test of power difference across subjects per spinal cord source point)
group_analysis_classical_fixed concatenates across subjects and looks for fixed effect

Extra control analyses:
bootstrap_peak_loc: bootstraps to see if peak location moves

Connectivity:
getCoh - sensor level coherence and power modulations brain muscle
getPowMod - spectral power modulations (both directions) brain channels and EMG. also looks at power corr. single subject, saves for group analysis
This also does CVA RN

Have added
beamforming_FT
and min_norm_FT 
using infinite volume conductor