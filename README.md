# Psilocybin-shapes-propagation-of-brain-activity
Codes corresponding to the article "Psilocybin shapes the slow, global propagation of brain activity over the cortical layout of 5HT2a receptors"

The code for detection of travelling waves is provided as the BrainWaves Toolbox: https://github.com/VeronicaMaki-Marttunen/BrainWavesToolbox

The code for functional connectivity analyses is included here.

The code has been developed in MATLAB R2024b (Mathworks). In addition, FieldTrip should be in the path, as well as the BrainSpace Toolbox for the gradient analysis.

The code loads the fMRI data which should be in surface space.

run_conn.m: creates FC matrix; calls calc_FC_function.m
run_FCchange.m: calculates FC change between sessions
run_gradients.m: calculates spatial gradients

If you use this code, please cite the source publication:
Mäki-Marttunen, V. (2026) Communications Biology
