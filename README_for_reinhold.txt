This README pertains to Reinhold et al., 2023.
This points to high-level scripts for running analyses and generating figures, and also any example datasets.



GENERAL METHODS AND ANALYSES
1. Design of behavior rig
Description of behavior rig and instructions to build rig here: https://github.com/kimerein/reaching-behavior-rig

2. Control of Arduino behavior rig
All scripts to control the Arduino here: https://github.com/kimerein/behaviorRig

3. GLM
Code here: https://github.com/kimerein/k-glm
Start at kim_run_glm.py.
Relevant example datasets: in forglm_trainingSet_wreach, behEvents.mat, neuron_data_matrix.mat, timepoints.mat, unitnames.mat
Ran in VSCode on Windows 10 Pro, 64 bit.
Expected output: See code. Expected run time, given example datasets: <10 mins.

4. Tensor regression
Code here: https://github.com/kimerein/tensor_regression
Start at script demo_tensorRegression_forKim.ipynb.
Relevant example datasets: allLabels.mat, tensor.mat, timepoints_for_tensor.mat at example datasets to run code\tensor_regression
Ran in VSCode on Windows 10 Pro, 64 bit.
Expected output: See code. Expected run time, given example datasets: <10 mins.

5. Behavior video analysis
All code in repos: https://github.com/kimerein/reach-behavior-analysis and https://github.com/kimerein/reachBehavior
Start at prepForOrchestra_analyzeReachVideo.m and analyzeReachVideo.m.
Relevant example datasets: in example_video_behavior_data, 2011-07-28 02-57-15-C.AVI and 2011-07-28 02-57-15-C_OUTPUT.TXT, and 2011-07-28 03-27-03-C.AVI and 2011-07-28 03-27-03-C_OUTPUT.TXT
Ran on Linux, Python, Matlab R2021b.
Expected output: See code. Expected run time, given example datasets: <30 mins.

6. After processing videos, further analyze behavior data
All code here: https://github.com/kimerein/integrate-phys-and-beh
Start at script_for_reaching_rate_analysis.m. 
Note that only the top sections of this script were used. Bottom of script contains defunct code.
Relevant example datasets: in behavior_analysis
To run, change exptDataDir and mouseDBdir at beginning of script_for_reaching_rate_analysis.m to point to folder alltbt28Aug2023152404 and Combo Behavior Log - Slimmed down w old mice added.csv.
Ran in Matlab R2021b on Windows 10 Pro, 64 bit.
Expected output: See code. Expected run time, given example datasets: <10 mins.

7. Photometry acquisition
Code here: https://github.com/kimerein/photometry

8. Photometry analysis
Code here: https://github.com/kimerein/integrate-phys-and-beh
See the function processPhotometry.m.
Relevant example datasets: in example datasets to run code\photometry
Pass this directory (example datasets to run code\photometry) as first argument to processPhotometry.m. Second argument, 'Marci'.
Expected output: See code. Expected run time, given example datasets: <10 mins.

9. Single unit data analysis
All code here: https://github.com/kimerein/integrate-phys-and-beh
Start at process_units_for_reaching_behavior.m.
Ran in Matlab R2021b on Windows 10 Pro, 64 bit.

10. Figure 1
Code here: https://github.com/kimerein/integrate-phys-and-beh
script_for_reaching_rate_analysis.m 
Ran in Matlab R2021b on Windows 10 Pro, 64 bit.

11. Figure 2
Code here: https://github.com/kimerein/integrate-phys-and-beh
script_for_reaching_rate_analysis.m 
Ran in Matlab R2021b on Windows 10 Pro, 64 bit.

12. Figure 3
Code here: https://github.com/kimerein/integrate-phys-and-beh
script_for_reaching_rate_analysis.m 
Ran in Matlab R2021b on Windows 10 Pro, 64 bit.

13. Figure 4
Code here: https://github.com/kimerein/integrate-phys-and-beh
scriptToMakeOutcomeFigure4.m 
Ran in Matlab R2021b on Windows 10 Pro, 64 bit.

14. Figure 5
Code here: https://github.com/kimerein/integrate-phys-and-beh
makeFig5figs.m
getErrorBars_onShuffleConsensus.m
getErrorBars_onShuffleConsensus_trialbytrial.m
Ran in Matlab R2021b on Windows 10 Pro, 64 bit.



OTHER CODE FROM OTHER RESEARCHERS
1. UltraMegaSort for spike sorting
UltraMegaSort is described elsewhere: https://neurophysics.ucsd.edu/lab/UltraMegaSort2000%20Manual.pdf

2. DeepLabCut
Described elsewhere: https://github.com/DeepLabCut
Relevant example datasets: in example_video_behavior_data\high_speed, vid_2020-12-28-175235-0000.avi, vid_2020-12-28-175235-0001.avi, vid_2020-12-28-175235-0002.avi

3. Physiology acquisition using SpikeGLX
Described elsewhere: https://billkarsh.github.io/SpikeGLX/



INSTALLATION GUIDE
All source code can be accessed at https://github.com/kimerein. Also required: Python, Matlab