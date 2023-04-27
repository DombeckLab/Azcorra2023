"# Azcorra2023 - Fiber photometry data analysis" 

Azcorra, Gaertner et al. Nat Neuro 2023 custom MATLAB scripts and data.

In our manuscript (Azcorra, Gaertner et al. Nat Neuro 2023), we analyzed fiber photometry recordings from head-fixed mice running on a treadmill while receiving rewards and aversive air puffs in order to analyze differences or similarities in functional responses of different subtypes of striatonigral dopamine neurons. Custom code necessary fiber photometry data pre-processing (bleaching correction, DF/F calculation, behavioral parsing...) and for data analysis/figure generation in the paper are provided.


Code is organized in 3 folders:

1) "Data pre processing"- scripts necessary to obtain the 'data6' table used in all scripts in the following folder. This 'data6' is also available on Zenodo (see below)
	'addData_paper.m'	runs all other scripts in a row, but each can be run independently on existing data. 
	'concatPhotom405.m'	concatenates raw data collected on Picoscope, separates the fluorescence due to 405 vs 470 nm illumination (405 nm is GCaMP's isosbestic point and thus serves as a movement control, and re-bins the data from 4 kHz to 100 Hz.
	'data_processing.m'	calculates DF/F from raw fluorescence by correcting for fiber/setup/tissue autofluorescence and bleaching, and then normalizes transients from 0 to 1 (saving the normalization value to allow reverting it). It also converts velocity from Volts as read out from the rotary encoder to m/s.
	'selectSignals_paper.m'	selects behavioral events used to obtain triggered averages (rewards, air puffs, accelerations...)
	'exclusionCriteria_allFigs.m'	determines exclusion criteria for each recording to decide which will be included for analysis (running time, signal-to-noise ratio)	

2) "Data analysis for Azcorra 2023 figures" - scripts to generate all the Figures and Extended Data Figures in the manuscript.
	'LocomotionAnalysis.m'	all figures analyzing locomotion signaling of dopamine subtypes (Fig. 2, 3A, and Extended Data Fig. 1F,H, 6, 7, 9H-K).
	'RewardAnalysis.m'	all figures analyzing reward/air puff signaling of dopamine subtypes and those analyzing both locomotion and reward/air puff (Fig. 4, 5C-E, 6, and Extended Data Fig. 1D,E,G,I, 8, 9A-G).
	'SNcVsStrAnalysis.m'	all figures comparing signaling in SNc vs striatum (Fig. 7 and Extended Data Fig. 10)
	'AccDecSizeAnalysis.m'	analysis of the response to accelerations/decelerations of different sizes (Fig. 5A-B)
	'FiberLocations.m'	all figures showing recording locations (Fig. 2E, 3A, 4J-K and Extended Data Fig. 6A,J, 8J-K)

3) "Code bits - necessary for main scripts" - custom functions used in several of the main scripts in the previous two folders. 


The data necessary to generate the figures found in the manuscript ('data6.mat') is available at Zenodo (DOI: 10.5281/zenodo.7871982). This data has already been pre-processed (using the code in the "Data pre processing" folder) to obtain normalized DF/F traces (though de-normalized traces can be easily recovered using the 'norm' variable in the table). This dataset also includes metadata not present in the raw data files, as well as XYZ coordinates for most recordings obtained from histology. 
The raw data used to generate data6 is also available at Zenodo (DOI: 10.5281/zenodo.7871634)


This scripts was written for MATLAB R2021a (compatibility with newer versions is not guaranteed). 
This code requires the following toolboxes:
	'Signal Processing Toolbox'
	'Image Processing Toolbox'
	'Statistics and Machine Learning Toolbox'
	'Curve Fitting Toolbox'
	'Econometrics Toolbox'
This code also requires the following add-ons:
	'fastsmooth.m' 		Tom O'Haver (2023). Fast smoothing function (https://www.mathworks.com/matlabcentral/fileexchange/19998-fast-smoothing-function), MATLAB Central File Exchange. Retrieved April 27, 2023.
	'customcolormap.m'	Víctor Martínez-Cagigal (2023). Custom Colormap (https://www.mathworks.com/matlabcentral/fileexchange/69470-custom-colormap), MATLAB Central File Exchange. Retrieved April 27, 2023.
	'plotSpread.m'		Jonas (2023). plot spread points (beeswarm plot) (https://www.mathworks.com/matlabcentral/fileexchange/37105-plot-spread-points-beeswarm-plot), MATLAB Central File Exchange. Retrieved April 27, 2023.
	

	




