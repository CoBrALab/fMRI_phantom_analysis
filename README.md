# fMRI_phantom_analysis
# **Purpose**
This script performs automated analysis of a timeseries acquired using a phantom. The purpose is to assess the temporal stability of the scanner/sequence prior to the start of an fMRI experiment, this is accomplished by calculating the SNR, percent fluctuation, drift, Fourier spectrum.

 Many of the metrics are taken from (Friedman & Glover, 2006), but we additionally incorporate a temporal Principal Component Analysis (PCA) for a more in-depth examination of separate components within the timeseries that explain the maximum amount of variance. By plotting the timeseries, Fourier spectrum, and spatial pattern of each component, the user can better diagnose the origin of unwanted signals that are corrupting their timeseries.

 Note, this script was primarily designed for use with a mouse phantom. If using with a human-sized phantom, we recommend manually drawing an ROI of an appropriate size to provide as input.

# **Usage** 

## Create the environment
`conda env create --file environment.yml`

## Decide whether you want to manually draw an ROI
If no roi.mnc file is provided in the arguments, the script will automatically define a 10x10 ROI in the center-slice of the phantom. When it is doing so, it assumes that the phantom is approximately in the center of the FOV. The summary metrics shown will be for this ROI. If you wish to extract summary metrics from a specific slice/region, then you can draw your own ROI and save it as a minc file. Drawing your own ROI is also useful if you want it to be located in the same spot within the phantom across different sessions. Important: your ROI should only be drawn in a single slice.

## Run the analysis
`execute_analysis.sh timeseries roi output_filepath slice_to_plot TR `

The roi and the slice_to_plot are optional. If not specifying an ROI, and if you don't wish to look at a particular slice, then provide an empty string. 
* The timeseries should be a nifti file, and should be provided as a string. 
* The roi should be a minc file, provided as a string. 
* The output_filepath is the path and name of a pdf file where the outputs will be stored. Provide as a string.
* The slice_to_plot is an integer slice value. 
* The TR is a float (in seconds). 

Example with an ROI and slice_to_plot provided:

`execute_analysis.sh "phantom_timeseries_4D.nii.gz" "roi_3D.mnc" "/home/phantom_analysis_results.pdf" 12 1.0`

Example with no ROI and slice_to_plot provided:

`execute_analysis.sh "phantom_timeseries_4D.nii.gz" "" "/home/phantom_analysis_results.pdf" "" 1.0`

## Outputs
The script outputs a single pdf file, where each page contains a figure with the analysis results.

# **References**
Friedman & Glover. Report on a multicenter fMRI quality assurance protocol. 2006.
