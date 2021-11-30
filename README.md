# fMRI_phantom_analysis
# **Purpose**
This script performs a temporal stability analysis of an fMRI timeseries acquired using a phantom. The outputs aid in the diagnosis of issues with scanner stability. A quality assessment of scanner stability should be performed prior to starting an fMRI study, particularly if the study is longitudinal or multi-site

The outputs include standard scanner temporal quality assessment metrics from Friedman & Glover, 2006, but we additionally incorporate a temporal Principal Component Analysis (PCA) for a more in-depth examination of separate components within the timeseries that explain the maximum amount of variance. By plotting the timeseries, Fourier spectrum, and spatial pattern of each component, the user can better diagnose the origin of unwanted signals that are corrupting their timeseries.

Note, this script was primarily designed for use with a mouse-sized agar phantom, which is reflected in the default options. If using with a human-sized phantom, we recommend changing the defaults.

# **Installation** 
Download from github into the desired folder:
`git clone https://github.com/CoBrALab/fMRI_phantom_analysis`

Create the environment from the provided environment file:`conda env create --file environment.yml`

# **Usage** 

## Data acquisition and processing
Construct a phantom with a similar size and relaxation properties to the brain of the species that you intend to scan. Scan the phantom using the same fMRI sequence that you intend to use in your study. 

It is recommended to run the phantom analysis on both the raw phantom timeseries, as well as on the timeseries after it is processed with the script that will be used in the fMRI study, then compare the results.

## The command line interface
Below are described the --help outputs from the command line interface
```
Usage: execute_analysis.sh [--repetition_time <arg>] [--input_roi <arg>] [--desired_slice <arg>] [--weisskoff_max_roi_width <arg>] [-h|--help] <input_epi> <output_path>
       
        <input_epi>: The phantom fMRI timeseries, as a 4D nifti.
        <output_path>: Path for where the output pdf report will be saved.
        --repetition_time: The acquisition TR, in seconds. (default: '1.0')
        --input_roi: A manually drawn single-slice ROI, as a 3D nifti. A 10x10 ROI in the center slice is used by default. (default: 'None')
        --desired_slice: The slice to be plotted in the report, as an integer. Center slice is computed and used by default. (default: 'None')
        --weisskoff_max_roi_width: The width of the largest ROI that is analyzed during Weisskoff analysis, in pixels. It should be as large as possible without extending outside the phantom (default: '20')
        -h, --help: Prints help
```
## Examples

No options:

`bash execute_analysis.sh phantom_timeseries_4D.nii.gz /home/phantom_analysis_results.pdf`

With options:

`bash execute_analysis.sh phantom_timeseries_4D.nii.gz /home/phantom_analysis_results.pdf --repetition_time 1.5 --input_roi phantom_custom_roi_3D.nii.gz --desired_slice 30 --weisskoff_max_roi_width 40`


## Guidance for choosing the options
#### When to provide a manually drawn input roi?

If your phantom is much bigger in diameter than 10 voxels (eg if you are using a human-sized phantom), then the default ROI of 10x10 may not provide an accurate representation of what's happening in the phantom, and you can draw your own. Alternatively, if you have a weirdly-shaped phantom, an unusually small phantom, or if you want to investigate a certain slice in particular, draw your own. The ROI should only be drawn in a SINGLE SLICE. *Note, the PCA is performed on the whole FOV, so if that's the only output you care about, you don't need to worry about ROIs.

#### When to provide a value for --desired_slice?

The script outputs a spatial map of the SNR, temporal noise, etc across a single slice of the phantom. If you wish the view a different slice other than the center one, you can specify it. This option only affects visualization.

#### When to provide a value for --weisskoff_max_roi_width?

For proper Weisskoff analysis, ROIs of increasing diameter are created in the phantom, going up to the largest possible ROI that is still entirely within the phantom. If you are using a large or small phantom, you should change the default max width. You will be able to see if your ROIs are an appropriate size in figure 3.

# Outputs
The script outputs a single pdf file, where each page contains a figure with the analysis results.

#### Spatial maps of signal and noise in slices of the phantom

Explanation of the figure: The figure displays a single slice in the phantom, you can change which slice is displayed with the --desired_slice option. The signal image is the mean across time within each voxel (before detrending). The temporal fluctuation noise image is the standard deviation across time (after detrending with second order polynomial).  The SFNR image is the signal image divided by the temporal fluctuation noise image (signal to fluctuation noise ratio). The static spatial image is the sum of the odd volumes, minus the sum of the even volumes (Friedman & Glover, 2006; Kayvanrad et al., 2021) The location of the ROI shows the positioning of the ROI that is used in future analyses - if the positioning does not look good, you should manually draw an ROI and provide it as input.

How to interpret the results: The signal, temporal fluctuation and SFNR are useful for comparing coils and scanners to decide which configuration is best for your fMRI study. It also helps to be aware of any spatial patterns in signal (for example, the CryoProbe exhibits a slight drop-off in signal and SFNR along the dorsal-ventral axis, but no variation in temporal noise). In the static spatial noise, there should not be any obvious spatial structure. 

![github1](https://user-images.githubusercontent.com/47565996/143964315-2a1c026d-1e4b-4f38-bb57-53959d9332fd.png)

#### Metrics within the ROI

Explanation of the figure: The timeseries was averaged across voxels within the ROI (shown in previous figure). This timeseries is shown in the first panel. A second order polynomial (first panel, orange) is fitted to the timeseries and subtracted to obtain residuals. The residuals are used to obtain the Fourier spectrum and histogram. Average summary values within the ROI are displayed at the bottom, such as the drift (change in the signal intensity over the scan, as a percentage of the mean signal) and the percentage fluctuation (standard deviation of the residuals as a percentage of the mean signal).

How to interpret the results:
Ideally, there will be no evident peak in the Fourier spectrum (unlike in the example figure). The histogram should approximate a Gaussian shape - with the QQ correlation plot quantifying the agreement between the data and a Gaussian distribution (if the points show a straight line and the correlation is high then the data is Gaussian). Also, the percent fluctuation value should be much smaller than the expected percent changes in the BOLD signal (much less than 1%, ideally <0.25%), otherwise you will not be able to distinguish BOLD activation from temporal scanner noise in your experiment.

*Note: the ROI provided should be one single slice! A multislice ROI might hide any peaks in the Fourier spectrum since different slices are acquired at slightly different times (especially in interleaved acquisitions).

![github2](https://user-images.githubusercontent.com/47565996/143964847-43fcd0cc-e4ca-42a6-b09c-21f086e718f8.png)

#### Ghosting analysis

Explanation of the figure: A mask of the phantom is created, then shifted by N/2 voxels to capture Nyquist ghosts. The ghost intensity in each voxel at the last timepoint is shown in the ghosting image, this can be compared to the signal in the background shown in the background image. The mean ghost to background ratio is shown within each timepoint.

![github3](https://user-images.githubusercontent.com/47565996/143964848-69e9a741-124f-4d95-b03a-ad0153c1b6f9.png)

#### Weisskoff analysis
Explanation of the figures: The location of each ROI used in the Weisskoff analysis is plotted in the first figure, the caption indicates the width of each ROI in voxels. If the default sizes don't suit your phantom, use this figure to estimate how much you should adjust the option --weisskoff_max_roi_width. The Weisskoff plot indicates the point at which a voxel loses statistical independence from neighbouring voxels due to low spatial frequency correlations caused by scanner instabilities. This radius of decorrelation is the ROI width at which the measured line deviates from the theoretical. the See Weisskoff, 1996 for detailed explanation.

How to interpret the results: A low radius of decorrelation (eg ~5 voxels) can be an indicator of low spatial frequency correlations caused by scanner instabilities. 
![github4](https://user-images.githubusercontent.com/47565996/143964849-7f8d8112-fc35-4132-aa05-9c669a069a27.png)
![github5](https://user-images.githubusercontent.com/47565996/143964853-3cac9735-f81c-4dca-984f-ab7bcfa6bf46.png)

#### Principal components analysis (PCA)
Explanation of the figures: Temporal PCA is applied on the entire timeseries (all voxels) to decompose the data into orthogonal components that explain the maxmium variance across time. In the upper figure, we plot the timeseries of the top 6 components that explain the most variance, as well as the Fourier spectrum of each timeseries. The amount of variance explained is indicated as a percentage next to the component number. In the lower figure, the spatial pattern of each component is plotted across all slices, with the color bar indicating how strongly weighted that component is in a given voxel.

How to interpret the results: The PCA is the most useful way to identify issues with scanner stability. It is immediately obvious from the timeseries and Fourier spectrum of each component, whether or not it represents random noise or suspicious scanner activity. For example, components 0-4 are suspicious since they have either a slowly-varying timeseries or a peak in the Fourier spectrum, whereas components 5-6 appear random. The cause of components 0-4 can be more easily diagnosed by examining their spatial patterns.
![github6](https://user-images.githubusercontent.com/47565996/143964854-3ca6cb3b-5fcc-47a5-993c-6bba30ba49b5.png)
![github7](https://user-images.githubusercontent.com/47565996/143964840-fe606782-4607-4784-8f52-93e0d922a484.png)


# **References**
Friedman, Lee, and Gary H. Glover. 2006. “Report on a Multicenter fMRI Quality Assurance Protocol.” Journal of Magnetic Resonance Imaging: JMRI 23 (6): 827–39.

Kayvanrad, Aras, Stephen R. Arnott, Nathan Churchill, Stefanie Hassel, Aditi Chemparathy, Fan Dong, Mojdeh Zamyadi, et al. 2021. “Resting State fMRI Scanner Instabilities Revealed by Longitudinal Phantom Scans in a Multi-Center Study.” NeuroImage 237 (August): 118197.

Weisskoff, R. M. 1996. “Simple Measurement of Scanner Stability for Functional NMR Imaging of Activation in the Brain.” Magnetic Resonance in Medicine: Official Journal of the Society of Magnetic Resonance in Medicine / Society of Magnetic Resonance in Medicine 36 (4): 643–45.
