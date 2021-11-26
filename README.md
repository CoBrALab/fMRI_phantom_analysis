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

## Data acquisition
Construct a phantom with a similar size and relaxation properties to the brain of the species that you intend to scan. Scan the phantom using the same fMRI sequence that you intend to use in your study.

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

# **References**
Friedman & Glover. Report on a multicenter fMRI quality assurance protocol. 2006.
