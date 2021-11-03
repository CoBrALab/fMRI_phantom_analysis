#!/bin/bash
source activate phantom_analysis

#define arguments, if an empty string is provided (e.g. for the manually defined ROI), the default argument is interpreted as 'None' in python
input_epi=$1
input_roi=${2:-None}
output_path=$3
desired_slice=${4:-None}
TR=$5
weisskoff_max_roi_width=${6:-20}

#get the path to the folder where the current script is located
wdir="$PWD"; [ "$PWD" = "/" ] && wdir=""
case "$0" in
  /*) scriptdir="${0}";;
  *) scriptdir="$wdir/${0#./}";;
esac
scriptdir="${scriptdir%/*}"

python3.8 $scriptdir/phantom_analysis_functions.py $input_epi $input_roi $output_path $desired_slice $TR $weisskoff_max_roi_width
