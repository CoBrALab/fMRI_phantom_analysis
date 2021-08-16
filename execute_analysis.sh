#!/bin/bash
source activate phantom_analysis

#define arguments, if an empty string is provided (e.g. for the manually defined ROI), the default argument is interpreted as 'None' in python
input_epi=$1
input_roi=${2:-None}
output_path=$3
desired_slice=${4:-None}
TR=$5

#call the python script
python3.8 phantom_analysis_functions.py $input_epi $input_roi $output_path $desired_slice $TR
