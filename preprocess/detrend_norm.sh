#!/bin/bash

if [[ $# -lt 2 ]]; then
    echo "Outputs normalized and detrended volume as <vol>_dt_norm.nii.gz"
    echo ""
    echo "usage: $0 <vol> <mask>"
    exit
fi

vol=$1
mask=$2

base=${vol//.nii*}
tmp=$(mktemp -dt $(basename $base)) || exit 1

# check da args yo
if [[ ! -e $vol ]]; then
    echo "Error: volume not found: $vol"
    exit 1
fi

if [[ ! -e $mask ]]; then
    echo "Error: Mask not found: $mask"
    exit 1
fi

# num of time points
nvols=$(fslnvols "$vol")
# make linear sequence
echo $(seq 1 $nvols) > ${tmp}_trend

# detrend linear sequence
echo "Detrending"
vol_dt=${base}_dt.nii.gz
fsl_glm -i "$vol" -d ${tmp}_trend --out_data=$vol_dt -m $mask #--demean

# compute mean
echo "Computing mean"
vol_mean=${tmp}_mean.nii.gz
fslmaths "$vol_dt" -Tmean $vol_mean -odt float

# center and normalize variance
echo "Normalizing variance"
vol_dt_norm=${base}_dt_norm.nii.gz 
fslmaths "$vol_dt" -Tstd ${tmp}_std.nii.gz -odt float
fslmaths "$vol_dt" -sub "$vol_mean" -div ${tmp}_std.nii.gz -mas "$mask" "$vol_dt_norm" -odt float 