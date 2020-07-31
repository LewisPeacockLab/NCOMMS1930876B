#!/bin/bash
 
# ./imp_map_correction.sh 10 20 inputDir basename 
# $1: cluster size
# $2: smoothing fwhm x*x*x
# $3: input directory (peak_impmap:top_10_norest)
# $4: file basename

###########################################
# arguments
###########################################
export FREESURFER_HOME=/Applications/freesurfer
source /Applications/freesurfer/FreeSurferEnv.sh

args=("$@")

XCLUSTER=${args[0]}
XFWHM=${args[1]}
XDIR=${args[2]}
XBASE=${args[3]}

echo "(+) cluster correction: $XCLUSTER"
echo "(+) FWHM size: $XFWHM"
echo "(+) impmap dir: $XDIR"
echo "(+) basename: $XBASE"

###########################################
# file/directory variables
###########################################
# find all importance map in the decoding folder
cd ${XDIR}
IMPMAP_LIST=$(ls -d ${XBASE}*)

###########################################
# smoothing/correction
###########################################

for x in $IMPMAP_LIST; do
	
	echo $x
	xname=$(echo "$x" | cut -f 1 -d '.')

    # clusterizing
    echo "cluster correction: $XCLUSTER"
    XCOUT=$(printf '%s/c%d_%s.nii.gz' "$XDIR" "$XCLUSTER" "$XFWHM" "$xname")
    3dmerge -dxyz=1 -1clust 1 $XCLUSTER -prefix $XCOUT $x

    # cluster correction table
    echo "cluster correction table"
    XCOUT_TX=$(printf '%s/out_c%d_%s.txt' "$XDIR" "$XCLUSTER" "$XFWHM" "$xname")
    touch $XCOUT_TX
    3dclust -NN1 $XCLUSTER $XCOUT > $XCOUT_TX

    # smoothing
    echo "smooting: fwhm = $XFWHM"
    XSOUT=$(printf '%s/s%d_c%d_%s.nii.gz' "$XDIR" "$XFWHM" "$XCLUSTER" "$xname")
    3dBlurToFWHM -FWHM $XFWHM -prefix $XSOUT -input $XCOUT -acf $(printf '_s%d_c%d_%s' "$XFWHM" "$XCLUSTER" "$xname")
    # 3dBlurToFWHM -FWHM $XFWHM -automask -prefix $XSOUT -input $x

done

