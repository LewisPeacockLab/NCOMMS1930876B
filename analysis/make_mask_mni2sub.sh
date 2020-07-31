#!/bin/bash

if [[ $# -lt 2 ]]; then
    echo "usage: $0 <fsl_dir> <subject_dir> <mni_mask>"
    exit
fi

FSLDIR=$1
SUBDIR=$2
ST_MASK=$3

###########################################
# environment setup
###########################################
. ${FSLDIR}/etc/fslconf/fsl.sh
PATH=${FSLDIR}/bin:${PATH}
export FSLDIR PATH

###########################################
# directories
###########################################
FUNC_DIR=${SUBDIR}/bold
FUNC_AVG_DIR=${FUNC_DIR}/avg_func_ref

fRI=${FUNC_AVG_DIR}/bold_avg_mcf_brain.nii.gz

if [[ ! -e $fRI ]]; then
    echo "Cannot find file: $fRI. Was the subject preprocessed?"
    exit
fi

MASK=$SUBDIR/masks/$(basename $ST_MASK)

###########################################
#compute inverse transform (standard to MPRAGE)
mni2sub=${FUNC_AVG_DIR}/mni2sub

#set the fRI to MPRAGE mat file 
sub2ref=${FUNC_AVG_DIR}/sub2ref
 
#concatenate both mat and warp files to achieve fRI to standard
echo "Applying MNI->subject transform to mask"
echo "applywarp --ref=${fRI} --in=${ST_MASK} --warp=${mni2sub} --postmat=${sub2ref} --out=${MASK}"
applywarp --ref=${fRI} \
          --in=${ST_MASK} \
          --warp=${mni2sub} \
          --postmat=${sub2ref} \
          --out=${MASK}

# binarized version
# mri_binarize --i $MASK --o ${MASK//.nii*}_bin.nii.gz --min 0.0001
# fslmaths $MASK -thr 1 -bin ${MASK//.nii*}_bin.nii.gz

