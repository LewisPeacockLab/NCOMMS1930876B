#!/bin/bash

# nomalize all epis to mni
# $1: cluster 'blanca' | 'local'
# $2: subject's directory
# $3: MNI template directory
# $4: smooth fwhm

###########################################
# arguments
###########################################

args=("$@")

CLUSTER=${args[0]}
SUBDIR=${args[1]}
MASK_DIR=${args[2]}
XFWHM=${args[3]}

# registration
ARGS_REG=false 

echo "... running in: $CLUSTER"
echo "... sub dir: $SUBDIR"

if [ "${CLUSTER}" = "local" ]; then
	export FREESURFER_HOME=/Applications/freesurfer
	source ${FREESURFER_HOME}/FreeSurferEnv.sh
elif [ "${CLUSTER}" = "blanca" ]; then
	export FREESURFER_HOME=/projects/ics/software/freesurfer/5.3.0
	source ${FREESURFER_HOME}/FreeSurferEnv.sh

fi

###########################################
# file/directory variables
###########################################
ANAT_DIR=${SUBDIR}/anatomy
FUNC_DIR=${SUBDIR}/bold
FUNC_AVG_DIR=${FUNC_DIR}/avg_func_ref
fRI=${FUNC_AVG_DIR}/bold_avg_mcf_brain
WARP_PARAM=${SUBDIR}/bold/avg_func_ref/warp_param

mkdir -p ${WARP_PARAM}

# MNI templete
ST_TEMPLATE_HEAD="${MASK_DIR}/MNI152_T1_3mm_brain.nii.gz"
ST_TEMPLATE="${MASK_DIR}/MNI152_T1_3mm.nii.gz"

# collect all parameters for transformation: ref:EPI avg,sub:T1
# functional: 3x3x3mm, T1: .8x.8x.8, mni-mask: 2x2x2
t1=${ANAT_DIR}/highres001
t1_brain=${ANAT_DIR}/highres001_brain
func2struct=${WARP_PARAM}/func2struct.mat
struct2mni=${WARP_PARAM}/struct2mni_3mm
mni2struct=${WARP_PARAM}/mni2struct_3mm
struct2func=${WARP_PARAM}/struct2func
struct2mni_transf=${WARP_PARAM}/affine_sub2mni_4fnirt_3mm.mat

###########################################
# registration
###########################################
# warpres = config (2mm) * invert mm (3mm)
# if [[ $ARGS_REG ]]; then
# 	echo "(+) functional epi (avg) to T1 MPRAGE_brain"
# 	flirt -ref ${t1_brain} -in ${fRI} -dof 7 -omat ${func2struct}
# 	echo "(+) creating AFFINE transform for non-linear registration"
# 	flirt -ref ${ST_TEMPLATE_HEAD} -in ${t1_brain} -omat ${struct2mni_transf}
# 	echo "(+) FNIRT in progress, non-linear registration of MPRAGE to MNI152_3mm"
# 	fnirt --in=${t1} --aff=${struct2mni_transf} --cout=${struct2mni} --config=T1_2_MNI152_2mm --warpres=6,6,6
# fi

###########################################
# transfer bold -> in standard space / smoothing
###########################################

# find all bold folders
SCAN_LIST=$(ls -d ${FUNC_DIR}/clearmem*localizer*)

for x in $SCAN_LIST; do
	cd $x;
	echo "... convert epi to MNI space"
	echo $x

	XNORM=bold_mcf_brain_hpass_dt_mni.nii.gz
	      
	# to standard	  
	applywarp --ref=${ST_TEMPLATE} \
		  	  --in=bold_mcf_brain_hpass_dt.nii.gz \
		  	  --warp=${struct2mni} \
		  	  --premat=${func2struct} \
		  	  --out=${XNORM}

	# smoothing
	echo "... smooting: fwhm = ${XFWHM}"

	XSNORM=$(printf 'bold_mcf_brain_hpass_dt_mni_s%d.nii.gz' "${XFWHM}")
	fslmaths ${XNORM} -s ${XFWHM} ${XSNORM}
done		  

echo "(-) Done normalizing/smoothing localizer epi in standard space"

# find all bold folders
SCAN_LIST=$(ls -d ${FUNC_DIR}/clearmem*study*)

for x in $SCAN_LIST; do
	cd $x;
	echo "... convert epi to MNI space"
	echo $x

	XNORM=bold_mcf_brain_hpass_dt_mni.nii.gz
	
	# to standard	  
	applywarp --ref=${ST_TEMPLATE} \
		  	  --in=bold_mcf_brain_hpass_dt.nii.gz \
		  	  --warp=${struct2mni} \
		  	  --premat=${func2struct} \
		  	  --out=${XNORM}

	# smoothing
	echo "... smooting: fwhm = $XFWHM"

	XSNORM=$(printf 'bold_mcf_brain_hpass_dt_mni_s%d.nii.gz' "$XFWHM")
	fslmaths ${XNORM} -s ${XFWHM} ${XSNORM}
done		  

echo "(-) Done normalizing/smoothing study epi in standard space"
