#!/bin/bash

# $1: cluster 'blanca' | 'local'
# $2: MVPA directory 
# $3: subject's directory
# $4: importance map basename

###########################################
# arguments
###########################################

args=("$@")

CLUSTER=${args[0]}
IMPMAP_DIR=${args[1]}
SUBDIR=${args[2]}
IMPMAP=${args[3]}
MASK_DIR=${args[4]}

echo "(+) sub dir: $SUBDIR"
echo "(+) impmap: $IMPMAP_DIR"
echo "(+) importance map for: $IMPMAP"
echo "(+) running in: $CLUSTER"

if [ "${CLUSTER}" = "local" ]; then
	export FREESURFER_HOME=/Applications/freesurfer
	source ${FREESURFER_HOME}/FreeSurferEnv.sh
elif [ "${CLUSTER}" = "blanca" ]; then
	export FREESURFER_HOME=/projects/ics/software/freesurfer/5.3.0/bin/freesurfer
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

# find all importance map in the decoding folder
cd ${IMPMAP_DIR}
IMPMAP_LIST=$(ls -d ${IMPMAP}*)

###########################################
# reorient
###########################################

for x in $IMPMAP_LIST; do
	
	echo $x
	
	echo "(+) reorient: copy s-form to q-form"
	fslorient -copysform2qform $x
	
done

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
echo "(+) functional epi (avg) to T1 MPRAGE_brain"
flirt -ref ${t1_brain} -in ${fRI} -dof 7 -omat ${func2struct}
echo "(+) creating AFFINE transform for non-linear registration"
flirt -ref ${ST_TEMPLATE_HEAD} -in ${t1_brain} -omat ${struct2mni_transf}
echo "(+) FNIRT in progress, non-linear registration of MPRAGE to MNI152_3mm"
fnirt --in=${t1} --aff=${struct2mni_transf} --cout=${struct2mni} --config=T1_2_MNI152_2mm --warpres=6,6,6

###########################################
# transfer importance map -> in standard space
###########################################

for x in $IMPMAP_LIST; do
	
	echo $x
	
	x_norm_impmap=$(printf 'norm_%s' "$x")
	echo $x_norm_impmap
	
	#concatenate both mat and warp files to achieve importance map to standard
	echo "(+) Applying subject->MNI transform to importance map"		  
	applywarp --ref=${ST_TEMPLATE} \
		  	  --in=$x \
		  	  --warp=${struct2mni} \
		  	  --premat=${func2struct} \
		  	  --out=$x_norm_impmap
done		  

echo "(+) Done normalizing importance map in standard space"
