#!/bin/bash

# subject's directory
SUBDIR=$1

# reference functional run (by name)
MIDRUN=$2


###########################################
# arguments
###########################################

readopt='getopts $opts opt;rc=$?;[ $rc$opt == 0? ]&&exit 1;[ $rc == 0 ]||{ shift $[OPTIND-1];false; }'
opts=amds:r:

# Enumerating options
while eval $readopt; do
    eval "arg_${opt}=${OPTARG-true}"
done

# default: run both motion and alignment
if [[ -z "$arg_a" ]] && [[ -z "$arg_m" ]] && [[ -z "$arg_d" ]]; then
    arg_a=true
    arg_m=true
    arg_d=true
fi

# required: subject_dir
# required: reference run (if running motion-correction)
if [[ -z "$arg_s" ]] || ( [[ $arg_m ]] && [[ -z "$arg_r" ]] ); then
    echo "usage: $0 [-a,-m,-d] -s <subject_dir> -r <ref_func_run_name>"
    echo "where: -a run alignment"
    echo "       -m run motion correction"
    echo "       -d detrend and normalize"
    echo "   (default is to run all)"
    exit
fi

SUBDIR=$arg_s
MIDRUN=$arg_r

echo "sub dir: $SUBDIR"
echo "ref run: $MIDRUN"

# check if subject directory is real
if [[ ! -d "$SUBDIR" ]]; then
    echo "Error: subject directory does not exist: $SUBDIR"
    exit
fi

##########################################
# current directory 
##########################################

SOURCE="${BASH_SOURCE[0]}"
while [ -h "$SOURCE" ]; do
  # resolve $SOURCE until the file is no longer a symlink
  DIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"
  SOURCE="$(readlink "$SOURCE")"
  # if $SOURCE was a relative symlink, we need to resolve it relative 
  # to the path where the symlink file was located
  [[ $SOURCE != /* ]] && SOURCE="$DIR/$SOURCE"
done
SCRIPT_DIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"


###########################################
# file/directory variables
###########################################

START_DIR=$PWD

MASK_DIR=${SUBDIR}/masks
FUNC_DIR=${SUBDIR}/bold
FMAP_DIR=$SUBDIR/fieldmap
ANAT_DIR=${SUBDIR}/anatomy

FUNC_AVG_DIR=${FUNC_DIR}/avg_func_ref
MIDRUNDIR=${FUNC_DIR}/${MIDRUN}
MIDRUN_BOLD=${MIDRUNDIR}/bold_mid_mcf
pre_fRI=${MIDRUNDIR}/bold_avg_mcf
fRI=${FUNC_AVG_DIR}/bold_avg_mcf_brain
fRI_head=${FUNC_AVG_DIR}/bold_avg_mcf

ref=${ANAT_DIR}/highres001_brain
ref_head=${ANAT_DIR}/highres001

ST_TEMPLATE="$FSL_DIR/data/standard/MNI152_T1_1mm_brain.nii.gz"
ST_TEMPLATE_HEAD="$FSL_DIR/data/standard/MNI152_T1_1mm.nii.gz"
ST_TEMPLATE_MASK="$FSL_DIR/data/standard/MNI152_T1_1mm_brain_mask_dil.nii.gz"

# find all functional scans

SCAN_LIST=$(ls -d ${FUNC_DIR}/* | grep -v $FUNC_AVG_DIR)

###########################################
# make directories
###########################################

cd ${SUBDIR}
mkdir -p ${MASK_DIR}
mkdir -p ${FUNC_AVG_DIR}

###########################################
# motion correction
###########################################

function motion_correction() {
    #Create Middle Run Mean = this is your functional reference image fRI
    # take the bold.nii from the middle run, mcflirt and mean
    echo "movement correction on middle run"
    mcflirt -in ${MIDRUNDIR}/bold.nii.gz -out ${MIDRUN_BOLD} -mats -plots -stages 4 -sinc_final;
    echo "fslmaths ${MIDRUN_BOLD} -Tmean ${pre_fRI}"
    fslmaths ${MIDRUN_BOLD} -Tmean ${pre_fRI} 
    echo "copy middle run average to functional average directory and brain extraction"
    cp ${pre_fRI}.nii.gz ${FUNC_AVG_DIR}
    echo "extract brain from reference image"
    bet ${pre_fRI} ${fRI} -R -m

    # motion correct (coreg) all functionals the to average of the middle run;

    for x in $SCAN_LIST; do
        cd $x; 
        echo "Processing: $x";
        echo " mcflirt -in bold.nii.gz"
        mcflirt -in bold.nii.gz -out bold_mcf -reffile ${fRI_head} -mats -plots
        echo " bet bold_mcf.nii.gz bold_mcf_brain.nii.gz -F"
        bet bold_mcf.nii.gz bold_mcf_brain.nii.gz -F
        echo ' creating temporary mean'
        fslmaths bold_mcf_brain.nii.gz -Tmean tempmean.nii.gz
        echo ' high pass temporal filter'
        fslmaths bold_mcf_brain.nii.gz -bptf 64 -1 -add tempmean.nii.gz bold_mcf_brain_hpass.nii.gz
        rm tempmean.nii.gz;
        cd $START_DIR
    done 
}

###########################################
# alignment (mni -> reference -> t1)
###########################################

function alignment() {
    # extrain brain from anatomical
    if [[ ! -e ${ref_head}.nii.gz ]]; then
        mprage=$(ls $ANAT_DIR/*/mprage.nii.gz | xargs | awk '{print $1}')
        ln -s $mprage ${ref_head}.nii.gz
    fi

    echo "bet ${ref_head} ${ref} -R -B"
    bet ${ref_head} ${ref} -R -B

    #Register your fRI (bold_avg_mcf_brain) to the anatomical scan (MPRAGE)
    echo "epi_reg --epi=${fRI} --t1=${ref_head} --t1brain=${ref}"
    echo "        --out=${FUNC_AVG_DIR}/bold_co_avg_mcf_brain.nii.gz"
    epi_reg --epi=${fRI} \
            --t1=${ref_head} \
            --t1brain=${ref} \
            --out=${FUNC_AVG_DIR}/bold_co_avg_mcf_brain.nii.gz 
    #Go to each functional Run and align with mean image of the middle run (fRI)
    #Obtain the names for each functional run and run mcflirt 

    #you already have the fRI to MPRAGE from epi_reg prior. 
    #Coregister the MPRAGE to the standard (T1 MNI152)

    cd ${FUNC_AVG_DIR}
    echo "creating AFFINE transform for non-linear registration"
    flirt -ref ${ST_TEMPLATE} -in ${ref} -omat affine_sub2mni_4fnirt.mat
    echo "FNIRT in progress, non-linear registration of MPRAGE to MNI152_1mm.nii"
    fnirt --ref=${ST_TEMPLATE_HEAD} \
          --in=${ref_head} \
          --refmask=${ST_TEMPLATE_MASK} \
          --aff=affine_sub2mni_4fnirt.mat \
          --cout=sub2mni \
          --iout=highres_sub2mni_image  

    #compute inverse transform (standard to MPRAGE)
    mni2sub=${FUNC_AVG_DIR}/mni2sub
    #set the fRI to MPRAGE mat file 
    ref2sub=${FUNC_AVG_DIR}/bold_co_avg_mcf_brain.mat
    sub2ref=${FUNC_AVG_DIR}/sub2ref
    sub2mni=${FUNC_AVG_DIR}/sub2mni
     
    echo "compute inverse transform MNI to T1"
    invwarp --ref=${ref_head} --warp=${sub2mni} --out=mni2sub
    #compute inverse transform (T1 to fRI)
    echo "compute inverse transform T1 to fRI"
    convert_xfm -omat sub2ref -inverse ${ref2sub}
}


function detrend() {
    for x in $SCAN_LIST; do
        # for each scan, detrend (using external script)
        cd $x; 
        echo "Processing: $x";
        "$SCRIPT_DIR"/detrend_norm.sh bold_mcf_brain_hpass.nii.gz bold_mcf_brain_mask.nii.gz;
        cd $START_DIR
    done
}

if [[ $arg_m ]]; then
    motion_correction;
fi 

if [[ $arg_a ]]; then
    alignment;
fi

if [[ $arg_d ]]; then
    detrend;
fi

# return to birth
cd $START_DIR
