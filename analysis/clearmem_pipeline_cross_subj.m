function[] = clearmem_pipeline_cross_subj(xcluster, xstep, xmask, xpenalty, subject_id)

% clearmem_pipeline_cross_subj('tacc','mvpa_operation','whole',1),'clearmem_v1_sub001')
% xpenalty: step1: 1-8, step2: 1-10
% xcluster = 'tacc'; xstep = 'mvpa_operation'; xmask = 'whole'; xpenalty = 1; subject_id = 'clearmem_v1_sub001';

args.cluster  = xcluster;
args.xstep    = xstep;
args.analysis = pwd;
args.level    = 'category';

args.selected_category    = 1:3;
args.selected_subcategory = 1:3;
args.class_selecting      = 0;
args.subclass_selecting   = 0;

operation_mvpa_decode = 0;

%% ============= SETUP DIRECTORIES & ARGS
if strcmp(xcluster, 'local')
    curr_dir       = pwd; cd ~
    dirs.home      = pwd; cd(curr_dir);
    dirs.project   = fullfile(dirs.home,'lewpealab_dpbx','Hyojeong','STUDY','Clearmem');
    dirs.fmri_data = fullfile(dirs.project,'imaging_data','utaustin');%'utaustin'
    
    %*************** FSL PATH
    if strcmp(dirs.home, '/Users/hk9643')
        dirs.fsl = '/Applications/fsl';
    elseif strcmp(dirs.home, '/Users/hyojeongkim')
        dirs.fsl = '/usr/local/fsl';
    end
    
elseif strcmp(xcluster, 'blanca')
    dirs.fmri_data = '/work/ics/data/projects/banichlab/studies/wmem/fmri/mvpa/utaustin/data';
    %*************** FSL PATH
    dirs.fsl = '/projects/ics/software/fsl/5.0.10';
    
elseif strcmp(xcluster, 'tacc')
    dirs.work = '/work/03358/hjkim/lonestar/clearmem';%processed data
    dirs.fmri_data = '/work/03358/hjkim/lonestar/clearmem';%raw data
%     dirs.fmri_data = '/scratch/03358/hjkim/clearmem';%raw data
    
    %*************** FSL PATH
    dirs.fsl = '/work/IRC/ls5/opt/apps/fsl-5.0.9';
end

dirs.mni_mask = fullfile(dirs.fmri_data, 'mni-masks');

%% ============= SETUP TOOLBOX & DIRECTORIES
%*************** SPM / MVPA PRINCETON TOOLBOX
if strcmp(xcluster, 'local')
    addpath ~/github/mvpa_custom/
    addpath(genpath('~/github/GroupICATv4.0b/icatb/'))

    %*************** SPM TOOLBOX
    addpath(genpath('/Applications/spm12'))
    
elseif strcmp(xcluster, 'blanca')
    %*************** MVPA PRINCETON TOOLBOX
    addpath('/work/ics/data/projects/banichlab/studies/wmem/fmri/mvpa/utaustin/cu_src/mvpa')
    addpath(genpath('/work/ics/data/projects/banichlab/studies/wmem/fmri/mvpa/utaustin/cu_src/GroupICATv4.0b/icatb'))

elseif strcmp(xcluster, 'tacc')
    %*************** MVPA PRINCETON TOOLBOX
    addpath('/home1/03358/hjkim/mvpa_custom')
    addpath(genpath('/corral-repl/utexas/lewpealab/archived/software/GroupICATv4.0b/icatb/'))
    addpath(genpath('/corral-repl/utexas/lewpealab/archived/software/spm12'))
end

mvpa_add_paths;
spm('Defaults','fMRI'); spm_jobman('initcfg');

%*************** HOMEMADE CODES
addpath ./homemade/
addpath ./spm_custom/

%% ============= SETUP SUBJECTS
if nargin < 5
    clear w_list
    xlist = load('subj_lists.mat');
    
    for i = 1:length(xlist.list)
        subject_list(i).name = sprintf('clearmem_v1_sub%0.3d', xlist.list(i)); %#ok<*AGROW>
        w_list{i,1} = subject_list(i).name;
    end
    fid = fopen('subject_lists.txt','w');
    fprintf(fid,'%s\n', w_list{:});
    fclose(fid);
    
    args.n_sub           = length(subject_list);
    args.g_sub           = 1:args.n_sub;
    args.filtered_subs   = args.g_sub;%
    xrefer_subj          = 'clearmem_v1_sub004';
    
else
    subject_list(1).name = subject_id;
    args.n_sub           = length(subject_list);
    args.g_sub           = 1;
    args.filtered_subs   = args.g_sub;
    
    xrefer_subj          = 'clearmem_v1_sub004';
end

disp(subject_list')
args.subject_list = subject_list;

%% ============= CLASSIFIER ARGUMENTS
%*************** REGRESOR PARAMETER
args.shift_TRs = 10;%str2double(xshift);% 10, 14
args.xtrim     = 0;

%*************** COLOR SETUP
args.cond_color{1} = [189, 0, 6]/255;
args.cond_color{2} = [253, 81, 0]/255;
args.cond_color{3} = [255, 155, 36]/255;
args.cond_color{4} = [0, 149, 137]/255;
args.cond_color{5} = [0, 134, 196]/255;

args.base_color    = [144, 144, 144]/255;

%*************** for func > mni (3x3x3 mm)
xsample_vox    = '3';
args.mni       = fullfile(dirs.mni_mask, sprintf('MNI152_T1_%smm.nii.gz', xsample_vox));
args.mni_brain = fullfile(dirs.mni_mask, sprintf('MNI152_T1_%smm_brain.nii.gz', xsample_vox));
args.mni_mask  = fullfile(dirs.mni_mask, sprintf('MNI152_T1_%smm_brain_mask.nii.gz', xsample_vox));
args.mni_grey  = fullfile(dirs.mni_mask, sprintf('MNI152_T1_%smm_brain_grey_bin0.5.nii.gz', xsample_vox));
  
if ~exist(args.mni_grey, 'file')
    setenv('FSLDIR', dirs.fsl);  % this to tell where FSL folder is
    setenv('FSLOUTPUTTYPE', 'NIFTI_GZ'); % this to tell what the output type would be
    
    xref_vol   = sprintf('%s/data/standard/MNI152_T1_1mm.nii.gz', dirs.fsl); 
    xref_brain = sprintf('%s/data/standard/MNI152_T1_1mm_brain.nii.gz', dirs.fsl);
    xref_mask  = sprintf('%s/data/standard/MNI152_T1_1mm_brain_mask.nii.gz', dirs.fsl);
    xref_grey  = sprintf('%s/MNI152_T1_1mm_brain_grey.nii.gz', dirs.mni_mask);
    
    system(sprintf('%s/bin/flirt -in %s -ref %s -out %s -applyisoxfm %d -usesqform', ...
        dirs.fsl, xref_vol, xref_vol, args.mni, str2double(xsample_vox)));
    
    system(sprintf('%s/bin/flirt -in %s -ref %s -out %s -applyisoxfm %d -usesqform', ...
        dirs.fsl, xref_brain, xref_brain, args.mni_brain, str2double(xsample_vox)));
    
    system(sprintf('%s/bin/flirt -in %s -ref %s -out %s -applyisoxfm %d -usesqform', ...
        dirs.fsl, xref_mask, xref_brain, args.mni_mask, str2double(xsample_vox)));
    
    system(sprintf('%s/bin/flirt -in %s -ref %s -out %s -applyisoxfm %d -usesqform', ...
        dirs.fsl, xref_grey, xref_brain, args.mni_grey, str2double(xsample_vox)));
    
end

switch(xmask)
    case'whole'
        args.wholebrain = 1;
        args.mask_name  = sprintf('MNI152_T1_%smm_brain_grey', xsample_vox);
        args.whole_mask = sprintf('MNI152_T1_%smm_brain_grey', xsample_vox);
end

%% ============= SETUP PARAMETERS / ARGS STRUCTURE
[~,user] = system('whoami'); [~,host] = system('hostname');

%*************** SCRIPT PARAMETER
args.script        = mfilename;
args.user          = strtrim(user);
args.hostname      = strtrim(host);
args.date.start    = sprintf('%s_%s', datetime, args.hostname);

%*************** DATA PARAMETER
args.experiment    = 'clearmem';
args.midrun        = 'avg_func_ref';
args.phase_name    = {'localizer','study','operation'};
args.epiext        = 'nii';
args.epi_name      = 'bold_mcf_brain_hpass_dt_mni';%'bold_mcf_brain_hpass_dt';%'bold_mcf_brain';%

%*************** FEATURE SELECTION + CLASSIFIER PARAMETER
args.featSel       = 1;
args.featSelThresh = '0.05';%'0.05';%'0.05';%
args.featVox       = 0;%feature selection based on top voxel                    
args.fsVosNum      = 1000;%number of voxels selected                            

args.penalty       = [0 0.01 0.1 1.0 10 100 1000 10000];%[0 1 10:10:70 100:100:1000 10000];
args.it_penalty    = xpenalty;%for step1:1-8, for step2:1-10
args.n_pen         = 10;%step2 narrow penalty numbers
args.classifier    = 'L2logreg';%'L2logreg' & 'ridge' 'logreg'

%*************** operation
args.four_oper_regress = 1;% 1_using 4 operations (remove all repSubcate)
args.reset_regs{1} = 1;%run new penalty
args.reset_regs{2} = 1;%run new penalty
args.reset_regs{3} = 1;%run new penalty
args.reset_6tr{3}  = 0;%run new regressor-based MVPA (only operation 6TRs) + new penalty
args.stw{3}        = 0;%sliding-time-window

%% ============= SETUP ANALYSES PROCEDURE

switch(xstep)
    case 'mvpa_operation'
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% WORKING MEMORY OPERATION DECODING
        % step 1: MNI normalization + smoothing in subject-level
        % step 2: concatenate data + feature selection
        % step 3: penalty check: cross-subject classification
        % step 4: classification (max)
        % step 5: parse the results
        
        operation_mvpa_decode = [1 0 0 0 0];
        
        args.cs_type = 'anatomy_align';%'localizer_align'|'operation_align'
        
        % step 1: normalization in 1st level
        args.subj_norm   = 0;%normalization to MNI for each subject
        args.subj_smooth = 0;%in mm for smoothing after MNI normalization
        
        % step 3: penalty check
        args.penalty_step = 1;%1_broad/2_narrow
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%============= WORKING MEMORY OPERATION DECODING (TRAIN:STUDY + TEST:STUDY)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if sum(operation_mvpa_decode)
    xph                     = 3;%1_localizer, 2_study, 3_operation
    args.xphase             = xph;
    args.rest               = 'norest';%xrest;
    
    %================= study arguments
    args.train_phase        = 'study';
    args.train_regress_type = 'shift';
    args.test_phase         = 'study';
    args.test_regress_type  = 'shift';
    
    args.cross_valid        = 1;%1 is default: n-1 cross-validation 
    args.cross_decoding     = 0;%1 cross validation + decoding on all testing timepoints (non-shifted)
    
    %*************** ANALYSIS PARAMETERS
    args.regress_type       = args.train_regress_type;
    
    %*************** MVPA_CLASSIFICATION
    if operation_mvpa_decode(1)
        %*************** first level
        for xsub = args.filtered_subs
            %*************** setup subject & directories
            args.subject_num = xsub;
            args.subject_id  = subject_list(xsub).name;
            dirs             = setup_directory(dirs, args);
            
            %*************** mkdir
%             if ~isdir(dirs.mvpa.cs.home), mkdir(dirs.mvpa.cs.home);end
%             if ~isdir(dirs.mvpa.cs.subj), mkdir(dirs.mvpa.cs.subj); end
            
            if args.subj_norm
                %*************** set FSL environment
                setenv('FSLDIR', dirs.fsl);  % this to tell where FSL folder is
                setenv('FSLOUTPUTTYPE', 'NIFTI_GZ'); % this to tell what the output type would be
                
                %*************** CONVERT SUBJECT TO STANDARD (MNI)
                % # $1: cluster 'blanca' | 'local'
                % # $2: subject's directory
                % # $3: MNI template directory
                % # $4: smooth fwhm
                
                fprintf('(+) normalize subject bold images: %s\n', datestr(now, 0));
                
                system(sprintf('./MNI_inversetrans_3mm.sh %s %s %s %s', ...
                    args.cluster, dirs.subj_home, ...
                    fullfile(dirs.fmri_data,'mni-masks'), num2str(args.subj_smooth)))
            end
            
            %*************** operation_mvpa_decode
            if operation_mvpa_decode(1), clearmem_operation_cs_subj(args, dirs); end
        end
    end
    
    %********************************
    % 2nd level group analyses
    %********************************
    if sum(operation_mvpa_decode(2:5))
        
        %*************** setup subject & directories
        args.subject_id = xrefer_subj;
        dirs            = setup_directory(dirs, args);
        
        %*************** mkdir
        if ~isdir(dirs.mvpa.cs.home), mkdir(dirs.mvpa.cs.home);end
        if ~isdir(dirs.mvpa.cs.scratch), mkdir(dirs.mvpa.cs.scratch);end
        if ~isdir(dirs.mvpa.cs.out), mkdir(dirs.mvpa.cs.out);end
        if ~isdir(dirs.mvpa.cs.penalty), mkdir(dirs.mvpa.cs.penalty); end
        
        %*************** operation_mvpa_decode
        if operation_mvpa_decode(2), clearmem_operation_cs_01(args, dirs); end
        
    end
end

%% ============= RM PATH
%*************** SPM / MVPA PRINCETON TOOLBOX
if strcmp(xcluster, 'local')
    rmpath(genpath('~/github/mvpa_custom/'))
    rmpath(genpath('~/github/GroupICATv4.0b/icatb/'))
elseif strcmp(xcluster, 'blanca')
    %*************** MVPA PRINCETON TOOLBOX
    rmpath(genpath('/work/ics/data/projects/banichlab/studies/wmem/fmri/mvpa/utaustin/cu_src/mvpa'))
    rmpath(genpath('/work/ics/data/projects/banichlab/studies/wmem/fmri/mvpa/utaustin/cu_src/GroupICATv4.0b/icatb'))
elseif strcmp(xcluster, 'tacc')
    rmpath('/home1/03358/hjkim/mvpa_custom')
    rmpath(genpath('/corral-repl/utexas/lewpealab/archived/software/GroupICATv4.0b/icatb/'))
end

%*************** HOMEMADE CODES
rmpath(fullfile(pwd,'homemade'))
rmpath(fullfile(pwd,'spm_custom'))

