function[] = clearmem_pipeline(xcluster, xstep, xlevel, xtrim, xrest, xmask, subject_id)

%**********************************************************
    %%%%% all variables should be in string
    %%%%% subject_id: if input, individual analysis; else, g_sub
    %%%%% xstep:
    %%%%%   1. 'presetup': motion-correction check, refine mask
    %%%%%   2. 'mvpa_localizer': localizer mvpa
    %%%%%   3. 'mvpa_WM': working memory mvpa
    %%%%%   4. 'mvpa_operation': operation
    %%%%%   5. 'rsa_WM': working memory rsa
    %%%%%   6. 'mvpa_operation_stw': sliding-time-window operation
    %%%%% xcluster: 'local' / 'blanca' (uc)
    %%%%% xlevel: 'category' / 'subcategory'
    %%%%% xtrim: '0' / '1'
    %%%%% xrest: 'rest' / 'norest'
    %%%%% xmask: 'whole' / 'vvs'
    %%%%% 'category' first: presetup & operation decoding are turned off for 'subcategory'
    
    %**** testing args
    % xcluster='local'; xstep='mvpa_localizer'; xlevel='category'; xtrim='0'; xrest = 'norest'; xmask = 'vvs'; subject_id='clearmem_v1_sub001';
    % xcluster='local'; xstep='mvpa_localizer'; xlevel='subcategory'; xtrim='0'; xrest = 'norest'; xmask = 'whole'; subject_id='clearmem_v1_sub001'
    % xcluster='local'; xstep='mvpa_operation'; xlevel='category'; xtrim='0'; xrest = 'norest'; xmask = 'whole'; subject_id='clearmem_v1_sub080';
    % xcluster='local'; xstep='mvpa_WM'; xlevel='category'; xtrim='0'; xrest = 'norest'; xmask='vvs'; subject_id='clearmem_v1_sub001';
    % xcluster='local'; xstep='mvpa_WM'; xlevel='subcategory'; xtrim='0'; xrest = 'norest'; xmask='whole'; subject_id='clearmem_v1_sub001';
    % xcluster='local'; xstep='rsa_WM'; xlevel='category'; xtrim='0'; xrest = 'norest'; xmask='whole'; subject_id='clearmem_v1_sub032';
    
    % clearmem_pipeline('local', 'mvpa_localizer', 'category', '0', 'norest', 'vvs','clearmem_v1_sub001')
    % clearmem_pipeline('local', 'mvpa_localizer', 'subcategory', '0', 'norest', 'whole')
    % clearmem_pipeline('local', 'mvpa_WM', 'subcategory', '0', 'norest', 'whole')
    % clearmem_pipeline('local', 'mvpa_WM', 'category', '0', 'norest', 'vvs')
    % clearmem_pipeline('local', 'rsa_WM', 'category', '0', 'norest', 'vvs')
    % clearmem_pipeline('local', 'rsa_WM', 'item', '0', 'norest', 'vvs','clearmem_v1_sub001')
    % clearmem_pipeline('local', 'rsa_WM', 'item', '0', 'norest', 'hc','clearmem_v1_sub001')
    % clearmem_pipeline('local', 'mvpa_operation', 'category', '0', 'norest', 'whole','clearmem_v1_sub001')
    % clearmem_pipeline('local', 'mvpa_operation_stw', 'category', '0', 'norest', 'whole'),'clearmem_v1_sub001')
    
    % xcluster='blanca'; xstep='presetup'; xlevel='category'; xtrim='0'; xrest ='norest'; xmask ='vvs'; subject_id='clearmem_v1_sub050';
    % xcluster='blanca'; xstep='mvpa_localizer'; xlevel='category'; xtrim='0'; xrest ='norest'; xmask ='vvs'; subject_id='clearmem_v1_sub004';
    % xcluster='blanca'; xstep='mvpa_localizer'; xlevel='subcategory'; xtrim='0'; xrest ='norest'; xmask ='whole'; subject_id='clearmem_v1_sub044';
    % xcluster='blanca'; xstep='mvpa_operation'; xlevel='category'; xtrim='0'; xrest = 'norest'; xmask = 'whole'; subject_id='clearmem_v1_sub001';
    % xcluster='blanca'; xstep='mvpa_WM'; xlevel='category'; xtrim='0'; xrest = 'norest'; xmask='vvs'; subject_id='clearmem_v1_sub001';
    % xcluster='blanca'; xstep='mvpa_WM'; xlevel='subcategory'; xtrim='0'; xrest = 'norest'; xmask='whole'; subject_id='clearmem_v1_sub084';
    % xcluster='blanca'; xstep='rsa_WM'; xlevel='category'; xtrim='0'; xrest = 'norest'; xmask='whole'; subject_id='clearmem_v1_sub001';
    
    % clearmem_pipeline('blanca', 'mvpa_localizer', 'category', '0', 'norest', 'vvs','clearmem_v1_sub001')
    % clearmem_pipeline('blanca', 'mvpa_localizer', 'subcategory', '0', 'norest', 'whole','clearmem_v1_sub004')
    % clearmem_pipeline('blanca', 'rsa_WM', 'category', '0', 'norest','vvs','clearmem_v1_sub001')
    % clearmem_pipeline('blanca', 'rsa_WM', 'category', '0', 'norest', 'hc')
    % clearmem_pipeline('blanca', 'rsa_WM', 'item', '0', 'norest','vvs','clearmem_v1_sub001')
    % clearmem_pipeline('blanca', 'rsa_WM', 'item', '0', 'norest','vvs')
    % clearmem_pipeline('blanca', 'mvpa_WM', 'category', '0', 'norest', 'vvs'),'clearmem_v1_sub001')
    % clearmem_pipeline('blanca', 'mvpa_WM', 'subcategory', '0', 'norest', 'whole'), 'clearmem_v1_sub077')
    % clearmem_pipeline('blanca', 'mvpa_operation', 'category', '0', 'norest', 'whole'),'clearmem_v1_sub001')
    % clearmem_pipeline('blanca', 'presetup', 'category', '0', 'norest', 'vvs')
    % clearmem_pipeline('blanca', 'mvpa_operation_stw', 'category', '0', 'norest', 'whole'),'clearmem_v1_sub001')
%**********************************************************

presetup = 0; localizer_classification = 0; mvpa_decode = 0; 
operation_mvpa_decode = 0; rsa_decode = 0; rsa_spm_result = 0;
operation_mvpa_stw = 0;

%*************** PICK CATEGORY
args.class_selecting   = 0;% select category
args.selected_category = [1 3];%1:3;

args.subclass_selecting   = 0;
args.selected_subcategory = 1:3;%1:3;

if ~(args.class_selecting); args.selected_category = 1:3; end
if ~(args.subclass_selecting); args.selected_subcategory = 1:3; end

args.four_oper_regress = 1;% 1_using 4 operations (remove all repSubcate)

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
    addpath(genpath('/projects/ics/software/spm/spm12_v6470'))
    mvpa_add_paths;
elseif strcmp(xcluster, 'tacc')
    %*************** MVPA PRINCETON TOOLBOX
    addpath('/home1/03358/hjkim/mvpa_custom')
    addpath('/corral-repl/utexas/lewpealab/archived/software/GroupICATv4.0b/icatb/icatb_io_data_functions')
    addpath('/corral-repl/utexas/lewpealab/archived/software/GroupICATv4.0b/icatb')
    addpath('/corral-repl/utexas/lewpealab/archived/software/GroupICATv4.0b/icatb/icatb_spm8_files/icatb_nifti')
    mvpa_add_paths;
end

%*************** SPM / MVPA PRINCETON TOOLBOX
spm('Defaults','fMRI'); spm_jobman('initcfg');
mvpa_add_paths;

%*************** HOMEMADE CODES
addpath ./homemade/
addpath ./spm_custom/

args.xcluster = xcluster;
args.analysis = pwd;

%% ============= SETUP DIRECTORIES & ARGS
if strcmp(xcluster, 'local')
    curr_dir       = pwd; cd ~
    dirs.home      = pwd; cd(curr_dir);
    dirs.project   = fullfile(dirs.home,'Clearmem');
    dirs.fmri_data = fullfile(dirs.project,'imaging_data','utaustin');%'utaustin'
    
elseif strcmp(xcluster, 'blanca')
    dirs.fmri_data =  '/work/ics/data/projects/banichlab/studies/wmem/fmri/mvpa/utaustin/data';
elseif strcmp(xcluster, 'tacc')
    dirs.fmri_data =  '/work/03358/hjkim/lonestar/clearmem';
end

%% ============= SETUP SUBJECTS
if nargin < 7
    subject_lists         = dir(fullfile(dirs.fmri_data,'clearmem_v1_sub*'));
    args.n_sub           = length(subject_lists);
    args.g_sub           = 1:args.n_sub;
    
    excluded_subs        = unique([2, 7, 8, 16, 18]); %[2, 9, 10, 21, 23];% excluded subs
    auc_cutoff_subs      = [4, 10, 14, 24, 40];%[4, 10, 14, 24, 32, 55];%[4, 10, 14, 20, 24, 32, 35, 55];%n=8 [10, 14, 20, 24, 32, 55]//4(005), 35(052)
    
    args.g_sub           = args.g_sub(~ismember(args.g_sub, excluded_subs));
    args.filtered_subs   = args.g_sub(~ismember(args.g_sub, auc_cutoff_subs));
    xrefer_subj          = 'clearmem_v1_sub004';
    
    filtered_subs        = args.filtered_subs;
    save(fullfile(dirs.fmri_data, 'filtered_subjects.mat'), 'filtered_subs');
else
    subject_lists(1).name = subject_id;
    args.n_sub           = length(subject_lists);
    args.g_sub           = 1;
    args.filtered_subs   = args.g_sub;
    
    xrefer_subj          = 'clearmem_v1_sub004';
end

args.subject_list = subject_lists;

%% ============= CLASSIFIER ARGUMENTS
switch(xmask)
    case'whole'
        args.wholebrain = 1;
        args.mask_name = 'highres001_brain_grey_3mm_bin0.2';
    case 'vvs'
        args.wholebrain = 0;  
        args.mask_name  = 'mask_vvps_temporal_fined';
    case 'hc'
        args.wholebrain = 0;  
        args.mask_name  = 'hippocampus';
    case 'ppa'
        args.wholebrain = 0;  
        args.mask_name  = 'PPAsphere_10_bin';
end

args.rest              = xrest;
args.trim_trs          = str2double(xtrim); %0_non_trimmed, 1_trimmed
if args.trim_trs, args.xtrim = 0;
else              args.xtrim = 10; end% # of TR trimmed 
args.cluster           = xcluster;
args.level             = xlevel;%'category';%'subcategory';%

%*************** REGRESOR PARAMETER
args.LSS               = 0;
args.shift_TRs         = 10;%str2double(xshift);% 10, 14

%*************** COLOR SETUP
args.cond_color{1} = [189, 0, 6]/255;
args.cond_color{2} = [253, 81, 0]/255;
args.cond_color{3} = [255, 155, 36]/255;
args.cond_color{4} = [0, 149, 137]/255;
args.cond_color{5} = [0, 134, 196]/255;

args.cate_color{1} = [243, 116, 56]/255;%face
args.cate_color{2} = [185, 35, 83]/255;%fruit
args.cate_color{3} = [25, 97, 105]/255;%scene

args.base_color    = [144, 144, 144]/255;
args.onset_color   = 'r';

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
args.whole_mask    = 'highres001_brain_grey_3mm_bin0.2';%'bold_avg_mcf_brain_mask';%epi masks
args.epi_name      = 'bold_mcf_brain_hpass_dt';%'bold_mcf_brain_hpass_dt';%'bold_mcf_brain';%
args.spm_epi_name  = 'bold_mcf_brain';%
    
%*************** FEATURE SELECTION + CLASSIFIER PARAMETER
args.featSel       = 1;
args.featSelThresh = '0.05';%'0.05';%'0.05';%
args.featVox       = 0;%feature selection based on top voxel                    
args.fsVosNum      = 1000;%number of voxels selected                            

args.penalty       = [0 0.01 0.1 1.0 10 100 1000 10000];%[0 1 10:10:70 100:100:1000 10000];
args.n_pen         = 10;%step2 narrow penalty numbers
args.classifier    = 'L2logreg';%'L2logreg' & 'ridge' 'logreg'
args.repetitions   = 1;%n_iterations

%*************** IMPORTANCE MAP
args.imp_type      = 'mcduff';
args.peak_thresh   = 3;%10;%2;%0.1;% top 2sd
if strcmp(xcluster, 'local')
    dirs.fsl = '/usr/local/fsl';
elseif strcmp(xcluster, 'blanca')
    dirs.fsl = '/projects/ics/software/fsl/5.0.10';
end

%*************** for func > mni (3x3x3 mm)
xsample_vox    = '3';
args.mni       = fullfile(dirs.fmri_data, 'mni-masks', sprintf('MNI152_T1_%smm.nii.gz', xsample_vox));
args.mni_brain = fullfile(dirs.fmri_data, 'mni-masks', sprintf('MNI152_T1_%smm_brain.nii.gz', xsample_vox));

if ~exist(args.mni_brain, 'file')
    setenv('FSLDIR', dirs.fsl);  % this to tell where FSL folder is
    setenv('FSLOUTPUTTYPE', 'NIFTI_GZ'); % this to tell what the output type would be
    
    xref_vol       = sprintf('%s/data/standard/MNI152_T1_1mm.nii.gz', dirs.fsl); 
    xref_brain     = sprintf('%s/data/standard/MNI152_T1_1mm_brain.nii.gz', dirs.fsl);
    
    system(sprintf('%s/bin/flirt -in %s -ref %s -out %s -applyisoxfm %d -usesqform', ...
        dirs.fsl, xref_vol, xref_vol, args.mni, str2double(xsample_vox)));
    
    system(sprintf('%s/bin/flirt -in %s -ref %s -out %s -applyisoxfm %d -usesqform', ...
        dirs.fsl, xref_brain, xref_brain, args.mni_brain, str2double(xsample_vox)));
    
end

%*************** TIMECOURSE
args.tc_tr         = 38;
args.dur_sync      = 16;
args.tc_tr_disp    = 30;%for stats/display
args.dur_sync_disp = 12;
args.n_tr_blks     = 3;%6
args.alpha         = 0.05;

args.item_peak     = [11 16; 22 27];%for replace: 1_before, 2_after switch
args.baseline_trs  = 6;
args.ttest_trs     = 7:21;%24
args.anova_trs     = 7:21;%1:args.tc_tr_disp;

%*************** operation
args.reset_regs{1} = 1;%run new penalty
args.reset_regs{2} = 1;%run new penalty
args.reset_regs{3} = 1;%run new penalty
args.reset_6tr{3}  = 0;%run new regressor-based MVPA (only operation 6TRs) + new penalty
args.stw{3}        = 0;%sliding-time-window
args.xstep         = xstep;

%% ============= SETUP ANALYSES PROCEDURE
% NOTE: args.load_pattern{phase}: EPI patterns are loaded and saved in the scratch
% so, if the pattern has been loaded before put 0 to just load the
% patterns from subject structure

switch(xstep)
    case 'presetup'
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PRE-SETUP
        presetup                   = 1;
        cutting_off_mask_outvoxels = 0;%#ok<*NASGU> %checking the mask if it's within mean-whole-brain
        motion_plotting            = 0;%creating motion correction plots
        creat_param                = 0;
        motion_confound            = 1;
        
    case 'mvpa_localizer'
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LOCALIZER
        % step 1: read EPI patterns + regressors/selectors + feature selection
        % step 2: classification (with all penalties if args.penalty_check = 1)
        %         if args.penalty_check = 0, load maximum penalty
        % step 3: parse the mvpa outcome / AUC 1st
        % step 4: AUC 2nd
        % step 5: importance map 1st + 2nd
        localizer_classification   = [0 0 0 1 0];%[1 1 1 0 0];
        
        % step 1
        args.load_pattern{1}       = 1;% load the pattern at the first time
        args.read_mat{1}           = 0;% read pattern from whole pattern mat
        
        % step 2
        args.fixed_penalty{1}      = 0;
        if args.fixed_penalty{1}, args.xpenalty = 50; end
        args.penalty_check{1}      = 1;% 0_load max penalty
        args.penalty_step{1}       = 1;% 1_taking two steps penalty check
        args.run_selecting         = 0;% 1_run only selected run: args.selected_run
        
        % step 3
        args.permutation{1}        = 0;% random permutation for baseline
        args.n_iteration           = 1000;
        
        % step 4: auc analyses
        args.mvpa_out{1}           = 1;% 1_save mvpa results in group structure
        args.grp_auc{1}            = 1;% 1_save 1st-level AUC into group structure
        
        % step 5: importance map
        args.subj_impmap{1}        = 1;% 1_run 1st-level importance map
        args.group_mat{1}          = 1;% 1_save 1st-level importance map into group structure
        args.group_impmap{1}       = 1;% 1_2nd-level
        args.group_impmap_diff{1}  = 0;% 1_differences of subcategory_only for subcategory level
        
        if strcmp(args.level, 'category'), args.group_impmap_diff{1} = 0; end
       
    case 'mvpa_WM'
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% WORKING MEMORY DECODING
        % step 1: train: EPI patterns + regressors + feature selection
        %         test: EPI patterns under the feature-selected mask with training data
        %               + regressors
        %         decoding setup: combine train & test: EPI + regressors/selectors
        % step 2: load penalty from localizer classifier + classification
        % step 3: parse the mvpa outcome + 1st level analysis
        % step 4: 2nd level analysis + timecourse
        % step 5: AUC 2nd
        mvpa_decode                    = [0 0 0 1 0];
        
        args.load_pattern{2}           = 1;% load study patterns under ROI
        args.read_mat{2}               = 0;% read pattern from operation pattern mat
        
        % step 3
        args.parse_mvpa{2}             = 1;% parse mvpa again
        args.add_mvpa{2}               = 1;% add additional parse(add to mvpaout)
        
        if sum(mvpa_decode(4:5))
            % step 4
            args.group_mvpa{2}         = 1;% 1_save 1st-level parsed mvpaout into group structure
            args.filtered_subs         = args.filtered_subs;%args.g_sub;
            args.parse_timecourse      = 0;
            args.parse_timecourse_cate = 1;
            
            %*************** BASELINE CORRECTION
            args.bin_windows           = [11 16; 22 27];%for replace: 1_before, 2_after switch
            
            % step 5: auc analyses (note: higher than Matlab R2014b)
            args.mvpa_out{2}           = 1;% 1_save mvpa results in group structure
            args.grp_auc{2}            = 1;% 1_save 1st-level AUC into group structure
        end
        
        args.load_mvpa_cat{2}          = 1;% for separated subcategory mvpa
        
    case 'mvpa_operation'
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% WORKING MEMORY OPERATION DECODING
        % step 1: read EPI patterns + regressors/selectors + feature selection
        %         decoding setup: decode all timepoints for the testing run
        % step 2: classification (penalty: 50)
        % step 3: parse the mvpa outcome/ AUC 1st level
        % step 4: 2nd level analysis
        % step 5: AUC/ROC 2nd
        % step 6: importance map 1st + 2nd
        
        operation_mvpa_decode          = [0 0 0 0 1 0 0];%[1 1 1 0 0 0];
        
        % step 1
        args.load_pattern{3}           = 0;% load study patterns under ROI (whole brain)
        args.read_mat{3}               = 0;
        
        % step 2
        args.fixed_penalty{3}          = 0;
        if args.fixed_penalty{3}, args.xpenalty = 50; end
        args.penalty_check{3}          = 1;% 0_load max penalty
        args.penalty_step{3}           = 1;% 1_taking two steps penalty check
        % step 3
        args.parse_mvpa{3}             = 1;% parse mvpa again
        args.permutation{3}            = 0;% random permutation for baseline
        args.n_iteration               = 1000;
        
        if sum(operation_mvpa_decode(4:7))
            % step 4
            args.group_mvpa{3}         = 1;% 1_save 1st-level parsed mvpaout into group structure
            args.mvpa_results{3}       = 1;% 1_save 1st-level classifier results into group structure
            args.filtered_subs         = args.filtered_subs;%args.g_sub;
            
            % step 5: auc analyses (note: higher than Matlab R2014b)
            args.grp_auc{3}            = 1;% 1_save 1st-level AUC into group structure
            args.se_filtered           = 0;% cutoff
            args.se_cutoff             = 10;
         
            % step 6: importance map
            args.subj_impmap{3}        = 1;% 1_run 1st-level importance map
            args.group_mat{3}          = 1;% 1_save 1st-level importance map into group structure
            args.group_impmap{3}       = 1;% 1_2nd-level
            args.group_impmap_diff{3}  = 0;% 1_differences of subcategory_only for subcategory level
            args.impmap_1st_reset      = 0;%reset 1st level maps
            
            % step 7: cross-subject
            args.subj_norm   = 1;%normalization to MNI for each subject
            args.subj_smooth = 2;%in mm for smoothing after MNI normalization
        end
        
        %*************** turn off presetup & operation decoding
        if strcmp(xlevel, 'subcategory')
            operation_mvpa_decode = [0 0 0 0 0 0 0];
        end
    case 'rsa_WM'
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% WORKING MEMORY DECODING
        % step 1. ROI extraction (GLM): localizer / 128 hf applied
        % step 2. create template by averaging the extracted pattern
        % step 3: RSA timecourse: test on study: 1st + 2nd level
        % step 4: item RSA bootstrap: N = N+1
        
        %*************** RESULT PARAMETERS
        rsa_decode               = [0 0 1 0];   
        rsa_spm_result           = 0;
        args.new_roi             = 0;
        
        % step 1
        args.spm_p_thresh        = 0.05;% t threshold of uncorrected p < 0.05
        if strcmp(args.mask_name, 'hippocampus') || strcmp(args.mask_name, 'PPAsphere_10_bin')
            args.spm_v_extent    = 0;% extent voxel
        else
            args.spm_v_extent    = 10;% extent voxel
        end
        args.spm_thresh_type     = 'none';%uncorrected 'FWE'; 'FDR';
        args.item_within         = 0;%contrast within category
        
        % step 2
        args.rsa_mask            = args.mask_name;
        args.blk                 = 1;%category, 1_block TR window, 0_stimulus on window (2TR)
        
        args.load_pattern_rsa{1} = 1;% localizer
        args.load_pattern_rsa{2} = 1;% study
        args.rsa_pattern         = [1 1];
        
        % step 3
        args.grp_pattern         = 0;%saving individual pattern into grp structure
        args.additional_rsa      = 0;
        args.percept_win         = 11:16;%peak for presented item in study
        args.bin_trs{1}          = 13:18;%n-targ
        args.bin_trs{2}          = 11:16;%in synched (from 22trs)
        args.bin_trs{3}          = 5:10;%baseline in synched (from 22trs)
        
        % step 4
        args.boot_strap_rsa      = 0;
        args.n_iteration         = 1000;
        args.n_sample            = 15;% sample trials
    
    case 'mvpa_operation_stw'
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SLIDING-TIME-WINDOW WM OPERATION DECODING
        % step 1: load EPI patterns + regressors/selectors
        %         + feature selection + classification (penalty: 0.5)
        % step 2: parse the mvpa outcome
        % step 3: 2nd level analysis
        
        operation_mvpa_stw       = [0 0 1];
        
        args.stw{3}              = 1;
        args.reset_regs{3}       = 0; 
        args.reset_6tr{3}        = 0;
        args.four_oper_regress   = 0;
        
        % step 1
        args.shift_TRs           = 0;
        args.stw_win             = 5;%time-window for stw
        args.xpenalty            = 50; 
        args.classification      = 0;%run only cross-validation
                
        if sum(operation_mvpa_stw(3))
            % step 3
            args.group_mvpa{3}   = 1;% 1_save 1st-level parsed mvpaout into group structure
            args.filtered_subs   = args.filtered_subs;%args.g_sub;
        end
end

%% ============= PRE SETUP
if presetup
    for xsub = args.g_sub
        args.xphase     = 1:2;
        args.wholebrain = 0;
        args.mask_name  = 'mask_vvps_temporal';
        
        %*************** setup subject & directories
        args.subject_id = subject_lists(xsub).name;
        dirs            = setup_directory(dirs, args);
        
        %*************** mkdir
        if ~isdir(dirs.motion), mkdir(dirs.motion); end
            
        %*************** creating fined mask
        if cutting_off_mask_outvoxels, create_fined_mask(args, dirs); end %#ok<*UNRCH>
        
        %*************** creating motion plot
        if motion_plotting, create_motion_plot(args, dirs); end
        
        %*************** creating parameter index
        if creat_param, clearmem_params_creator(args, dirs); end
    end
    
    if motion_confound
        clearmem_motion_confound(args, dirs);
    end
end

%% ============= LOCALIZER CLASSIFICATION 
if sum(localizer_classification)
    xphase                  = 1;%1_localizer, 2_study
    args.xphase             = xphase;
    args.rest               = xrest;%'rest';% including rest or not
    
    args.train_phase        = 'localizer';
    args.train_regress_type = 'shift';%'beta'
    args.test_phase         = args.train_phase;
    args.test_regress_type  = args.train_regress_type;
    
    args.cross_valid        = 1;%1 is default
    args.cross_decoding     = 0;%0 is default
    
    %*************** base name
    if strcmp(args.rest, 'norest')
        args.analysis_basename = sprintf('%s_%s_%s_%s_shift%str_%s', ...
            args.phase_name{xphase}, args.mask_name, args.epi_name, ...
            args.level, num2str(args.shift_TRs), args.rest);
    elseif strcmp(args.rest, 'rest')
        args.analysis_basename = sprintf('%s_%s_%s_%s_shift%str', ...
            args.phase_name{xphase}, args.mask_name, args.epi_name, ...
            args.level, num2str(args.shift_TRs));
    end

    if args.class_selecting
        args.analysis_basename = sprintf('cate%s_%s', sprintf('%d', args.selected_category), args.analysis_basename);
        
        if args.subclass_selecting 
            args.analysis_basename = sprintf('subcate%s_%s', sprintf('%d', args.selected_subcategory), args.analysis_basename);
        end
    end
    
    %*************** fixed penalty    
    if args.fixed_penalty{xphase}
        args.analysis_basename = sprintf('%s_fixpen%s', args.analysis_basename, num2str(args.xpenalty));
    end
    
    if sum(localizer_classification(1:3))
        for xsub = args.filtered_subs
            %*************** setup subject & directories
            args.subject_id     = subject_lists(xsub).name;
            dirs                = setup_directory(dirs, args);
            
            %*************** make directories
            if ~isdir(dirs.mvpa.home{xphase}), mkdir(dirs.mvpa.home{xphase}); end
            if ~isdir(dirs.mvpa.output{xphase}), mkdir(dirs.mvpa.output{xphase}); end
            if ~isdir(dirs.mvpa.scratch{xphase}), mkdir(dirs.mvpa.scratch{xphase}); end
            if ~isdir(dirs.mvpa.selected_mask{xphase}), mkdir(dirs.mvpa.selected_mask{xphase}); end
            if ~isdir(dirs.mvpa.regressor{xphase}), mkdir(dirs.mvpa.regressor{xphase}); end
            if ~isdir(dirs.mvpa.parse{xphase}), mkdir(dirs.mvpa.parse{xphase}); end
            if ~isdir(dirs.param), mkdir(dirs.param); end
            
            %*************** design parameters
            clear xindex
            fname = fullfile(dirs.param, sprintf('localizer_parameters_%s.mat', args.subject_id));%,'index'
            if exist(fname, 'file'), xindex{1} = load(fname); args.index{1} = xindex{1}.index;
            else args.index{1} = create_design_index_localizer(args, dirs); end
            
            %*************** regressors
            args.regs{1} = create_mvpa_regressor_localizer(args, dirs);
            
            %*************** running classification
            if localizer_classification(1), clearmem_localizer_classification_01(args, dirs); end
            if localizer_classification(2), clearmem_localizer_classification_02(args, dirs); end
            if localizer_classification(3), clearmem_localizer_classification_03(args, dirs); end
        end
    end
    
    %*************** 1st + 2nd level importance map analyses
    if sum(localizer_classification(4:5))
        %*************** setup subject & directories
        args.subject_id   = subject_lists(3).name;
        dirs              = setup_directory(dirs, args);
        
        %*************** design parameters
        clear xindex
        if strcmp(xcluster,'local')
            fname = fullfile(dirs.param, sprintf('localizer_parameters_%s.mat', args.subject_id));%,'index'
            if exist(fname, 'file'), xindex{1} = load(fname); args.index{1} = xindex{1}.index;
            else args.index{1} = create_design_index_localizer(args, dirs); end
        else
            args.index{1} = create_design_index_localizer(args, dirs);
        end
        
        %*************** regressors
        clear xregs
        args.regs{1}  = create_mvpa_regressor_localizer(args, dirs);
        
        %********************************
        % ROC (Receiver operating characteristic)
        %********************************
        if localizer_classification(4)
            %*************** create mvpa outcome directory
            if ~isdir(dirs.mvpa.group.home{xphase}), mkdir(dirs.mvpa.group.home{xphase}); end
            if ~isdir(dirs.mvpa.group.out{xphase}), mkdir(dirs.mvpa.group.out{xphase}); end
            if ~isdir(dirs.mvpa.group.mvpa_result{xphase}), mkdir(dirs.mvpa.group.mvpa_result{xphase}); end
            if ~isdir(dirs.mvpa.group.auc{xphase}), mkdir(dirs.mvpa.group.auc{xphase}); end
%             if ~isdir(dirs.mvpa.group.auc_baseline{xphase}), mkdir(dirs.mvpa.group.auc_baseline{xphase}); end
            
            %*************** group analysis
            clearmem_localizer_classification_04(args, dirs);
        end
        
        %********************************
        % IMPORTANCE MAP
        %********************************
        if localizer_classification(5)
            %*************** create mvpa outcome directory
            if ~isdir(dirs.mvpa.group.imp_map{xphase}), mkdir(dirs.mvpa.group.imp_map{xphase}); end
            
            %*************** importance map analysis
            args.refer_mask   = fullfile(dirs.ref_mask, sprintf('%s.%s', args.mask_name, args.epiext));
            args.subject_list = subject_lists;
            clearmem_localizer_importance_map(args, dirs);
        end
    end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%============= DECODING WORKING MEMORY (TRAIN:LOCALIZER + TEST:STUDY)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if sum(mvpa_decode)
    xphase                  = 2;%1_localizer, 2_study, 3_operation
    args.xphase             = xphase;
    args.rest               = xrest;
    
    args.train_phase        = 'localizer';
    args.train_regress_type = 'shift';
    args.test_phase         = 'study';
    args.test_regress_type  = 'shift';
        
    args.cross_valid        = 0;%0 is default
    args.cross_decoding     = 0;%0 is default
    
%     if strcmp(args.level, 'category')
%         args.class_selecting = 0; it_sels = 1;
%     elseif strcmp(args.level, 'subcategory')
%         args.class_selecting = 1; it_sels = 1:3; 
%     end

    %*************** base name
    if strcmp(args.rest, 'norest')
        args.analysis_basename = sprintf('%s_%s_%s_%s_shift%str_%s', ...
            args.phase_name{xphase}, args.mask_name, args.epi_name, ...
            args.level, num2str(args.shift_TRs), args.rest);
    elseif strcmp(args.rest, 'rest')
        args.analysis_basename = sprintf('%s_%s_%s_%s_shift%str', ...
            args.phase_name{xphase}, args.mask_name, args.epi_name, ...
            args.level, num2str(args.shift_TRs));
    end
    
    if strcmp(args.level, 'subcategory') && (args.class_selecting)
        args.analysis_basename = sprintf('cate%s_%s', sprintf('%d', args.selected_category), args.analysis_basename);
    end
    
    %*************** ANALYSIS PARAMETERS
    args.cross_valid       = 0;% 1_if n-1 cross-validation (within localizer only) <--------------------!!!
    
    %*************** MVPA_CLASSIFICATION
    if sum(mvpa_decode(1:3))
        for xsub = args.filtered_subs
            %*************** setup subject & directories
            args.subject_num    = xsub;
            args.subject_id     = subject_lists(xsub).name;
            dirs                = setup_directory(dirs, args);
            args.mask_dir       = dirs.mvpa.selected_mask;
            
            fprintf('...sub %s: %s\n', num2str(xsub), args.subject_id);
            
            %*************** make directories
            if ~isdir(dirs.mvpa.home{xphase}), mkdir(dirs.mvpa.home{xphase}); end
            if ~isdir(dirs.mvpa.output{xphase}), mkdir(dirs.mvpa.output{xphase}); end
            if ~isdir(dirs.mvpa.scratch{xphase}), mkdir(dirs.mvpa.scratch{xphase}); end
            if ~isdir(dirs.mvpa.selected_mask{xphase}), mkdir(dirs.mvpa.selected_mask{xphase}); end
            if ~isdir(dirs.mvpa.regressor{xphase}), mkdir(dirs.mvpa.regressor{xphase}); end
            if ~isdir(dirs.mvpa.parse{xphase}), mkdir(dirs.mvpa.parse{xphase}); end
            if ~isdir(dirs.param), mkdir(dirs.param); end
            
            %*************** design parameters
            clear xindex
            fname = fullfile(dirs.param, sprintf('localizer_parameters_%s.mat', args.subject_id));%,'index'
            if exist(fname, 'file'), xindex{1} = load(fname); args.index{1} = xindex{1}.index;
            else args.index{1} = create_design_index_localizer(args, dirs); end
            
            fname = fullfile(dirs.param, sprintf('study_parameters_%s.mat', args.subject_id));%,'index'
            if exist(fname, 'file'), xindex{2} = load(fname); args.index{2} = xindex{2}.index;
            else args.index{2} = create_design_index_study(args, dirs); end
            
            %*************** regressors: localizer
            args.regs{1} = create_mvpa_regressor_localizer(args, dirs);
            
            %*************** regressors: study
            args.regs{2} = create_mvpa_regressor_study(args, dirs);
            
            %*************** classification
            if mvpa_decode(1), clearmem_mvpa_decode_01(args, dirs); end
            if mvpa_decode(2), clearmem_mvpa_decode_02(args, dirs); end
            if mvpa_decode(3), clearmem_mvpa_decode_03(args, dirs); end
        end
    end
    
    %*************** 2nd level group analyses
    if sum(mvpa_decode(4:5))
        %*************** setup subject & directories
        args.subject_id   = subject_lists(3).name;
        dirs              = setup_directory(dirs, args);
        
        %*************** design parameters
        clear xindex
        xindex{1} = load(fullfile(dirs.param, sprintf('localizer_parameters_%s.mat', args.subject_id)));
        xindex{2} = load(fullfile(dirs.param, sprintf('study_parameters_%s.mat', args.subject_id)));
        for i = 1:2, args.index{i} = xindex{i}.index; end
        
        %*************** regressors
        clear xregs
        %*************** regressors: localizer
        args.regs{1} = create_mvpa_regressor_localizer(args, dirs);
        
        %*************** regressors: study
        args.regs{2} = create_mvpa_regressor_study(args, dirs);
        
        %********************************
        % 2nd level/ TIMECOURSE
        %********************************
        if mvpa_decode(4)
            %*************** create mvpa outcome directory
            if ~isdir(dirs.mvpa.group.home{xphase}), mkdir(dirs.mvpa.group.home{xphase}); end
            if ~isdir(dirs.mvpa.group.out{xphase}), mkdir(dirs.mvpa.group.out{xphase}); end
            if ~isdir(dirs.mvpa.group.mvpa_result{xphase}), mkdir(dirs.mvpa.group.mvpa_result{xphase}); end
            
            %*************** group analysis
            clearmem_mvpa_decode_04(args, dirs);
        end
        
        %********************************
        % AUC/ ACCURACY (Receiver operating characteristic)
        %********************************
        if mvpa_decode(5)
            %*************** create mvpa outcome directory
            if ~isdir(dirs.mvpa.group.home{xphase}), mkdir(dirs.mvpa.group.home{xphase}); end
            if ~isdir(dirs.mvpa.group.out{xphase}), mkdir(dirs.mvpa.group.out{xphase}); end
            if ~isdir(dirs.mvpa.group.mvpa_result{xphase}), mkdir(dirs.mvpa.group.mvpa_result{xphase}); end
            if ~isdir(dirs.mvpa.group.auc{xphase}), mkdir(dirs.mvpa.group.auc{xphase}); end
            
            %*************** group analysis
            args.subject_list = subject_lists;
            clearmem_mvpa_decode_05(args, dirs);
        end
    end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%============= OPERATION DECODING (TRAIN:STUDY + TEST:STUDY)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if sum(operation_mvpa_decode)
    xphase                  = 3;%1_localizer, 2_study, 3_operation
    args.xphase             = xphase;
    args.rest               = 'norest';%xrest;
    
    %================= study arguments
    args.train_phase        = 'study';
    args.train_regress_type = 'shift';
    args.test_phase         = 'study';
    args.test_regress_type  = 'shift';
    
    args.cross_valid        = 1;%1 is default: n-1 cross-validation 
    args.cross_decoding     = 1;%1 cross validation + decoding on all testing timepoints (non-shifted)
    
    %*************** base name
    args.analysis_basename = sprintf('%s_%s_%s_shift%str_%s', ...
        args.phase_name{xphase}, args.mask_name, args.epi_name, ...
        num2str(args.shift_TRs), args.rest);
    
    %*************** fixed penalty
    if args.fixed_penalty{xphase}
        args.analysis_basename = sprintf('%s_fixpen%s', args.analysis_basename, num2str(args.xpenalty));
    end
    
    %*************** ANALYSIS PARAMETERS
    args.regress_type       = args.train_regress_type;
    
    %*************** MVPA_CLASSIFICATION
    if sum(operation_mvpa_decode(1:3))
        for xsub = args.filtered_subs%args.g_sub
            %*************** setup subject & directories
            args.subject_num = xsub;
            args.subject_id  = subject_lists(xsub).name;
            dirs             = setup_directory(dirs, args);
            
            %*************** make directories
            if ~isdir(dirs.mvpa.home{xphase}), mkdir(dirs.mvpa.home{xphase}); end
            if ~isdir(dirs.mvpa.output{xphase}), mkdir(dirs.mvpa.output{xphase}); end
            if ~isdir(dirs.mvpa.scratch{xphase}), mkdir(dirs.mvpa.scratch{xphase}); end
            if ~isdir(dirs.mvpa.selected_mask{xphase}), mkdir(dirs.mvpa.selected_mask{xphase}); end
            if ~isdir(dirs.mvpa.regressor{xphase}), mkdir(dirs.mvpa.regressor{xphase}); end
            if ~isdir(dirs.mvpa.parse{xphase}), mkdir(dirs.mvpa.parse{xphase}); end
            if ~isdir(dirs.param), mkdir(dirs.param); end
            if ~isdir(dirs.mvpa.group.auc_baseline{xphase}), mkdir(dirs.mvpa.group.auc_baseline{xphase}); end
            
            %*************** design parameters
            clear xindex
            args.index{xphase} = create_design_index_study(args, dirs);
            
            %*************** regressors
            if args.four_oper_regress
                args.regs{xphase} = create_mvpa_regressor_operation_four(args, dirs);
            else
                args.regs{xphase} = create_mvpa_regressor_operation(args, dirs);
            end
            
            %*************** classification
            if operation_mvpa_decode(1),     clearmem_mvpa_operation_01(args, dirs); end
            if operation_mvpa_decode(2), clearmem_mvpa_operation_02(args, dirs); end
            if operation_mvpa_decode(3), clearmem_mvpa_operation_03(args, dirs); end
        end
    end
    
    %********************************
    % 2nd level group analyses
    %********************************
    if sum(operation_mvpa_decode(4:7))
        
        %*************** setup subject & directories
        args.subject_id = xrefer_subj;
        dirs            = setup_directory(dirs, args);
        
        %*************** design parameters
        clear xindex
        xindex = load(fullfile(dirs.param, sprintf('study_parameters_%s.mat', args.subject_id)));
        args.index{xphase} = xindex.index;
        
        %*************** regressors
        %*************** regressors
        if args.four_oper_regress
            args.regs{xphase} = create_mvpa_regressor_operation_four(args, dirs);
        else
            args.regs{xphase} = create_mvpa_regressor_operation(args, dirs);
        end
        
        %********************************
        % 2nd level group analyses
        %********************************
        if operation_mvpa_decode(4)
            %*************** create mvpa outcome directory
            if ~isdir(dirs.mvpa.group.home{xphase}), mkdir(dirs.mvpa.group.home{xphase}); end
            if ~isdir(dirs.mvpa.group.out{xphase}), mkdir(dirs.mvpa.group.out{xphase}); end
            if ~isdir(dirs.mvpa.group.mvpa_result{xphase}), mkdir(dirs.mvpa.group.mvpa_result{xphase}); end
            
            %*************** group analysis
            clearmem_mvpa_operation_04(args, dirs);
        end
        
        %********************************
        % AUC/ ACCURACY (Receiver operating characteristic)
        %********************************
        if operation_mvpa_decode(5)
            %*************** create mvpa outcome directory
            if ~isdir(dirs.mvpa.group.home{xphase}), mkdir(dirs.mvpa.group.home{xphase}); end
            if ~isdir(dirs.mvpa.group.out{xphase}), mkdir(dirs.mvpa.group.out{xphase}); end
            if ~isdir(dirs.mvpa.group.mvpa_result{xphase}), mkdir(dirs.mvpa.group.mvpa_result{xphase}); end
            if ~isdir(dirs.mvpa.group.auc{xphase}), mkdir(dirs.mvpa.group.auc{xphase}); end
            if ~isdir(dirs.mvpa.group.auc_baseline{xphase}), mkdir(dirs.mvpa.group.auc_baseline{xphase}); end
            
            %*************** group analysis
            clearmem_mvpa_operation_05(args, dirs);
        end
    
        %********************************
        % IMPORTANCE MAP
        %********************************
        if operation_mvpa_decode(6)
            args.refer_mask = fullfile(dirs.ref_mask, sprintf('%s.%s', args.mask_name, args.epiext));
            
            %*************** first level
            if args.subj_impmap{xphase}
                for xsub = args.filtered_subs
                    %*************** setup subject & directories
                    args.subject_num = xsub;
                    args.subject_id  = subject_lists(xsub).name;
                    dirs             = setup_directory(dirs, args);
                    
                    if ~isdir(dirs.mvpa.imp_map{xphase}), mkdir(dirs.mvpa.imp_map{xphase}); end
                    
                    %*************** 1st-level importance map analysis
                    clearmem_operation_importance_map_subj(args, dirs);
                end
            end
            
            %*************** create mvpa outcome directory
            if ~isdir(dirs.mvpa.group.imp_map{xphase}), mkdir(dirs.mvpa.group.imp_map{xphase}); end
            
            %*************** importance map analysis
            clearmem_operation_importance_map(args, dirs);
        end

        %********************************
        % CROSS-SUBJECT CLASSIFICATION
        %********************************
        if operation_mvpa_decode(7)
            % normalize to MNI: args.mni_brain > smoothing 6mm
            args.epi_name = 'bold_mcf_brain_hpass_dt_mni';
            
            %*************** first level
            if args.subj_norm
                for xsub = args.filtered_subs
                    %*************** setup subject & directories
                    args.subject_num = xsub;
                    args.subject_id  = subject_lists(xsub).name;
                    dirs             = setup_directory(dirs, args);
                    
                    clearmem_operation_cross_subject(args, dirs);
                    
                end
            end
        end
    end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%============= RSA WORKING MEMORY (ROI:LOCALIZER + RSA:STUDY)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if sum(rsa_decode) || (rsa_spm_result)
    xphase                  = 2;%1_localizer, 2_study, 3_operation
    args.xphase             = xphase;
    args.rest               = xrest;
    
    args.train_phase        = 'localizer';
    args.train_regress_type = 'boxcar';
    args.test_phase         = 'study';
    
    %*************** base name
    args.analysis_basename = sprintf('%s_%s_%s_%s_shift%str_%s', ...
        args.phase_name{xphase}, args.rsa_mask, args.epi_name, ...
        args.level, num2str(args.shift_TRs), args.rest);
    
    %*************** TEMPATE PATTERN EXTRACTION
    % step 1. ROI extraction: localizer 
    % step 2. mean the extracted pattern
    if rsa_decode(1)
        if strcmp(xlevel, 'category')
            args.epi_name = args.spm_epi_name;
        end
        for xsub = args.filtered_subs
            %*************** setup subject & directories
            args.subject_num    = xsub;
            args.subject_id     = subject_lists(xsub).name;
            dirs                = setup_directory(dirs, args);
            
            %*************** make directories
            if ~isdir(dirs.rsa.home); mkdir(dirs.rsa.home); end
            if ~isdir(dirs.rsa.roi.home); mkdir(dirs.rsa.roi.home); end
            if ~isdir(dirs.rsa.roi.regressor); mkdir(dirs.rsa.roi.regressor); end
            if ~isdir(dirs.rsa.roi.spm); mkdir(dirs.rsa.roi.spm); end
            
            %*************** design parameters
            clear xindex
            if strcmp(xcluster, 'local')
                args.index{1} = create_design_index_localizer(args, dirs, 0);
            else
                args.index{1} = create_design_index_localizer(args, dirs, 1);
            end
            
            %*************** regressors: localizer
            args.regs{1}  = create_spm_regressor_localizer(args, dirs);
            
            %*************** classification
            clearmem_rsa_decode_01(args, dirs);
            
        end
    elseif rsa_decode(2)
        for xsub = args.filtered_subs
            %*************** setup subject & directories
            args.subject_num    = xsub;
            args.subject_id     = subject_lists(xsub).name;
            dirs                = setup_directory(dirs, args);
            
            %*************** make directories
            if ~isdir(dirs.rsa.pattern); mkdir(dirs.rsa.pattern); end
            if ~isdir(dirs.rsa.scratch); mkdir(dirs.rsa.scratch); end
            if ~isdir(dirs.rsa.group.home); mkdir(dirs.rsa.group.home); end
            if ~isdir(dirs.rsa.group.pattern); mkdir(dirs.rsa.group.pattern); end
            if ~isdir(dirs.rsa.group.parse); mkdir(dirs.rsa.group.parse); end
        
            %*************** design parameters
            clear xindex
            xindex{1} = load(fullfile(dirs.param, sprintf('localizer_parameters_%s.mat', args.subject_id)));
            xindex{2} = load(fullfile(dirs.param, sprintf('study_parameters_%s.mat', args.subject_id)));
            for i = 1:2, args.index{i} = xindex{i}.index; end
            
            %*************** regressors
            clear xregs
            %*************** regressors: localizer
            args.regs{1} = create_mvpa_regressor_localizer(args, dirs);
            
            %*************** regressors: study
            args.regs{2} = create_mvpa_regressor_study(args, dirs);
            
            %*************** rsa
            if ~strcmp(args.level, 'item')
                clearmem_rsa_decode_02(args, dirs);
            else
                clearmem_rsa_decode_02_item(args, dirs);
            end
            
        end
    end
    
    %********************************
    % 2nd level group analyses
    %********************************
    if sum(rsa_decode(3:4))
        %*************** setup subject & directories
        args.subject_id   = subject_lists(3).name;
        dirs              = setup_directory(dirs, args);
        
        %*************** design parameters
        clear xindex
        xindex{1} = load(fullfile(dirs.param, sprintf('localizer_parameters_%s.mat', args.subject_id)));
        xindex{2} = load(fullfile(dirs.param, sprintf('study_parameters_%s.mat', args.subject_id)));
        for i = 1:2, args.index{i} = xindex{i}.index; end
        
        %*************** regressors
        clear xregs
        %*************** regressors: localizer
        args.regs{1} = create_mvpa_regressor_localizer(args, dirs);
        
        %*************** regressors: study
        args.regs{2} = create_mvpa_regressor_study(args, dirs);
        
        %*************** create mvpa outcome directory
        if ~isdir(dirs.rsa.group.home); mkdir(dirs.rsa.group.home); end
        if ~isdir(dirs.rsa.group.pattern); mkdir(dirs.rsa.group.pattern); end
        if ~isdir(dirs.rsa.group.parse); mkdir(dirs.rsa.group.parse); end
        
        %*************** group analysis
        if rsa_decode(3)
            if ~strcmp(args.level, 'item')
                clearmem_rsa_decode_03(args, dirs);
            else
                clearmem_rsa_decode_03_item(args, dirs);
            end
            
        elseif rsa_decode(4) && (strcmp(args.level, 'item'))
            clearmem_rsa_decode_04_bootstrap(args, dirs)
        end
    end

    %********************************
    % 1st level spm result after 1st step
    %********************************
    if rsa_spm_result
        if strcmp(xlevel, 'category') 
            args.epi_name = args.spm_epi_name;
        end
        for xsub = args.filtered_subs
            %*************** setup subject & directories
            args.subject_num    = xsub;
            args.subject_id     = subject_lists(xsub).name;
            dirs                = setup_directory(dirs, args);
            
            %*************** make directories
            if ~isdir(dirs.rsa.home); mkdir(dirs.rsa.home); end
            if ~isdir(dirs.rsa.roi.home); mkdir(dirs.rsa.roi.home); end
            if ~isdir(dirs.rsa.roi.regressor); mkdir(dirs.rsa.roi.regressor); end
            if ~isdir(dirs.rsa.roi.spm); mkdir(dirs.rsa.roi.spm); end
            
            %*************** classification
            clearmem_rsa_decode_01_tmp(args, dirs);
        end
    end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%============= WORKING MEMORY OPERATION DECODING (TRAIN:STUDY + TEST:STUDY)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if sum(operation_mvpa_stw)
    xphase                  = 3;%1_localizer, 2_study, 3_operation
    args.xphase             = xphase;
    args.rest               = 'norest';%xrest;
    
    %================= study arguments
    args.train_phase        = 'study';
    args.train_regress_type = 'shift';
    args.test_phase         = 'study';
    args.test_regress_type  = 'shift';
    
    args.cross_valid        = 1;%1 is default: n-1 cross-validation 
    args.cross_decoding     = 0;%1 cross validation + decoding on all testing timepoints (non-shifted)
    
    %*************** base name
    args.analysis_basename = sprintf('stw_%s_%s_%s_shift%str_%s', ...
        args.phase_name{xphase}, args.mask_name, args.epi_name, ...
        num2str(args.shift_TRs), args.rest);
    
    %*************** ANALYSIS PARAMETERS
    args.regress_type       = args.train_regress_type;
    
    %*************** MVPA_CLASSIFICATION
    if sum(operation_mvpa_stw(1:2))
        for xsub = args.filtered_subs
            %*************** setup subject & directories
            args.subject_num = xsub;
            args.subject_id  = subject_lists(xsub).name;
            dirs             = setup_directory(dirs, args);
            
            %*************** make directories
            if ~isdir(dirs.mvpa.home{xphase}), mkdir(dirs.mvpa.home{xphase}); end
            if ~isdir(dirs.mvpa.output{xphase}), mkdir(dirs.mvpa.output{xphase}); end
            if ~isdir(dirs.mvpa.scratch{xphase}), mkdir(dirs.mvpa.scratch{xphase}); end
            if ~isdir(dirs.mvpa.selected_mask{xphase}), mkdir(dirs.mvpa.selected_mask{xphase}); end
            if ~isdir(dirs.mvpa.parse{xphase}), mkdir(dirs.mvpa.parse{xphase}); end
            
            %*************** design parameters
            clear xindex
            xindex = load(fullfile(dirs.param, sprintf('study_parameters_%s.mat', args.subject_id)));
            args.index{xphase} = xindex.index;
            
            %*************** classification
            if operation_mvpa_stw(1), clearmem_stw_operation_01(args, dirs); end
            if operation_mvpa_stw(2), clearmem_stw_operation_02(args, dirs); end
        end
    end
    
    %********************************
    % 2nd level group analyses
    %********************************
    if sum(operation_mvpa_stw(3))
        
        %*************** setup subject & directories
        args.subject_id = subject_lists(3).name;
        dirs            = setup_directory(dirs, args);
        
        %*************** design parameters
        clear xindex
        xindex = load(fullfile(dirs.param, sprintf('study_parameters_%s.mat', args.subject_id)));
        args.index{xphase} = xindex.index;
       
        %********************************
        % 2nd level group analyses
        %********************************
        %*************** create mvpa outcome directory
        if ~isdir(dirs.mvpa.group.home{xphase}), mkdir(dirs.mvpa.group.home{xphase}); end
        if ~isdir(dirs.mvpa.group.out{xphase}), mkdir(dirs.mvpa.group.out{xphase}); end
        if ~isdir(dirs.mvpa.group.mvpa_result{xphase}), mkdir(dirs.mvpa.group.mvpa_result{xphase}); end
        
        %*************** group analysis
        clearmem_stw_operation_03(args, dirs);
    end
end