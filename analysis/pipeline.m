function[] = clearmem_pipeline(xstep, xlevel, xmask)

%**********************************************************
    %%%%% xstep:
    %%%%%   'mvpa_localizer': localizer mvpa
    %%%%%   'mvpa_WM': working memory mvpa
    %%%%%   'mvpa_operation': operation
    %%%%%   'rsa_WM': working memory rsa
    %%%%%   'mvpa_operation_stw': sliding-time-window operation
    %%%%% xlevel: 'category' / 'subcategory'
    %%%%% xmask: 'whole' / 'vvs'
%**********************************************************

localizer_classification = 0; 
mvpa_decode = 0; 
operation_mvpa_decode = 0; 
rsa_decode = 0; 
operation_mvpa_stw = 0;

args.four_oper_regress = 1;% 1_using 4 operations (remove all repSubcate)

%% ============= SETUP TOOLBOX & DIRECTORIES
%*************** MVPA PRINCETON TOOLBOX
addpath ~/github/mvpa/
mvpa_add_paths;
addpath(genpath('~/github/GroupICATv4.0b/icatb/'))

%*************** SPM TOOLBOX
addpath(genpath('/Applications/spm12'))
spm('Defaults','fMRI'); spm_jobman('initcfg');

%*************** HOMEMADE CODES
addpath ./homemade/
addpath ./spm_custom/
args.analysis = pwd;

%% ============= SETUP DIRECTORIES & ARGS
curr_dir       = pwd; cd ~
dirs.home      = pwd; cd(curr_dir);
dirs.project   = fullfile(dirs.home,'Clearmem');
dirs.fmri_data = fullfile(dirs.project,'imaging_data');

%% ============= SETUP SUBJECTS
subject_lists      = dir(fullfile(dirs.fmri_data,'sub*'));
args.subject_list  = subject_lists;
args.n_sub         = length(subject_lists);
args.g_sub         = 1:args.n_sub;

excluded_subs      = unique([2, 7, 8, 16, 18]); % subjects who slept
auc_cutoff_subs    = [4, 10, 14, 24, 40]; % operation classifier ACU < .5

args.g_sub         = args.g_sub(~ismember(args.g_sub, excluded_subs));
args.filtered_subs = args.g_sub(~ismember(args.g_sub, auc_cutoff_subs));
xrefer_subj        = 'sub004';

%% ============= CLASSIFIER ARGUMENTS
switch(xmask)
    case'whole'
        args.mask_name  = 'highres001_brain_grey_3mm_bin0.2';
    case 'vvs'
        args.mask_name  = 'mask_vvps_temporal_fined';
end

args.rest          = 'norest';
args.trim_trs      = 10;% trimming first # of TRs 
args.level         = xlevel;%'category';%'subcategory';%

%*************** REGRESOR PARAMETER
args.shift_TRs     = 10;%str2double(xshift);% 10, 14

%*************** COLOR SETUP
args.cond_color{1} = [0, 149, 137]/255;
args.cond_color{2} = [0, 134, 196]/255;
args.cond_color{3} = [51, 153, 153]/255;
args.cond_color{4} = [189, 0, 6]/255;
args.cond_color{5} = [253, 81, 0]/255;

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
args.regress_type  = 'shift';
args.featSelThresh = '0.05';%'0.05';%'0.05';%
args.penalty       = [0 0.01 0.1 1.0 10 100 1000 10000];%[0 1 10:10:70 100:100:1000 10000];
args.n_pen         = 10;%step2 narrow penalty numbers
args.classifier    = 'L2logreg';%'L2logreg' & 'ridge' 'logreg'
args.repetitions   = 1;%n_iterations

%*************** IMPORTANCE MAP
args.imp_type      = 'mcduff';
args.peak_thresh   = 3;%10;%2;%0.1;% top 2sd
dirs.fsl = '/usr/local/fsl';

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
args.tc_tr_disp    = 30;%for stats/display
args.n_tr_blks     = 3;%6
args.alpha         = 0.05;

args.item_peak     = [11 16; 22 27];%for replace: 1_before, 2_after switch
args.baseline_trs  = 6;
args.ttest_trs     = 7:21;%24
args.anova_trs     = 7:21;%1:args.tc_tr_disp;

%*************** operation
args.stw{3}        = 0;
args.xstep         = xstep;

%% ============= SETUP ANALYSES PROCEDURE
% NOTE: args.load_pattern{phase}: EPI patterns are loaded and saved in the scratch
% so, if the pattern has been loaded before put 0 to just load the
% patterns from subject structure

switch(xstep)   
    case 'mvpa_localizer'
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LOCALIZER
        % step 1: read EPI patterns + regressors/selectors + feature selection
        % step 2: classification (with all penalties if args.penalty_check = 1)
        %         if args.penalty_check = 0, load maximum penalty
        % step 3: parse the mvpa outcome
        % step 4: AUC 2nd
        localizer_classification   = [0 0 0 0];
        
        % step 2
        args.fixed_penalty{1}      = 0;
        if args.fixed_penalty{1}, args.xpenalty = 50; end
        args.penalty_check{1}      = 1;% 0_load max penalty
        
    case 'mvpa_WM'
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% WORKING MEMORY DECODING
        % step 1: train: EPI patterns + regressors + feature selection
        %         test: EPI patterns under the feature-selected mask with training data
        %               + regressors
        %         decoding setup: combine train & test: EPI + regressors/selectors
        % step 2: load penalty from localizer classifier + classification
        % step 3: parse the mvpa outcome + 1st level analysis
        % step 4: 2nd level analysis + timecourse
        mvpa_decode                    = [0 0 0 0];
        
    case 'mvpa_operation'
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% WORKING MEMORY OPERATION DECODING
        % step 1: read EPI patterns + regressors/selectors + feature selection
        %         decoding setup: decode all timepoints for the testing run
        % step 2: classification (penalty: 50)
        % step 3: parse the mvpa outcome
        % step 4: 2nd level analysis
        % step 5: AUC/ROC 2nd
        % step 6: importance map 1st + 2nd
        
        operation_mvpa_decode          = [0 0 0 0 0 0];
        
        % step 2
        args.fixed_penalty{3}          = 0;
        if args.fixed_penalty{3}, args.xpenalty = 50; end
        args.penalty_check{3}          = 1;% 0_load max penalty
        
        if sum(operation_mvpa_decode(4:7))            
            % step 6: importance map
            args.subj_impmap{3}        = 1;% 1_run 1st-level importance map
        end
    case 'rsa_WM'
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% WORKING MEMORY DECODING
        % step 1. ROI extraction (GLM): localizer / 128 hf applied
        % step 2. create template by averaging the extracted pattern
        % step 3: RSA timecourse: test on study: 1st + 2nd level
        % step 4: item RSA bootstrap: N = N+1
        
        %*************** RESULT PARAMETERS
        rsa_decode               = [0 0 0 0];   
        
        % step 1
        args.spm_p_thresh        = 0.05;% t threshold of uncorrected p < 0.05
        args.spm_v_extent        = 10;% extent voxel
        args.spm_thresh_type     = 'none';%uncorrected 'FWE'; 'FDR';
        
        % step 2
        args.rsa_mask            = args.mask_name; 
        
        % step 3
        args.percept_win         = 11:16;%peak for presented item in study
        args.bin_trs{1}          = 13:18;%n-targ
        args.bin_trs{2}          = 11:16;%in synched (from 22trs)
        args.bin_trs{3}          = 5:10;%baseline in synched (from 22trs)
        
        % step 4
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
                
        if sum(operation_mvpa_stw(3))
            % step 3
            args.group_mvpa{3}   = 1;% 1_save 1st-level parsed mvpaout into group structure
            args.filtered_subs   = args.filtered_subs;%args.g_sub;
        end
end

%% ============= LOCALIZER CLASSIFICATION 
if sum(localizer_classification)
    xphase                  = 1;%1_localizer, 2_study
    args.xphase             = xphase;
    args.cross_valid        = 1;%1 is default
    args.cross_decoding     = 0;%0 is default
    
    %*************** base name
    args.analysis_basename = sprintf('%s_%s_%s_%s_shift%str_%s', ...
        args.phase_name{xphase}, args.mask_name, args.epi_name, ...
        args.level, num2str(args.shift_TRs), args.rest);
    
    %*************** fixed penalty    
    if args.fixed_penalty{xphase}
        args.analysis_basename = sprintf('%s_fixpen%s', args.analysis_basename, num2str(args.xpenalty));
    end
    
    if sum(localizer_classification(1:3))
        for xsub = args.filtered_subs
            %*************** setup subject & directories
            args.subject_id = subject_lists(xsub).name;
            dirs            = setup_directory(dirs, args);
            
            %*************** make directories
            if ~isdir(dirs.mvpa.home{xphase}), mkdir(dirs.mvpa.home{xphase}); end %#ok<*ISDIR>
            if ~isdir(dirs.mvpa.output{xphase}), mkdir(dirs.mvpa.output{xphase}); end
            if ~isdir(dirs.mvpa.scratch{xphase}), mkdir(dirs.mvpa.scratch{xphase}); end
            if ~isdir(dirs.mvpa.selected_mask{xphase}), mkdir(dirs.mvpa.selected_mask{xphase}); end
            if ~isdir(dirs.mvpa.regressor{xphase}), mkdir(dirs.mvpa.regressor{xphase}); end
            if ~isdir(dirs.mvpa.parse{xphase}), mkdir(dirs.mvpa.parse{xphase}); end
            if ~isdir(dirs.param), mkdir(dirs.param); end
            
            %*************** design parameters
            args.index{1} = create_design_index_localizer(args, dirs);
            
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
        args.subject_id = subject_lists(3).name;
        dirs            = setup_directory(dirs, args);
        
        %*************** design parameters
        args.index{1} = create_design_index_localizer(args, dirs);
        
        %*************** regressors
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
            
            %*************** group analysis
            clearmem_localizer_classification_04(args, dirs);
        end
    end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%============= DECODING WORKING MEMORY (TRAIN:LOCALIZER + TEST:STUDY)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if sum(mvpa_decode)
    xphase                  = 2;%1_localizer, 2_study, 3_operation
    args.xphase             = xphase;
    args.cross_valid        = 0;%0 is default
    args.cross_decoding     = 0;%0 is default

    %*************** base name
    args.analysis_basename = sprintf('%s_%s_%s_%s_shift%str_%s', ...
        args.phase_name{xphase}, args.mask_name, args.epi_name, ...
        args.level, num2str(args.shift_TRs), args.rest);
    
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
            args.index{1} = create_design_index_localizer(args, dirs);
            args.index{2} = create_design_index_study(args, dirs);
            
            %*************** regressors
            args.regs{1} = create_mvpa_regressor_localizer(args, dirs);
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
        args.index{1} = load(fullfile(dirs.param, sprintf('localizer_parameters_%s.mat', args.subject_id)));
        args.index{2} = load(fullfile(dirs.param, sprintf('study_parameters_%s.mat', args.subject_id)));
        
        %*************** regressors
        args.regs{1} = create_mvpa_regressor_localizer(args, dirs);
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
    
    %================= study arguments
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
            args.index{xphase} = create_design_index_study(args, dirs);
            
            %*************** regressors
            args.regs{xphase} = create_mvpa_regressor_operation(args, dirs);
            
            %*************** classification
            if operation_mvpa_decode(1), clearmem_mvpa_operation_01(args, dirs); end
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
        args.regs{xphase} = create_mvpa_regressor_operation(args, dirs);
        
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
    end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%============= RSA WORKING MEMORY (ROI:LOCALIZER + RSA:STUDY)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if sum(rsa_decode)
    xphase                  = 2;%1_localizer, 2_study, 3_operation
    args.xphase             = xphase;
    
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
            args.index{1} = create_design_index_localizer(args, dirs, 0);
            
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
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%============= WORKING MEMORY OPERATION DECODING (TRAIN:STUDY + TEST:STUDY)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if sum(operation_mvpa_stw)
    xphase                  = 3;%1_localizer, 2_study, 3_operation
    args.xphase             = xphase;
    args.rest               = 'norest';%xrest;
    
    %================= study arguments
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