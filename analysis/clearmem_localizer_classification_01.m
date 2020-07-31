function[] = clearmem_localizer_classification_01(args, dirs)

xph = args.xphase;

%% ============= LOAD REGRESSORS & SELECTORS
param = args.index{xph}.param;

%% ============= SETUP ARGS STRUCTURE
args.regress_type = args.train_regress_type;

%% ============= LOGGING
diary off;
fprintf('(+) using scratch dir: %s\n', dirs.mvpa.scratch{xph});
%*************** turn on diary to capture analysis output
diary(sprintf('%s/diary.txt', dirs.mvpa.scratch{xph}));
fprintf('running code: %s\n', mfilename)
fprintf('#####################################################################\n\n');
disp(args);
fprintf('#####################################################################\n');

%% ============= MASK CHECK
xmask = fullfile(dirs.mask, sprintf('%s.nii.gz', args.mask_name));

if ~exist(xmask, 'file')
    %*************** set FSL environment
    setenv('FSLDIR', dirs.fsl);  % this to tell where FSL folder is
    setenv('FSLOUTPUTTYPE', 'NIFTI_GZ'); % this to tell what the output type would be
    
    %*************** nii to roi.mat
    xmni_mask = fullfile(dirs.mni_mask, sprintf('%s.nii.gz', args.mask_name));
        
    fprintf('... converting mask %s \n', args.mask_name);
    
    %*************** make_mask.sh
    % # $1: FSL_DIR
    % # $2: subject's directory
    % # $3: MNI mask
    
    system(sprintf('./make_mask_mni2sub.sh %s %s %s', ...
        dirs.fsl, dirs.subj_home, xmni_mask))                
end

%% ============= Initializing subj. structure:start by creating an empty subj structure
% summarize(subj): summarize all info in subj structure
% get_object/set_object/set_objfield/set_objsubfield
% get_mat/set_mat

subj = init_subj(args.experiment, args.subject_id);%identifier of the subj
fprintf('\n(+) %s phase data\n\n', args.phase_name{xph});

%% ============= 01: EPI PATTERNS
%*************** load mask + read in epis
%*************** zsoring epis
%*************** option: read epis in mask | wholebrain
%*************** define: args.wholebrain '0 | 1'

fprintf('\n#####################################################################\n');
fprintf('* STEP_01: load mask + read/zscore epis patterns\n');

% if strcmp(args.level,'subcategory') && (args.class_selecting) 
%     if args.selected_category > 1
%         xload = 0;
%     else
%         xload = args.load_pattern{1};
%     end
% else
%     xload = args.load_pattern{1};
% end
xload = args.load_pattern{1};

nm.info = 'object_name';
if xload, [ph1.subj, ph1.nm] = mvpa_ph01_patterns_gz(args, subj, nm, dirs); end

%*************** ph1. base filename
ph1.basename = sprintf('%s_%s_zscored_%s', args.phase_name{xph}, args.mask_name, args.epi_name); 

%*************** save phase 1.
fname = sprintf('%s/ph1_%s.mat', dirs.mvpa.scratch{xph}, ph1.basename);
if xload, ph1.args = args; save(fname,'ph1','-v7.3');
else fprintf('  ... load patterns from ph1\n'); load(fname); end% load ph1

%*************** merge subj. structure
subj = ph1.subj; nm = ph1.nm;
summarize(subj)

%% ============= 01.0: RESET EPI PATTERNS
%*************** cf. not for missing runs but selecting runs from full set
if args.run_selecting
    %*************** reset patterns
    fprintf('\n  ... reset z-scored patterns for selected runs: %s \n', num2str(args.selected_run));
    fprintf('      from subj.patterns{2}.mat\n');
    
    selected_pattern = subj.patterns{2}.mat(:, ...
        1:(param.n_volumes/param.n_runs) * length(args.selected_run));
    
    subj = set_mat(subj, 'pattern', 'localizer_patterns_z', selected_pattern);
    summarize(subj, 'objtype', 'pattern');
    
    %*************** reset selectors
    fprintf('\n  ... reset selected for selected runs: %s \n', num2str(args.selected_run));
    
    selected_selector = args.regs{xph}.selectors;
    subj = set_mat(subj, 'selector', 'localizer_runs', selected_selector);
    summarize(subj, 'objtype', 'selector');
end

%% ============= 02: REGRESSORS
%*************** define regressors + selectors
%*************** regressors: shift or beta

fprintf('\n#####################################################################\n');
fprintf('* STEP_02: define regressors + selectors\n');

[ph2.subj, ph2.nm] = mvpa_ph02_regressors(args, subj, nm);

%*************** ph2. base filename
if strcmp(args.regress_type, 'shift')
    ph2.basename = sprintf('%s_%s_%s%dtr_blk_%s',...
        ph1.basename, args.level, args.regress_type, ...
        args.shift_TRs, args.rest);
elseif strcmp(args.regress_type, 'beta')
    ph2.basename = sprintf('%s_%s_%s', ...
        ph1.basename, args.level, args.regress_type);
end

%*************** reset ph2. filename
if args.class_selecting
    ph2.basename = sprintf('cate%s_%s', sprintf('%d', args.selected_category), ph2.basename);
end

if args.subclass_selecting
    ph2.basename = sprintf('subcate%s_%s', sprintf('%d',args.selected_subcategory), ph2.basename);
end

%*************** merge subj. structure
subj = ph2.subj; nm = ph2.nm;
summarize(subj)

%% ============= 03: FEATURE SELECTION
%*************** creating cross-validation indices
%*************** feature selection anova.
%*************** options: shiftTRs, peakWindow, featSelThresh

fprintf('\n#####################################################################\n');
fprintf('* STEP_03:feature selection\n');

[ph3.subj, ph3.args, ph3.nm] = mvpa_ph03_featselection(args, subj, nm, dirs);

%*************** ph3. base filename
if args.featVox
    ph3.basename = sprintf('%s_featsel_%svox', ph2.basename, num2str(args.fsVosNum));
else
    ph3.basename = sprintf('%s_featsel_thresh%s', ph2.basename, num2str(args.featSelThresh));
end

%*************** save phase 3.
fname = sprintf('%s/ph3_%s.mat', dirs.mvpa.scratch{xph}, ph3.basename);
save(fname,'ph3','-v7.3');

%*************** merge subj. structure
subj = ph3.subj; 
summarize(subj)

end
