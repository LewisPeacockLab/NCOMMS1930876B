function[] = clearmem_mvpa_decode_01(args, dirs)
%*************** patterns + regressor + feature selection
% train localizer + test study setup
% step 1: train: EPI patterns + regressors + feature selection
%         test: EPI patterns under the feature-selected mask with training data
%               + regressors
%         decoding setup: combine train & test: EPI + regressors/selectors

%% ============= UNPACK ARGS. 
xph               = args.xphase;
mask_name         = args.mask_name;

%% ============= LOGGING
diary off;
fprintf('(+) using scratch dir: %s\n', dirs.mvpa.scratch{xph});
%*************** turn on diary to capture analysis output
diary(sprintf('%s/diary.txt', dirs.mvpa.scratch{xph}));
fprintf('running code: %s\n', mfilename)
fprintf('#####################################################################\n\n');
disp(args);
fprintf('#####################################################################\n');

%% ============= Initializing subj. structure:start by creating an empty subj structure
% summarize(subj): summarize all info in subj structure
% get_object/set_object/set_objfield/set_objsubfield
% get_mat/set_mat

fprintf('\n(+) %s phase data\n\n', args.train_phase);

%% ============= 01:TRAINING: LOCALIZER: EPI PATTERNS
%*************** load mask + read in epis
%*************** zsoring epis
%*************** option: read epis in mask | wholebrain
%*************** define: args.wholebrain '0 | 1'
args.xphase       = 1;
xph               = args.xphase;

fprintf('\n#####################################################################\n');
fprintf('* STEP_01: TRAIN: load mask + read/zscore epis patterns\n');

%*************** ph1. base filename
ph1.basename = sprintf('%s_%s_zscored_%s', args.phase_name{xph}, args.mask_name, args.epi_name);
    
%*************** load phase 1.
load(sprintf('%s/ph1_%s.mat', dirs.mvpa.scratch{xph}, ph1.basename));

%*************** merge subj. structure
subj = ph1.subj; nm = ph1.nm;
summarize(subj)

%% ============= 02:TRAINING: LOCALIZER: REGRESSORS
%*************** define regressors + selectors
%*************** shift regressors: only shift for training trials

fprintf('\n#####################################################################\n');
fprintf('* STEP_02: TRAIN: define/shifting regressors + selectors\n');

[ph2.subj, ph2.nm] = mvpa_ph02_regressors(args, subj, nm);

%*************** ph2. base filename
ph2.basename = sprintf('%s_%s_%s%dtr_blk_%s',...
    ph1.basename, args.level, args.regress_type, ...
    args.shift_TRs, args.rest);

%*************** merge subj. structure
subj = ph2.subj; nm = ph2.nm;
summarize(subj)

%% ============= 03:TRAINING: LOCALIZER: FEATURE SELECTION
%*************** creating cross-validation indices
%*************** feature selection anova.
%*************** options: shiftTRs, peakWindow, featSelThresh

fprintf('\n#####################################################################\n');
fprintf('* STEP_03: TRAIN: feature selection\n');

[ph3.subj, ph3.args, ph3.nm] = mvpa_ph03_featselection(args, subj, nm, dirs);

%*************** ph3. base filename
ph3.basename = sprintf('%s_featsel_thresh%s', ...
    ph2.basename, num2str(args.featSelThresh));

%*************** save phase 3.
fname    = sprintf('%s/ph3_%s.mat', dirs.mvpa.scratch{2}, ph3.basename);
ph3.args = args;
save(fname,'ph3','-v7.3');

%*************** merge subj. structure
subj = ph3.subj; nm = ph3.nm;
summarize(subj)

%% ============= TESTING: STUDY: SETUP ARGS STRUCTURE
args.xphase       = 2;
xph               = args.xphase;
args.regress_type = args.test_regress_type;
args.featSel      = 1;%from localizer

disp(nm.mask_selec{1})

args.mask_name = sprintf('%s_%s_%s_sh%d_blk_%s', mask_name, ...
    nm.mask_selec{1}, args.level, args.shift_TRs, args.rest);

%*************** move selected mask to study folder
movefile(fullfile(dirs.mvpa.selected_mask{1}, ...
    sprintf('%s.nii.gz', nm.mask_selected_file{1}{1})), dirs.mvpa.selected_mask{2});

%% ============= 04: TESTING: STUDY: EPI PATTERNS
%*************** load mask + read in epis
%*************** zsoring epis
%*************** option: read epis in mask | wholebrain
%*************** define: args.wholebrain '0 | 1'

%%============= DEFINE EPIS & REGRESSORS & SELECTORS
for i=1:2%1_localizer, 2_study
    args.epi_names{i} = fullfile(dirs.runs{i}, ...
        sprintf('%s.%s', args.epi_name, args.epiext)); %#ok<*AGROW>
end

fprintf('\n#####################################################################\n');
fprintf('* STEP_04: TESTING: load mask + read/zscore epis patterns\n');

%*************** ph4. base file name
ph4.basename = sprintf('%s_sh%d_%s_fselected_%s_%s_%s_%s_zepi', ...
    args.phase_name{xph}, args.shift_TRs, args.rest, mask_name, ...
    args.featSelThresh, args.level, args.epi_name);
fname = sprintf('%s/ph4_%s.mat', dirs.mvpa.scratch{xph}, ph4.basename);

%*************** load patterns
[ph4.subj, ph4.nm] = mvpa_ph01_patterns_gz(args, subj, nm, dirs);
ph4.args = args; save(fname,'ph4','-v7.3');

%*************** merge subj. structure
subj = ph4.subj; nm = ph4.nm; 
summarize(subj)

%% ============= 05: TESTING: STUDY: REGRESSORS
%*************** define regressors + selectors

fprintf('\n#####################################################################\n');
fprintf('* STEP_05: TESTING: define/shifting regressors + selectors\n');

[ph5.subj, ph5.nm] = mvpa_ph02_regressors(args, subj, nm);

%*************** ph5. base filename
ph5.basename = sprintf('%s_%s', ph4.basename, args.regress_type);

%*************** merge subj. structure
subj = ph5.subj; nm = ph5.nm; 
summarize(subj)

%% ============= 06: DECODING SETUP

fprintf('\n#####################################################################\n');
fprintf('* STEP_06: setup decoding\n');

%*************** combine train & test epi patterns
% nm.pats_z{1}: train: zscored localizer eip patterns under mask
% nm.pats_z{2}: test: zscored study eip patterns under seledted_mask
% nm.mask_selected{1}: feature selected mask from train-runs
% for decoding, use feature-selected maks from train- for test-runs

nm.combined_pats    = nm.pats_z{2};

fprintf('\n(+) combine train & test epi patterns: %s\n', nm.combined_pats);

masked_train_pats_z = get_masked_pattern(subj, nm.pats_z{1}, subj.masks{2}.name);
masked_test_pats_z  = get_mat(subj,'pattern', nm.pats_z{2});
decode_epis         = [masked_train_pats_z masked_test_pats_z];

%--------------- remove train epis: don't need it along anymore
subj                = remove_mat(subj,'pattern', nm.pats_z{1});

%--------------- replace the test_epis structure with the combined epis.
% so that the pattern has information of what mask it was created by.
subj                = set_mat(subj,'pattern', nm.pats_z{2}, decode_epis);

%*************** combine train & test phase regressors
% nm.conds_sh{1}: train: shifted regressor
% nm.conds{2}: test: regressor

nm.combined_regs    = sprintf('%s+%s', nm.conds_sh{1}, nm.conds{2});

fprintf('\n(+) combine train & test phase regressors: %s\n', nm.combined_regs);

shifted_train_regs  = get_mat(subj,'regressors', nm.conds_sh{1});
test_regs           = get_mat(subj,'regressors', nm.conds{2});%non_shifted
decode_regs         = [shifted_train_regs test_regs];

subj                = initset_object(subj,'regressors', nm.combined_regs, decode_regs);

%*************** combine train & test phase selectors
% nm.train_selector: train selector
% nm.runs{2}: study selector
% set train selector: 1 & test selector: 2

nm.combined_runs = sprintf('%s+%s', nm.train_selector, nm.runs{2});

fprintf('\n(+) combine train & test phase selectors: %s\n', nm.combined_runs);

train_selector   = get_mat(subj,'selector', nm.train_selector);
study_selector   = get_mat(subj,'selector', nm.runs{xph});
test_selector    = 2 * (study_selector > 0);

decode_runs      = [train_selector test_selector];

subj             = initset_object(subj,'selector', nm.combined_runs, ...
                   decode_runs, 'group_name', nm.combined_runs);

%*************** save subj. structure
ph6.subj = subj; ph6.nm = nm; ph6.args = args;

%*************** base filename
ph6.basename     = sprintf('%s_decoding_setup_train+test', ph5.basename);

%*************** save phase 6.
fname = sprintf('%s/ph6_%s.mat', dirs.mvpa.scratch{xph}, ph6.basename);
save(fname,'ph6','-v7.3');

end%function