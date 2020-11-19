function[] = clearmem_localizer_classification_01(args, dirs)
% step 1: read EPI patterns + regressors/selectors + feature selection

%% ============= UNPACK ARGS.
xph = args.xphase;

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

%% ============= Initializing subj. structure:start by creating an empty subj structure
% summarize(subj): summarize all info in subj structure
% get_object/set_object/set_objfield/set_objsubfield
% get_mat/set_mat

subj = init_subj(args.experiment, args.subject_id);%identifier of the subj
fprintf('\n(+) %s phase data\n\n', args.phase_name{xph});

%% ============= 01: EPI PATTERNS
%*************** load mask + read in epis
%*************** define selectors + zsoring epis

fprintf('\n#####################################################################\n');
fprintf('* STEP_01: load mask + read/zscore epis patterns\n');

nm.info = 'object_name';
[ph1.subj, ph1.nm] = mvpa_ph01_patterns_gz(args, subj, nm, dirs);

%*************** ph1. base filename
ph1.basename = sprintf('%s_%s_zscored_%s', args.phase_name{xph}, args.mask_name, args.epi_name); 

%*************** save phase 1.
fname = sprintf('%s/ph1_%s.mat', dirs.mvpa.scratch{xph}, ph1.basename);
ph1.args = args; save(fname,'ph1','-v7.3');

%*************** merge subj. structure
subj = ph1.subj; nm = ph1.nm;
summarize(subj)

%% ============= 02: REGRESSORS
%*************** define regressors
%*************** regressors: shift

fprintf('\n#####################################################################\n');
fprintf('* STEP_02: define regressors + selectors\n');

[ph2.subj, ph2.nm] = mvpa_ph02_regressors(args, subj, nm);

%*************** ph2. base filename
ph2.basename = sprintf('%s_%s_%s%dtr_blk_%s',...
    ph1.basename, args.level, args.regress_type, ...
    args.shift_TRs, args.rest);

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
ph3.basename = sprintf('%s_featsel_thresh%s', ph2.basename, num2str(args.featSelThresh));

%*************** save phase 3.
fname = sprintf('%s/ph3_%s.mat', dirs.mvpa.scratch{xph}, ph3.basename);
save(fname,'ph3','-v7.3');

%*************** merge subj. structure
subj = ph3.subj; 
summarize(subj)

end
