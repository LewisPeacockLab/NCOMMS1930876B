function[] = clearmem_mvpa_operation_01(args, dirs)
% step 1: read EPI patterns + regressors/selectors + feature selection
%         decoding setup: decode all timepoints for the testing run

%% ============= UNPACK ARGS.
xph = args.xphase;

%% ============= LOGGING
fprintf('(+) using scratch dir: %s\n', dirs.mvpa.scratch{xph});
%*************** turn on diary to capture analysis output
diary off;
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
%*************** define selectors / zsoring epis

fprintf('\n#####################################################################\n');
fprintf('* STEP_01: load mask + read/zscore epis patterns\n');

nm.info = 'object_name';
[ph1.subj, ph1.nm] = mvpa_ph01_patterns_gz(args, subj, nm, dirs);

%*************** ph1. base filename
ph1.basename = sprintf('%s_%s_%s', args.phase_name{xph}, args.mask_name, args.epi_name); 
ph1.args = args; 

%*************** save phase 1.
fname = sprintf('%s/ph1_%s.mat', dirs.mvpa.scratch{xph}, ph1.basename);
save(fname,'ph1','-v7.3');

%*************** merge subj. structure
subj = ph1.subj; nm = ph1.nm;
summarize(subj)

%% ============= 02: REGRESSORS
%*************** define regressors 
%*************** regressors: shift / norest

fprintf('\n#####################################################################\n');
fprintf('* STEP_02: define regressors + selectors\n');

[ph2.subj, ph2.nm] = mvpa_ph02_regressors(args, subj, nm);

%*************** ph2. base filename
ph2.basename = sprintf('%s_%s_tr%s_blk_%s',...
    ph1.basename, args.regress_type, ...
    num2str(args.shift_TRs), args.rest);

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
subj = ph3.subj;  nm = ph3.nm;
summarize(subj)

%% ============= 04: DECODING SETUP
%*************** change test regressors/selector for all timepoints

fprintf('\n#####################################################################\n');
fprintf('* STEP_04: setup decoding\n');

%*************** combine train & test phase selectors
% nm.train_selector: train selector
% nm.runs{2}: study selector
% set train selector: 1 & test selector: 2

nm.decoding_runs = sprintf('%s_decode', nm.run_xvalid{xph});

fprintf('\n(+) modify test phase selectors for all timepoints: %s\n', nm.decoding_runs);

ori_selector = get_mat(subj,'selector', nm.runs{xph});
n_iterations = numel(unique(ori_selector));

for xiter = 1:n_iterations
    xobj_name = sprintf('%s_%d', nm.decoding_runs, xiter);
    
    xselector = get_mat(subj,'selector', sprintf('%s_%d', nm.run_xvalid{xph}, xiter));
    xselector(ori_selector==xiter) = 2;% testing run
    
    subj = initset_object(subj,'selector', xobj_name, ...
        xselector, 'group_name', nm.decoding_runs);
end

%*************** save subj. structure
ph4.subj = subj; ph4.nm = nm; ph4.args = args;

%*************** base filename
ph4.basename = sprintf('%s_decoding_setup', ph3.basename);

%*************** save phase 4.
fname = sprintf('%s/ph4_%s.mat', dirs.mvpa.scratch{xph}, ph4.basename);
save(fname,'ph4','-v7.3');

fprintf('\n... saved decoding setup: ph4 of %s: %s\n', ...
    args.subject_id, fname);

end
