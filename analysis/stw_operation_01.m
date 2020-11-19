function[] = stw_operation_01(args, dirs)
% step 1: load EPI patterns + regressors/selectors
%         + feature selection + classification (penalty: 0.5)

%% ============= UNPACK ARGS
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

%% ============= 01: EPI PATTERNS
%*************** load mask + read in epis
%*************** zsoring epis
%*************** option: read epis in mask | wholebrain
%*************** define: args.wholebrain '0 | 1'

fprintf('\n#####################################################################\n');
fprintf('* STEP_01: load mask + epis patterns\n');

%*************** study ph1. base filename
ph1.basename = sprintf('%s_%s_%s', args.phase_name{xph}, args.mask_name, args.epi_name); 

src_args  = args;
src_args.reset_regs{xph} = 0;
src_args.stw{xph}        = 0;
src_dirs  = setup_directory(dirs, src_args);
src_fname = sprintf('%s/ph1_%s.mat', src_dirs.mvpa.scratch{xph}, ph1.basename);
load(src_fname);

%% ============= SETUP REGRESSORS
%*************** define regressors + selectors
%*************** regressors: shift or beta
xheader        = args.index{xph}.header;%from study
xmatrix        = args.index{xph}.matrix;
xparam         = args.index{xph}.param;

%*************** unpack parameters
condition_name = xparam.conds_names;
n_condition    = length(condition_name);
n_volumes      = xparam.n_cat_trim_vols;%size(args.index{xph}.matrix,2);%
n_runs         = xparam.n_runs;
n_trials       = xparam.n_trials;

%*************** sliding time window param
stw_win = floor(args.stw_win/2);%befor/after xtr in time-window for stw
it_wins = 1:(args.tc_tr_disp+1);

%*************** regressor
args.regs{xph}.selectors      = xmatrix(findCol(xheader, {'run'}), :);
args.regs{xph}.regressor_name = condition_name;

%% ============= STW 
%*************** regressors
for xwin = it_wins
    tic

    fprintf('\n#####################################################################\n');
    fprintf('(+) window %s out of %s: ', num2str(xwin), num2str(length(it_wins)));
    fprintf('\n#####################################################################\n\n');
    
    %*************** ph2. base filename
    ph2.basename = sprintf('win%s_%s_%s_tr%s_blk_%s',...
        num2str(xwin), ph1.basename, args.regress_type, num2str(args.shift_TRs), args.rest);
    
    %*************** ph3. base filename
    ph3.basename = sprintf('%s_featsel_thresh%s', ph2.basename, num2str(args.featSelThresh));
    
    %*************** merge subj. structure
    clear subj
    subj = ph1.subj; nm = ph1.nm;
    
    %% ============= 02: REGRESSORS
    clear regressors
    
    fprintf('* STEP_02: define regressors + selectors\n');
    
    %*************** selector
    regressors = zeros(n_condition, n_volumes);
    
    %*************** regressor
    for xrun = 1:n_runs
        for xtrial = 1:n_trials(xrun)
            tunit  = find(getDATA(xmatrix', xheader, ...
                {'run','trial'}, {xrun, xtrial}));
            
            xcond  = xmatrix(findCol(xheader, 'condition'), tunit(1));
            xunit  = (tunit(1) + (xwin-1) - stw_win):(tunit(1) + (xwin-1) + stw_win);
            xspike = xmatrix(findCol(xheader, 'spike'), xunit);
            xreg   = xunit(~xspike);
            
            regressors(xcond, xreg) = 1;
        end
    end
    
    args.regs{xph}.regressors = regressors;
    
    %*************** run ph2
    [ph2.subj, ph2.nm] = mvpa_ph02_regressors(args, subj, nm);
    
    %*************** merge subj. structure
    subj = ph2.subj; nm = ph2.nm;
    summarize(subj)
    
    %% ============= 03: FEATURE SELECTION
    %*************** creating cross-validation indices
    %*************** feature selection anova.
    %*************** options: shiftTRs, peakWindow, featSelThresh
    
    fprintf('\n* STEP_03:feature selection\n');
    
    [ph3.subj, ph3.args, ph3.nm] = mvpa_ph03_featselection(args, subj, nm, dirs);
    
    %*************** save phase 3.
    fname = sprintf('%s/ph3_%s.mat', dirs.mvpa.scratch{xph}, ph3.basename);
    save(fname,'ph3','-v7.3');
    
    %% ============= 04: MVPA CLASSIFICATION
    %*************** merge subj. structure
    subj = ph3.subj;  nm = ph3.nm;
    summarize(subj)
    
    %*************** definde penalty + classification
    fprintf('\n(+) subject:%s classifier: %s, max penalty: %s\n', ...
        args.subject_id, args.classifier, num2str(args.xpenalty));
    fprintf('... classification %s phase \n', args.phase_name{xph});
    
    [ph4.args, ph4.subj, ph4.results] = mvpa_ph04_classification(args, subj, nm);
    
    %*************** basename phase 4.
    ph4.basename  = sprintf('%s_penalty%s', ph3.basename, num2str(args.xpenalty));
    
    %*************** save phase 4.
    fname = sprintf('%s/%s.mat', dirs.mvpa.output{xph}, ph4.basename);
    save(fname,'ph4','-v7.3');
    
    fprintf('\n... saved results_penalty%s of %s: %s\n', ...
        num2str(args.xpenalty), args.subject_id, fname);
    
    summarize(subj)
    
    %*************** save phase 4.results
    xresults = ph4.results;
    fname = sprintf('%s/%s_result.mat', dirs.mvpa.output{xph}, ph4.basename);
    save(fname,'xresults','-v7.3');
    
    toc
    
end
end
