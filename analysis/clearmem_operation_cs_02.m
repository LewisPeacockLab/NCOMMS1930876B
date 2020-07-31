function[] = clearmem_operation_cs_02(args, dirs)

xph = args.xphase;

%*************** L2logreg
class_args.train_funct_name = 'train_L2_RLR';
class_args.test_funct_name  = 'test_L2_RLR';

%*************** getting rid of 1-of-n warnings using propval
cv_args.perfmet_args.ignore_1ofn = true;

%% ============= LOAD SUBJ: FEATURE SELECTED
fprintf('#############################################\n')
fprintf('[2]: penalty check: cross-subject classification\n')
fprintf('#############################################\n\n')

xbasename = sprintf('cs_operation_%s_%s_%s_fsthres%s_n%s', ...
    args.cs_type, args.epi_name, args.mask_name, ...
    args.featSelThresh, num2str(args.n_sub));

%*************** loading ph2 
xfname = fullfile(dirs.mvpa.cs.scratch, ...
    sprintf('ph2_%s.mat', xbasename));

fprintf('... loading ph2 concatenate data ...\n')
load(xfname);%'cs_ph2'
subj = cs_ph2.subj;
% nm   = cs_ph2.nm;

%% ============= PH3: MVPA CLASSIFICATION & PENALTY CHECK
nm.pat_z      = 'cross_subj_operation_patterns_z';
nm.reg        = 'cross_subj_operation_regs';
nm.reg_sh     = sprintf('%s_sh%d', nm.reg, args.shift_TRs);% subj.regressors{2}.name
nm.sel        = 'cross_subj_operation_selector';
nm.sel_norest = sprintf('%s_norest', nm.reg_sh);% norest selector
nm.sel_xval   = sprintf('%s_norest_xval', nm.sel);
xmask_name    = sprintf('%s_thresh%s', nm.pat_z, args.featSelThresh);

fprintf('\n(+) classification: cross-validation (leave 1-out)\n');

tic
if args.penalty_step == 1
    clear xpenalty
    xpenalty = args.penalty(args.it_penalty);
    class_args.lambda = xpenalty;
    
    [~, xresults] = cross_validation(subj, nm.pats_z, nm.reg_sh,...
        nm.sel_xval, xmask_name, class_args, cv_args);
    
    p_fname = fullfile(dirs.mvpa.cs.penalty, ...
        sprintf('perf_%s_p%s.mat', xbasename, num2str(xpenalty)));
    
    pen_check.penalty    = xpenalty;
    pen_check.total_perf = xresults.total_perf;
    
    save(p_fname, 'pen_check')
elseif args.penalty_step == 2
    %*************** checking the max penalty from the step 1
    for i = 1:length(args.penalty)
        clear p_fname
        xpenalty = args.penalty(i);
        
        p_fname = fullfile(dirs.mvpa.cs.penalty, ...
            sprintf('perf_%s_p%s.mat', xbasename, num2str(xpenalty)));
        
        load(p_fname);%pen_check
        
        cs_pen_check.penalty(i)    = pen_check.penalty;
        cs_pen_check.total_perf(i) = pen_check.total_perf;
    end
    
%     [xacc, whichmax] = max(pen_check.performance);
%     fprintf('... max penalty from step1: max acc: %4.4f\n', xacc);
%     pen_check.penalty(whichmax);
    
end
toc

end