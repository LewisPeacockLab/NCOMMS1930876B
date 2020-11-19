function[subj, nm] = mvpa_ph02_regressors(args, subj, nm)

    %---------------------------------------------------------------------------
    %*************** define regressors
    %*************** shift regressors: ShiftTR or Convolution
    %---------------------------------------------------------------------------
    
    xph        = args.xphase;
    regressors = args.regs{xph};
    
    %% ============= STRUCTURE REGRESSORS: subj.regressors{1}
    fprintf('\n(+) structrue regressors & selector: phase_%d\n', xph);
    fprintf('... # of regressors: %d\n', size(regressors.regressors,1));
    disp(regressors.regressor_name)
    fprintf('... # of selectors: %d\n', numel(unique(regressors.selectors)));

    nm.conds{xph} = sprintf('%s_conds', args.phase_name{xph});% regressors object name
    
    subj = initset_object(subj,'regressors', nm.conds{xph}, regressors.regressors);    
    subj = set_objfield(subj,'regressors', nm.conds{xph},'condnames', regressors.regressor_name);
    
    %% ============= SHIFTING THE REGRESSORS ALONG: subj.regressors{2}
    fprintf('\n(+) shifting regressors: %d TR\n', args.shift_TRs);
    
    subj = shift_regressors(subj, nm.conds{xph}, nm.runs{xph}, args.shift_TRs);% shift the regressor
    
    nm.conds_sh{xph} = sprintf('%s_sh%d', nm.conds{xph}, args.shift_TRs);% subj.regressors{2}.name
    
    % remove excessive movement timepoints
    xspike = find(args.regs{xph}.matrix(findCol(args.regs{xph}.header, {'spike'}), :));
    if ~isempty(xspike)
        xregs = get_mat(subj, 'regressors', nm.conds_sh{xph});
        xregs(:, xspike) = 0;
        subj = set_mat(subj, 'regressors', nm.conds_sh{xph}, xregs);
    end
    
    %% ============= REMOVE NO-USING TIME POINTS: subj.selectors{2}
    %*************** getting rid of unused trs (i.e., regressor==0)
    %*************** defined unused trs when making regressors
    subj                    = create_norest_sel(subj, nm.conds_sh{xph});% creating norest selector object
    nm.conds_sh_norest{xph} = sprintf('%s_norest', nm.conds_sh{xph});% norest selector
    
    fprintf('\n(+) removing norest from regressors\n ... regressors: %s\n ... selectors: %s\n\n', subj.r, subj.s);
    
    summarize(subj,'objtype','regressors','display_groups',false)
    summarize(subj,'objtype','selector','display_groups',false)
end
