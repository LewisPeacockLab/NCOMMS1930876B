function[regs] = create_mvpa_regressor_operation(args, dirs)
%% =============== LOAD DESIGN INDEX

xph            = 3;
regs.header    = args.index{xph}.header;%from study
regs.matrix    = args.index{xph}.matrix;
param          = args.index{xph}.param;

%***************** unpack parameters
condition_name = param.conds_names;
n_condition    = length(condition_name);
n_operation    = length(param.operation_names);
n_volumes      = param.n_cat_trim_vols;%size(args.index{xph}.matrix,2);%

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% GENERATE REGRESSOR FOR MVPA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%================ REGRESSOR/SELECTOR FOR SHIFTED EPI
%***************** selectors: runs
%***************** trial: 0_FIX, each trial: stim (6 tr) + manipulation (6tr) + iti (5:9 tr)
%***************** class: 5 operations + perception + rest | 5 operations
%***************** remove spike points
%***************** presentation: 1_stim, 2_manipulation, 0_fixation
%***************** manipulation: 0_stim, 1_instruction + fixation (operation)
regs.selectors = regs.matrix(findCol(regs.header, {'run'}), :);

%************* operations
regs.regressor_name   = condition_name;
regs.regressors       = zeros(n_condition, n_volumes);
regs.regressors_index = zeros(1, n_volumes);

if args.four_oper_regress
    for xcond = 1:n_operation
        xunit = find(getDATA(regs.matrix', regs.header, ...
            {'n_operation','manipulation','spike'}, ...
            {xcond, 1, 0}));
        
        regs.regressors(xcond, xunit) = 1;
        regs.regressors_index(xunit)  = xcond;
    end
else
    for xcond = 1:n_condition
        if args.reset_6tr{xph}
            xunit = find(getDATA(regs.matrix', regs.header, ...
                {'condition','presentation'}, ...
                {xcond, 2}));
        else
            xunit = find(getDATA(regs.matrix', regs.header, ...
                {'condition','manipulation'}, ...
                {xcond, 1}));
        end
        
        regs.regressors(xcond, xunit) = 1;
        regs.regressors_index(xunit)  = xcond;
    end
end

%% =========== rest
if strcmp(args.rest, 'rest')
    %************* presentation
    xregs = n_condition + 1;% 6: presentation
    regs.regressor_name{xregs} = 'perception';
    regs.regressors(xregs,:)   = zeros(1, n_volumes);
    
    xunit = find(getDATA(regs.matrix', regs.header, ...
        {'presentation','spike'}, {1, 0}));
    
    regs.regressors(xregs, xunit) = 1;
    regs.regressors_index(xunit)  = xregs;

    %************* rest
    xregs = length(regs.regressor_name) + 1;
    regs.regressor_name{xregs} = 'rest';
    regs.regressors(xregs,:)   = zeros(1, n_volumes);
    
    sum_regressors = sum(regs.regressors);
    regs.regressors(xregs, sum_regressors==0) = 1;
    regs.regressors_index(sum_regressors==0)  = xregs;
    
    %************* remove spike
    xunit = find(getDATA(regs.matrix', regs.header, {'spike'}, {1}));
    
    regs.regressors(xregs, xunit) = 0;
    regs.regressors_index(xunit)  = 0;
end

%% =========== SAVE REGRESSOR
if args.four_oper_regress
    save(fullfile(dirs.mvpa.regressor{xph}, ...
        sprintf('operation_classification_four_regressor_%s_%dtr_%s.mat', ...
        args.train_regress_type, args.shift_TRs, args.rest)),'regs');
else
    save(fullfile(dirs.mvpa.regressor{xph}, ...
        sprintf('operation_classification_regressor_%s_%dtr_%s.mat', ...
        args.train_regress_type, args.shift_TRs, args.rest)),'regs');
end

end