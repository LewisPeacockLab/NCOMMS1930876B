function[regs] = create_mvpa_regressor_operation_four(args, dirs)
%% =============== LOAD DESIGN INDEX
clear regs

xph            = args.xphase;
regs.header    = args.index{xph}.header;%from study
regs.matrix    = args.index{xph}.matrix;
param          = args.index{xph}.param;

%***************** unpack parameters
if args.four_oper_regress
    it_conds = [1 2 4 5];
else
    it_conds = 1:5;
end

n_condition    = length(it_conds);
condition_name = cell(1, n_condition);

for it = 1:n_condition
    xcond = it_conds(it);
    condition_name{it} = param.conds_names{xcond};
end

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

for it = 1:n_condition
    xcond = it_conds(it);
    if args.reset_6tr{xph}
        xunit = find(getDATA(regs.matrix', regs.header, ...
            {'condition','presentation'}, ...
            {xcond, 2}));
    else
        xunit = find(getDATA(regs.matrix', regs.header, ...
            {'condition','manipulation'}, ...
            {xcond, 1}));
    end
    
    regs.regressors(it, xunit)   = 1;
    regs.regressors_index(xunit) = xcond;
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
if ~isdir(dirs.mvpa.regressor{xph})
    mkdir(dirs.mvpa.regressor{xph});
end

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