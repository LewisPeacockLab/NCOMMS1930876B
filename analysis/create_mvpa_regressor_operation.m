function[regs] = create_mvpa_regressor_operation_four(args, dirs)
%% =============== LOAD DESIGN INDEX
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

n_volumes = param.n_cat_trim_vols;%size(args.index{xph}.matrix,2);%

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
    xunit = find(getDATA(regs.matrix', regs.header, ...
        {'condition','manipulation'}, ...
        {xcond, 1}));
    
    regs.regressors(it, xunit)   = 1;
    regs.regressors_index(xunit) = xcond;
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