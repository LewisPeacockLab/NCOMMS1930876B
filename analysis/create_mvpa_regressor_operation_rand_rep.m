function[regs] = create_mvpa_regressor_operation_new(args, dirs)
%% =============== LOAD DESIGN INDEX
% 4 operations: randomly select 1/2 trials from each replace

xph            = 3;
xheader        = args.index{xph}.header;%from study
xmatrix        = args.index{xph}.matrix;
param          = args.index{xph}.param;
n_header       = length(xheader);

%***************** unpack parameters
condition_name = param.operation_names;
n_operation    = length(param.operation_names);
n_volumes      = param.n_cat_trim_vols;%size(args.index{xph}.matrix,2);%
n_runs         = param.n_runs;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% GENERATE REGRESSOR FOR MVPA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%***************** 1/2 trial random selection for replaces
xheader{n_header + 1} = 'selected';%selected volume
xmatrix(findCol(xheader, {'selected'}), :) = 0;

it_conds = [2 3]; 
for xrun = 1:n_runs
    clear tt
    %***************** randomly pick trials per run
    for it = 1:length(it_conds)
        xcond = it_conds(it);
        xtrial = unique(getDATA(xmatrix', xheader, ...
            {'run','condition','manipulation','spike'}, ...
            {xrun, xcond, 1, 0}, findCol(xheader, {'trial'})));
        
        tt{it} = shuffle(xtrial)'; %#ok<*AGROW>
    end
    
    min_num = min(length(tt{1}), length(tt{2}));
    
    for it = 1:length(it_conds)
        sel_trials{it} = tt{it}(1:round(min_num/2));
    end
    
    %***************** setup selected volume
    for it = 1:length(it_conds)
        for xtrial = sel_trials{it}
            xcond = it_conds(it);
            xunit = getDATA(xmatrix', xheader, ...
                {'run','trial','condition','manipulation','spike'}, ...
                {xrun, xtrial, xcond, 1, 0});
            
            xmatrix(findCol(xheader, {'selected'}), xunit) = 1;
        end
    end
end

%% ================ REGRESSOR/SELECTOR FOR SHIFTED EPI
%***************** selectors: runs
%***************** trial: 0_FIX, each trial: stim (6 tr) + manipulation (6tr) + iti (5:9 tr)
%***************** class: 4 operations
%***************** remove spike points
%***************** presentation: 1_stim, 2_manipulation, 0_fixation
%***************** manipulation: 0_stim, 1_instruction + fixation (operation)

regs.selectors = xmatrix(findCol(xheader, {'run'}), :);

%************* operations
regs.regressor_name   = condition_name;
regs.regressors       = zeros(n_operation, n_volumes);
regs.regressors_index = zeros(1, n_volumes);

for xcond = 1:n_operation
    if xcond == 2
        xunit = find(getDATA(xmatrix', xheader, ...
            {'operation','manipulation','spike','selected'}, ...
            {xcond, 1, 0, 1}));
    else
        xunit = find(getDATA(xmatrix', xheader, ...
            {'operation','manipulation','spike'}, ...
            {xcond, 1, 0}));
    end
    
    regs.regressors(xcond, xunit) = 1;
    regs.regressors_index(xunit)  = xcond;
end

%% =========== SAVE REGRESSOR
regs.matrix = xmatrix;
regs.header = xheader;

save(fullfile(dirs.mvpa.regressor{xph}, ...
    sprintf('operation_classification_four_regressor_%s_%dtr_%s.mat', ...
    args.train_regress_type, args.shift_TRs, args.rest)),'regs');

end