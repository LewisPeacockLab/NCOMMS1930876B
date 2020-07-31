function[regs] = create_spm_regressor_study(args, dirs)
%%%%%%%%%%%%%%%%%% creating SPM regressors
%%%%%%%%%%%%%%%%%% Multiple condition: names, onsets, durations
%%%%%%%%%%%%%%%%%% Multiple regressors: motion correction

%% =============== LOAD DESIGN INDEX
% cf.: the matrix is for trimmed TR, but the epi is not trimmed
% So, I need to use not trimmed regressor for GLM analyses

xph           = 2;
xheader       = args.index{xph}.header;
xmatrix       = args.index{xph}.matrix_notrim';
xparam        = args.index{xph}.param;

%% *************** unpack parameters
n_category    = xparam.n_category;
n_subcategory = xparam.n_subcategory;
n_item        = xparam.n_item;
n_runs        = xparam.n_runs;% selected run
n_conditions  = length(xparam.conds_names);
it_class      = args.it_class;
xlevel        = ismember(args.levels, args.level);
n_class       = length(it_class{xlevel});

for xcate = 1:n_category
    for xsubcate = 1:n_subcategory
        xunit = getDATA(xmatrix, xheader, ...
            {'category', 'subcategory'}, {xcate, xsubcate});
        it_subcate = xsubcate + (n_subcategory * (xcate-1));
        xmatrix(xunit, findCol(xheader, {'subcategory'})) = it_subcate;
        
        for xitem = 1:n_item
            xunit = getDATA(xmatrix, xheader, ...
                {'category', 'subcategory','item'}, {xcate, xsubcate, xitem});
            it_item = xitem + (n_item * (xsubcate-1)) + ((n_item * n_subcategory) * (xcate-1));
            xmatrix(xunit, findCol(xheader, {'item'})) = it_item;
        end
    end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% GENERATE REGRESSOR FOR SPM GLM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% names, onsets, durations

n_regs    = n_conditions * n_class;
%***************** onsets & durations
names     = cell(1, n_regs);
onsets    = cell(1, n_regs);
durations = cell(1, n_regs);

for xcond = 1:n_conditions
    for xclass = 1:n_class
        xreg = xclass + (n_class * (xcond-1));
        names{xreg} = sprintf('cond%d_class%d', xcond, xclass);
        
        for xrun = 1:n_runs
            it_trials = unique(getDATA(xmatrix, xheader, ...
                {'run', 'condition', args.levels(xlevel)}, ...
                {xrun, xcond, xclass}, findCol(xheader, {'trial'})));
            
            for xtrial = it_trials'
                xunit = getDATA(xmatrix, xheader, ...
                    {'run', 'trial'}, {xrun, xtrial});
                
                xonset = xmatrix(xunit, findCol(xheader, {'it_volume'}));
                
                onsets{xreg}    = horzcat(onsets{xreg}, xonset(1));
                durations{xreg} = horzcat(durations{xreg}, 0);
            end
        end
    end
end

%***************** save in the structure
regs.names     = names;
regs.onsets    = onsets;
regs.durations = durations;

%***************** SAVE REGRESSOR
save(fullfile(dirs.uni.regressor, ...
    sprintf('study_spm_regressor_%s.mat', args.level)),'names','onsets','durations','-v7.3');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% MOTION REGRESSOR FOR SPM GLM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n_all_runs    = 6;
t_motion      = load(fullfile(dirs.motion, 'motion_study.txt'));
selected_runs = unique(xmatrix(:, findCol(xheader, {'it_run'})));
motion_matrix = t_motion(:, 1:6);
n_vols        = args.index{xph}.param.n_volumes;
vol_index     = ones(1, n_all_runs).* double(n_vols(1));
vol_index(selected_runs(~ismember(n_vols, n_vols(1)))) = ...
    n_vols(~ismember(n_vols, n_vols(1)));

xmotion = []; xrun_regs = [];
for xit = 1:length(selected_runs)
    xrun = selected_runs(xit);
    
    xunit   = (1:n_vols(xit)) + sum(n_vols(1:xit-1));
    yunit   = (1:vol_index(xrun)) + sum(vol_index(1:xrun-1));
    xmotion = vertcat(xmotion, motion_matrix(yunit, :)); %#ok<*AGROW>
    
    t_run_regs = zeros(length(find(xunit)), n_runs);
    t_run_regs(:, xit) = 1;
    
    xrun_regs  = vertcat(xrun_regs, t_run_regs);
end

%***************** save in the structure
regs.multi_regress = horzcat(xmotion, xrun_regs);

%***************** SAVE REGRESSOR
dlmwrite(fullfile(dirs.uni.regressor, 'study_multi_regressor.txt'),...
    regs.multi_regress);

end