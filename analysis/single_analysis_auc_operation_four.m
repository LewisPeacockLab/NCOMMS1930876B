function[] = single_analysis_auc_operation_four(xresults, args, dirs)

xph     = args.xphase;
xmatrix = args.index{xph}.matrix;
xheader = args.index{xph}.header;
param   = args.index{xph}.param;

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

%% ============= 1ST LEVEL

fprintf('%s roc: %s\n', args.phase_name{xph}, args.subject_id);

n_target = size(xresults.iterations(1).acts, 1);
n_iter   = size(xresults.iterations, 2);

if n_condition ~= n_target
    fprintf('!!!!! the number of targets are wrong!\n')
end

for xiter = 1:n_iter
    n_vol(xiter) = size(xresults.iterations(xiter).acts, 2);
end

%***************** reset regs: shifting TR
xregs.runs         = args.regs{xph}.selectors;
xregs.operation_sh = zeros(1, sum(n_vol));% only operation window

for it = 1:n_condition
    xcond = it_conds(it);
    xunit = getDATA(xmatrix', xheader, {'condition'}, {xcond});
    
    xregs.timecourse(xunit) = xcond;
    
    xunit = find(getDATA(xmatrix', xheader, ...
        {'condition','manipulation','spike'}, {xcond,1,0}));
    
    xregs.operation_sh(xunit + args.shift_TRs) = it;
end

%*************** extract MVPA performance
acts_array = [];
for xiter = 1:n_iter
    acts_array = horzcat(acts_array, xresults.iterations(xiter).acts); %#ok<*AGROW>
end

xunit     = (xregs.operation_sh ~=0);
xdesireds = xregs.operation_sh(xunit);
xacts     = acts_array(:,xunit);

%*************** ROC inputs: perfcurve(labels, scores, posclass)
% scores: classifier predictions: (evidence for that class)
% labels: given true class labels: xdesireds
% posclass: positive class label: (xcate)
% T: (mean, lower bound, upper bound)

for xcate = 1:n_target
    
    xlabels   = xdesireds;
    xscores   = xacts(xcate, :);
    xposclass = xcate;
    
    %*************** perfcurve    
    [~, ~, ~, AUC] = perfcurve(xlabels, xscores, xposclass);
    
    xroc_out{xcate} = AUC;
    
end

%% ============= RANDOM PERMUTATION FOR BASELINE
if args.permutation{3}    
    
    targ_array = 1:n_target;

    %% ============= EACH SUBJECT
    tic
    
    fprintf('%s baseline roc: %s\n', args.phase_name{xph}, args.subject_id);
    
    fprintf('...%s\n', args.subject_id);
    
    n_target = size(xresults.iterations(1).acts, 1);
    n_iter   = size(xresults.iterations, 2);
    for xiter = 1:n_iter
        n_vol(xiter) = size(xresults.iterations(xiter).acts, 2);
    end
    
    %*************** extract param array    
    for xrun = 1:n_iter
        t_trial = unique(getDATA(xmatrix', xheader, {'run'}, {xrun}, ...
            findCol(xheader, {'trial'})));
        n_trial(xrun) = length(t_trial(~ismember(t_trial, 0)));
    end
    
    for xrun = 1:n_iter
        for xtrial = 1:n_trial(xrun)
            xcond = unique(getDATA(xmatrix', xheader, ...
                {'run','trial'}, {xrun,xtrial}, ...
                findCol(xheader, {'condition'})));
            
            it_unit = xtrial + sum(n_trial(1:xrun-1));
            it      = find(it_conds==xcond);
            
            if sum(it)
                cond_array(it_unit) = it;
            else
                cond_array(it_unit) = 0;
            end
        end
    end
    
    %% ============= permutation
    for xboots = 1:args.n_iteration
        clear rand_array
        
        if mod(xboots, 10)==0, fprintf('%s.', num2str(xboots)); end
        if mod(xboots, 100)==0, fprintf('\n'); end
        
        rand_matrix = xmatrix;
        rand_header = [xheader 'rand_condition'];
        
        %***************** random labels matrix
        rand_array = zeros(1, size(cond_array, 2));
        for xcond = 1:n_target
            xunit    = find(cond_array == xcond);
            n_unit   = length(xunit);
            non_targ = targ_array(~ismember(targ_array, xcond));
            
            rand_xx  = shuffle(repmat(non_targ, [1, round(n_unit/(length(non_targ)))]));
            
            if length(xunit)==length(rand_xx)
                rand_array(xunit) = rand_xx;
            else
                tt = shuffle(non_targ);
                rand_xx = [rand_xx tt];
                rand_array(xunit) = rand_xx(1:length(xunit));
            end
        end
        
        for xrun = 1:n_iter
            for xtrial = 1:n_trial(xrun)
                xcond = unique(getDATA(xmatrix', xheader, ...
                    {'run','trial'}, {xrun,xtrial}, ...
                    findCol(xheader, {'condition'})));
                
                if ismember(xcond, it_conds)
                    
                    xunit = getDATA(xmatrix', xheader, ...
                        {'run','trial'}, {xrun, xtrial});
                    xtarg_unit = xtrial + sum(n_trial(1:xrun-1));
                    
                    rand_matrix(findCol(rand_header, {'rand_condition'}), xunit) = ...
                        rand_array(xtarg_unit);
                end
            end
        end
        
        %***************** random labels for baseline permutation test
        rand_operation_sh = zeros(1, sum(n_vol));
        
        for it = 1:n_condition            
            if args.reset_regs{xph}
                xunit = find(getDATA(rand_matrix', rand_header, ...
                    {'rand_condition','presentation','spike'}, ...
                    {it, 2, 0}));
            else
                xunit = find(getDATA(rand_matrix', rand_header, ...
                    {'rand_condition','manipulation','spike'}, ...
                    {it, 1, 0}));
            end
            
            rand_operation_sh(xunit + args.shift_TRs) = it;
        end
        
        xunit     = (rand_operation_sh ~=0);
        xdesireds = rand_operation_sh(xunit);
        xacts     = acts_array(:,xunit);
        
        %*************** ROC inputs: perfcurve(labels, scores, posclass)
        % scores: classifier predictions: (evidence)
        % labels: given true class labels: xdesireds
        % posclass: positive class label: (xcate)
        % T: (mean, lower bound, upper bound)
        for xcond = 1:n_target
            
            xlabels   = xdesireds;
            xscores   = xacts(xcond, :);
            xposclass = xcond;
            
            %*************** perfcurve
            [~, ~, ~, AUC] = perfcurve(xlabels, xscores, xposclass);
            
            baseline{xcond}(xboots) = AUC; %#ok<*NASGU>
        end
    end
    
    fprintf('\n')
    
    %*************** save mvpa_ROC
    fname = sprintf('%s/auc_baseline_four_oper_%siters_%s.mat', ...
        dirs.mvpa.group.auc_baseline{xph}, num2str(args.n_iteration), args.subject_id);
    save(fname, 'baseline', 'xroc_out', '-v7.3')
    
    fprintf('...%s: %3.4f seconds\n\n', args.subject_id, toc)
    
end

end