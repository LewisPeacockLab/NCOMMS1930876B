function[] = single_analysis_auc_localizer(xresults, args, dirs)

xph           = args.xphase;
param         = args.index{xph}.param;
tt_matrix     = args.index{xph}.matrix;
xheader       = args.index{xph}.header;
n_header      = length(xheader);
n_runs        = args.index{xph}.param.n_runs;

%*************** output basename
basename      = args.analysis_basename;
xoutput_dir   = dirs.mvpa.group.auc{xph};

n_subcategory = param.n_subcategory;
n_category    = param.n_category;

%% ============= RANDOM PERMUTATION FOR BASELINE

fprintf('%s baseline roc: SN %s\n', args.phase_name{xph}, args.subject_id);

t_matrix = [];
%*************** shift regressor
for xrun = 1:n_runs
    xunit = getDATA(tt_matrix',xheader, {'run'}, {xrun});
    it_matrix = tt_matrix(:, xunit);
    
    t_matrix = horzcat(t_matrix, ...
        [zeros(size(it_matrix,1), args.shift_TRs) it_matrix(:,1:end-args.shift_TRs)]); %#ok<*AGROW>
end

%*************** reset trials
for xrun = 1:n_runs
    t_trial = unique(getDATA(t_matrix',xheader, {'run'}, {xrun}, findCol(xheader, {'trial'})));
    t_trial = t_trial(~ismember(t_trial, 0));
    n_trials(xrun) = length(t_trial);
end

%*************** reset matrix
t_cate_array  = t_matrix(findCol(xheader, {'category'}), :);
t_spike_array = t_matrix(findCol(xheader, {'spike'}), :);

xmatrix = t_matrix(:,(t_cate_array~=0) & (t_spike_array~=1));
xheader{n_header + 1} = 'class';
xheader{n_header + 2} = 'rand_class';

xcate_array    = xmatrix(findCol(xheader, {'category'}), :);
xsubcate_array = xmatrix(findCol(xheader, {'subcategory'}), :);

if strcmp(args.level, 'category')
    xclass_array = xcate_array;
elseif strcmp(args.level, 'subcategory')
    xclass_array = xsubcate_array + (n_subcategory * (xcate_array-1));
end

xmatrix(findCol(xheader, {'class'}), :) = xclass_array;

%*************** extract param array
for xrun = 1:n_runs
    for xtrial = 1:n_trials(xrun)
        xcate = unique(getDATA(xmatrix', xheader, ...
            {'run','trial','spike'}, {xrun,xtrial,0}, ...
            findCol(xheader, {'category'})));
        xsubcate = unique(getDATA(xmatrix', xheader, ...
            {'run','trial','spike'}, {xrun,xtrial,0}, ...
            findCol(xheader, {'subcategory'})));
        
        if strcmp(args.level, 'category')
            xclass = xcate;
        elseif strcmp(args.level, 'subcategory')
            xclass = xsubcate + (n_subcategory * (xcate-1));
        end
        
        it_unit = xtrial + sum(n_trials(1:xrun-1));
        class_array(it_unit) = xclass;
    end
end

for xrun=1:n_runs
    n_acts(xrun) = length(find(getDATA(xmatrix', xheader,{'run','spike'}, ...
        {xrun,0},findCol(xheader, {'category'}))));
end

%*************** mvpaout
n_target   = size(xresults.iterations(1).acts, 1);
n_iter     = size(xresults.iterations, 2);
for xiter = 1:n_iter
    n_vol(xiter) = size(xresults.iterations(xiter).acts, 2);
end

targ_array = 1:n_target;

xacts = [];
for xiter = 1:n_iter
    %*************** pfmet: based on shifted regressor / (startTR:end)
    if n_acts(xiter) ~= n_vol(xiter)
        xacts = horzcat(xacts, xresults.iterations(xiter).acts(:,1:n_acts(xiter)));
    else
        xacts = horzcat(xacts, xresults.iterations(xiter).acts);
    end
end

%% ============= permutation
tic
for xboots = 1:args.n_iteration
    clear rand_array
    
    if mod(xboots, 10)==0, fprintf('%s.', num2str(xboots)); end
    if mod(xboots, 100)==0, fprintf('\n'); end
    
    rand_matrix = xmatrix; rand_header = xheader;
    
    %***************** random labels matrix
    for xclass = 1:n_target
        xunit    = find(class_array == xclass);
        n_unit   = length(xunit);
        non_targ = targ_array(~ismember(targ_array, xclass));
        
        t_rand   = shuffle(repmat(non_targ, [1, round(n_unit/length(non_targ))]));
        n_rand   = length(t_rand);
        
        if n_rand < n_unit
            for i = 1:(n_unit - n_rand)
                t_shuffle = shuffle(non_targ);
                t_rand(n_rand+i) = t_shuffle(i);
            end
        end
        
        rand_array(xunit) = t_rand(1:length(xunit));
    end
    
    for xrun = 1:n_runs
        for xtrial = 1:n_trials(xrun)
            xunit = getDATA(xmatrix', xheader, ...
                {'run','trial'}, {xrun, xtrial});
            xtarg_unit = xtrial + sum(n_trials(1:xrun-1));
            
            rand_matrix(findCol(rand_header, {'rand_class'}), xunit) = ...
                rand_array(xtarg_unit);
        end
    end
    
    %***************** random labels for baseline permutation test
    xdesireds = rand_matrix(findCol(rand_header, {'rand_class'}), :);
    
    %*************** ROC inputs: perfcurve(labels, scores, posclass)
    % scores: classifier predictions: (evidence)
    % labels: given true class labels: xdesireds
    % posclass: positive class label: (xcate)
    % T: (mean, lower bound, upper bound)
    for xclass = 1:n_target
        
        xlabels   = xdesireds;
        xscores   = xacts(xclass, :);
        xposclass = xclass;
        
        %*************** perfcurve
        [~, ~, ~, AUC] = perfcurve(xlabels, xscores, xposclass);
        
        baseline{xclass}(xboots) = AUC; %#ok<*NASGU>
    end
end

fprintf('\n')

%*************** save mvpa_ROC
fname = sprintf('%s/auc_baseline_%s_%siters_%s.mat', ...
    dirs.mvpa.group.auc_baseline{xph}, basename, num2str(args.n_iteration), args.subject_id);
save(fname, 'baseline', '-v7.3')

fprintf('...%s: %3.4f seconds\n\n', args.subject_id, toc)

end