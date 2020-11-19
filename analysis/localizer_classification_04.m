function[] = clearmem_localizer_classification_04(args, dirs)
% step 4: AUC 2nd

%% ============= UNPACK ARGS.
xph          = args.xphase;
param        = args.index{xph}.param;
subject_list = args.subject_list;
xsub_grp     = args.filtered_subs;

%*************** subject id num
for xsub = 1:args.n_sub
    sub_id(xsub) = str2double(subject_list(xsub).name(end-2:end));
end

%*************** filtered subject id num
for it = 1:length(xsub_grp)
    xsub = args.filtered_subs(it);
    sub_id_filtered(it) = str2double(subject_list(xsub).name(end-2:end));
end

%*************** output basename
basename      = args.analysis_basename;
xoutput_dir   = dirs.mvpa.group.auc{xph};

n_subcategory = param.n_subcategory;
n_category    = param.n_category;
it_categories = 1:n_category;

if strcmp(args.level, 'subcategory') && args.class_selecting
    it_categories = args.selected_category;
end

if strcmp(args.level, 'category')
    %*************** create tables
    xcate_name = param.category_name;
    n_class    = n_category;
elseif strcmp(args.level, 'subcategory')
    %*************** create tables
    for xcate = 1:length(it_categories)
        itcate = it_categories(xcate);
        for xsubcate = 1:n_subcategory
            xcate_name{xsubcate + n_subcategory*(xcate-1)} = ...
                param.subcategory_name{itcate}{xsubcate};
        end
    end
    
    n_class = n_category * n_subcategory;
end

if strcmp(args.rest, 'rest'), xcate_name{end + 1} = 'rest'; end

n_target    = length(xcate_name);
xsub_groups = args.filtered_subs;

xcond_color  = args.cond_color;
xcate_color  = args.cate_color;
xbase_color  = args.base_color;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%============= save MVPA results in a group structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ============= SETUP FILE NAMES
%*************** base filename
ph1.basename = sprintf('%s_%s_zscored_%s', args.phase_name{xph}, args.mask_name, args.epi_name);
ph2.basename = sprintf('%s_%s_%s%dtr_blk_%s',...
    ph1.basename, args.level, args.regress_type, ...
    args.shift_TRs, args.rest);
ph3.basename = sprintf('%s_featsel_thresh%s', ph2.basename, num2str(args.featSelThresh));
class_basename  = sprintf('classified_%s_%s', ph3.basename, args.classifier);

%% ============= LOAD 1ST LEVEL
for xsub = xsub_grp
    %*************** setup subject & directories
    args.subject_id = subject_list(xsub).name;
    dirs            = setup_directory(dirs, args);
    
    %% ============= LOAD LOCALIZER PENALTY
    if ~(args.fixed_penalty{xph})
        %*************** load phase 6.
        fprintf('\n\n#######################################\n')
        fprintf(sprintf('... loading penalty from localizer: s_%s_%s\n', num2str(xsub), args.subject_id));
        
        penalty_rest = args.rest;
        
        %*************** base filename
        penalty_ph2.basename = sprintf('%s_%s_%s%dtr_blk_%s',...
            ph1.basename, args.level, args.regress_type, ...
            args.shift_TRs, penalty_rest);
        penalty_ph3.basename = sprintf('%s_featsel_thresh%s', penalty_ph2.basename, num2str(args.featSelThresh));
        penalty_basename  = sprintf('classified_%s_%s', penalty_ph3.basename, args.classifier);
        
        %*************** load penalty
        loc_pen = load(sprintf('%s/penalty_check_%s.mat', dirs.mvpa.output{xph}, penalty_basename));%'pen_check'
        [xacc, whichmax] = max(loc_pen.pen_check.performance); %#ok<*NODEF>
        max_penalty      = loc_pen.pen_check.penalty(whichmax);
        args.penalty     = max_penalty;
        
        fprintf('... max_penalty: %s: acc: %1.4f\n', num2str(max_penalty), xacc);
        
        %*************** reset xpenalty
        args.xpenalty = args.penalty;
    end
    
    %*************** basename phase 4.
    ph4.basename  = sprintf('%s_penalty%s', class_basename, num2str(args.xpenalty));
    
    %*************** load phase 4.
    fname = sprintf('%s/%s.mat', dirs.mvpa.output{xph}, ph4.basename);
    load(fname);%'ph4'
    
    fprintf('... loaded classification results (penalty: %s) of %s: %s\n', ...
        num2str(args.xpenalty), args.subject_id, fname);
    
    %*************** save results in group structure
    mvpa_results{xsub}.results = ph4.results; %#ok<*NASGU>
    
end

%*************** save mvpa results
fname = sprintf('%s/mvpa_results_%s.mat', dirs.mvpa.group.mvpa_result{xph}, basename);
fprintf('(+) saving group mvpa_results: %s\n', fname);

save(fname, 'mvpa_results', '-v7.3')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%============= ROC for single subject + save mvpa results in group structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ============= SETUP
roc_matrix = []; perf_matrix = [];

%% ============= 1ST LEVEL
for xsub = xsub_grp
    fprintf('%s roc: %s\n', args.phase_name{xph}, subject_list(xsub).name);
    
    args.subject_id = subject_list(xsub).name;
    dirs = setup_directory(dirs, args);
    
    clear index xmatrix
    %*************** reset desired
    index     = create_design_index_localizer(args, dirs);
    n_runs    = index.param.n_runs;
    tt_matrix = index.matrix; xheader = index.header; n_header = length(xheader);
    
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
    %*************** unshifted spike array
    t_spike_array = tt_matrix(findCol(xheader, {'spike'}), :);
    
    xmatrix = t_matrix(:,(t_cate_array~=0) & (t_spike_array==0));
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
    
    for xrun=1:n_runs
        n_acts(xrun) = length(find(getDATA(xmatrix', xheader,{'run'}, ...
            {xrun},findCol(xheader, {'category'}))));
    end
    
    %*************** mvpaout
    xresults = mvpa_results{xsub}.results;
    
    n_target   = size(xresults.iterations(1).acts, 1);
    n_iter     = size(xresults.iterations, 2);
    for xiter = 1:n_iter
        n_vol(xiter) = size(xresults.iterations(xiter).acts, 2);
    end
    
    targ_array = 1:n_target;
    
    xacts = []; xdesireds = [];
    for xiter = 1:n_iter
        %*************** pfmet: based on shifted regressor / (startTR:end)
        xacts = horzcat(xacts, xresults.iterations(xiter).acts);
        xdesireds = horzcat(xdesireds, xresults.iterations(xiter).perfmet.desireds);
    end
    
    if ~isequal(xclass_array,xdesireds)
        fprintf('xdesired number is different: %s\n', args.subject_id);
    end
    
    n_target  = max(xdesireds);
    
    %% *************** ROC inputs: perfcurve(labels, scores, posclass)
    % scores: classifier predictions: (evidence)
    % labels: given true class labels: xdesireds(xguesses==xcate)
    % posclass: positive class label: (xcate)
    % T: (mean, lower bound, upper bound)
    
    for xcate = 1:n_target
        
        xlabels   = xdesireds;
        xscores   = xacts(xcate, :);
        xposclass = xcate;
        
        %*************** performance
        clear t_matrix
        
        t_matrix(:, 3:4) = [xlabels', xscores'];
        t_matrix(:, 1)   = xcate;
        t_matrix(:, 2)   = xsub;
        
        perf_matrix = vertcat(perf_matrix, t_matrix);
        
        %*************** perfcurve
        clear t_matrix
        
        [x_fa, y_hit, T, AUC] = perfcurve(xlabels, xscores, xposclass);
        
        t_matrix(:, 3:5) = [x_fa, y_hit, T];
        t_matrix(:, 1)   = xcate;
        t_matrix(:, 2)   = xsub;
        
        roc_matrix = vertcat(roc_matrix, t_matrix);
        
        mvpa_ROC{xph}.roc_out.rand_auc{xcate}(xsub, :) = [xsub, AUC];
    end
end

mvpa_ROC{xph}.perf.header    = {'category','subject','labels','scores'};
mvpa_ROC{xph}.perf.matrix    = perf_matrix;

mvpa_ROC{xph}.roc_out.header = {'category','subject','x_fa','y_hit','threshold'};
mvpa_ROC{xph}.roc_out.matrix = roc_matrix;

%*************** save mvpa_ROC
fname = sprintf('%s/mvpa_roc_%s.mat', xoutput_dir, basename);
save(fname, 'mvpa_ROC', '-v7.3')

%% ============= 1st LEVEL ROC PLOTTING
%*************** random effect
% mvpa_ROC{xph}.roc_out.rand_auc{xcate}(xsub,auc)

x_tick = 0:1;
x_lim  = [x_tick(1) x_tick(end)];
y_tick = x_tick; y_lim   = x_lim;
fig_w  = 1200; fig_h = fig_w;
lable_dist = 0.4;

xmatrix = mvpa_ROC{xph}.roc_out.matrix;
xheader = mvpa_ROC{xph}.roc_out.header;

for xcate = 1:n_target
    LogRegs_fig = figure;
    set(LogRegs_fig, 'Position', [0 0 fig_w fig_h])
    
    n_col = 10;
    if mod(length(args.g_sub),n_col) ~= 0, n_row = fix(length(args.g_sub)/n_col) + 1;
    else n_row = fix(length(args.n_sub)/n_col); end
    
    subplot(n_row, n_col, 1)
    plot([0 1], [0 1])
    xlabel('FA rate');
    ylabel('HIT rate');
    
    xh    = get(gca,'xlabel'); % handle to the label object
    xp    = get(xh,'position'); % get the current position property
    xp(2) = lable_dist * xp(2);% multiplt the distance, negative values put the label below the axis
    set(xh,'position', xp); % set the new position
    
    yh   = get(gca,'ylabel'); % handle to the label object
    yp    = get(yh,'position'); % get the current position property
    yp(1) = lable_dist * yp(1);% multiplt the distance, negative values put the label below the axis
    set(yh,'position', yp); % set the new position
    
    for xsub = args.filtered_subs
        xx_fa  = getDATA(xmatrix, xheader, {'category','subject'}, {xcate, xsub}, findCol(xheader,{'x_fa'}));
        yy_hit = getDATA(xmatrix, xheader, {'category','subject'}, {xcate, xsub}, findCol(xheader,{'y_hit'}));
        
        subplot(n_row, n_col, xsub)
        plot(xx_fa, yy_hit); hold on;
        plot([0 1], [0 1],'--k')
        
        set(gca,'xlim', x_lim, 'ylim', y_lim);
        set(gca,'XTick', x_tick, 'YTick', y_tick);
        
        if args.class_selecting
            title(sprintf('%s-%s (sub %s)', ...
                param.category_name{args.selected_category}, xcate_name{xcate}, num2str(sub_id(xsub))));
        else
            title(sprintf('%s (sub %s)', ...
                xcate_name{xcate}, num2str(xsub)));
        end
        
        xlabel('FA rate');
        ylabel('HIT rate');
        
        clear xh xp yh yp
        xh    = get(gca,'xlabel'); % handle to the label object
        xp    = get(xh,'position'); % get the current position property
        xp(2) = lable_dist * xp(2);% multiplt the distance, negative values put the label below the axis
        set(xh,'position', xp); % set the new position
        
        yh    = get(gca,'ylabel'); % handle to the label object
        yp    = get(yh,'position'); % get the current position property
        yp(1) = lable_dist * yp(1);% multiplt the distance, negative values put the label below the axis
        set(yh,'position', yp); % set the new position
    end
    hold off;
    
    %% save fig
    fig_fname = fullfile(xoutput_dir, sprintf('plot_mvpa_roc_%s_%s', basename, xcate_name{xcate}));
    
    savefig(LogRegs_fig, sprintf('%s.fig', fig_fname));
    set(gcf,'PaperUnits','inches','PaperPosition', [0 0 fig_w fig_h]/100)
    print('-djpeg', sprintf('%s.jpg',fig_fname), '-r100')
    saveas(LogRegs_fig, sprintf('%s.jpg',fig_fname), 'jpg')
    
    close(LogRegs_fig);
end

%% ============= 2nd LEVEL SUBJECT DATA: AUC / CLASSIFIER
%*************** random effect
% mvpa_ROC{xph}.roc_out.rand_auc{xcate}(xsub,auc)
n_condition = length(xcate_name);
n_subs      = length(xsub_groups);
xmeas_names = {'classifier_auc','classifier_accuracy','classifier_evidence'};

xout_txt    = fopen(sprintf('%s/classifier_stats_%s_n%s.txt', xoutput_dir, basename, num2str(n_subs)), 'w+');

fprintf(xout_txt, '####################################################################\n');
fprintf(xout_txt, '#### one-sample T test with normal distribution (mean: chance level)\n');
fprintf(xout_txt, '####################################################################\n\n');

for xmeas = 1:length(xmeas_names)%1_auc, 2_auc_baseline, 3_4_classifier accuracy, 4_classifier evidence
    clear xmatrix xnames xp xstats
    
    fprintf(xout_txt, '=========================================\n');
    fprintf(xout_txt, '======== %s\n\n', xmeas_names{xmeas});
    
    %*************** variable names
    xnames{1}    = 'subject';
    xnames{2}    = 'id';
    
    for xcate = 1:n_condition
        xnames{xcate + 2} = xcate_name{xcate};
    end
    
    %*************** data matrix
    xmatrix(:,1) = xsub_groups';
    xmatrix(:,2) = sub_id_filtered';
    
    if (xmeas == 1) || (xmeas == 2)
        for xcate = 1:n_condition
            if xmeas == 1, xdat = mvpa_ROC{xph}.roc_out.rand_auc{xcate}(:,2); end
            
            xmatrix = horzcat(xmatrix, xdat(xsub_groups));
        end
        
        xchance = 0.5;%for auc
        
    elseif (xmeas == 2) || (xmeas == 3)
        if xmeas == 2, xcsv = fullfile(dirs.mvpa.group.out{xph}, sprintf('grp_acc_matrix_%s_n%s.csv', basename, num2str(n_subs)));
        else,          xcsv = fullfile(dirs.mvpa.group.out{xph}, sprintf('grp_evi_matrix_%s_n%s.csv', basename, num2str(n_subs))); end
        
        xtable  = readtable(xcsv);
        tmatrix = table2array(xtable);
        xmatrix = horzcat(xmatrix, tmatrix);
        
        xchance = 1/n_condition;
    end
    
    %*************** total performance
    xmatrix = horzcat(xmatrix, mean(xmatrix(:,3:end),2));
    xnames{end + 1} = 'total';
    
    %*************** rest classification excluding rest in AUC
    if (xmeas==1) && (strcmp(args.rest, 'rest'))
        for xcate = 1:n_condition
            xmatrix = horzcat(xmatrix, mvpa_ROC{xph}.roc_out.rand_auc_ex_rest{xcate}(xsub_groups,2));
            xnames{xcate + 2 + n_condition} = sprintf('%s_ex_rest', xcate_name{xcate});
        end
    end
    
    %%============== one-sample ttest
    %*************** AUC baseline: 0.5
    for it_col = 3:(n_condition+3), [~, xp(it_col), ~, xstats{it_col}] = ttest(xmatrix(:,it_col), xchance); end
    
    for it_col = 3:(n_condition+3), fprintf(xout_txt, '%s\t', xnames{it_col}); end
    fprintf(xout_txt, '\n');
    for it_col = 3:(n_condition+3), fprintf(xout_txt, 'T(%s)=%4.4f ', num2str(xstats{it_col}.df), xstats{it_col}.tstat); end
    fprintf(xout_txt, '\n');
    for it_col = 3:(n_condition+3), fprintf(xout_txt, '%4.4f ', xp(it_col)); end
    fprintf(xout_txt, '\n');
    
    %*************** write tables to csv files
    xmatrix = vertcat(xmatrix, xp);
    xtable   = array2table(xmatrix, 'VariableNames', xnames);
    csv_name = sprintf('%s/%s_%s_n%s.csv', xoutput_dir, xmeas_names{xmeas}, basename, num2str(n_subs));
    writetable(xtable, csv_name,'WriteRowNames',true)
end
end