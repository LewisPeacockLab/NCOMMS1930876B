function[] = clearmem_localizer_classification_04(args, dirs)
% ROC/AUC analysis
% higher than matlab2014

%% ============= UNPACK ARGS.
xph               = args.xphase;
args.regress_type = args.train_regress_type;
param             = args.index{xph}.param;
subject_list      = args.subject_list;
xsub_grp          = args.filtered_subs;
n_iteration       = args.n_iteration;

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

if args.mvpa_out{xph}
    %% ============= LOAD EXISTING FILE
    %*************** load mvpa_results
    fname = sprintf('%s/mvpa_results_%s.mat', dirs.mvpa.group.mvpa_result{xph}, basename);
    if exist(fname, 'file'), load(fname); end % 'mvpa_results'
    
    %% ============= SETUP FILE NAMES
    %*************** ph1. base filename
    ph1.basename = sprintf('%s_%s_zscored_%s', args.phase_name{xph}, args.mask_name, args.epi_name);
    
    %*************** ph2. base filename
    if strcmp(args.regress_type, 'shift')
        ph2.basename = sprintf('%s_%s_%s%dtr_blk_%s',...
            ph1.basename, args.level, args.regress_type, ...
            args.shift_TRs, args.rest);
    elseif strcmp(args.regress_type, 'beta')
        ph2.basename = sprintf('%s_%s_%s', ...
            ph1.basename, args.level, args.regress_type);
    end
    
    %*************** reset ph2. filename
    if args.class_selecting
        ph2.basename = sprintf('cate%s_%s', sprintf('%d',args.selected_category), ph2.basename);
    end
    
    %*************** ph3. base filenames
    if args.featVox
        ph3.basename = sprintf('%s_featsel_%svox', ph2.basename, num2str(args.fsVosNum));
    else
        ph3.basename = sprintf('%s_featsel_thresh%s', ph2.basename, num2str(args.featSelThresh));
    end
    
    %*************** ph4. base filenames
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
            
            %*************** ph2. base filename
            if strcmp(args.regress_type, 'shift')
                penalty_ph2.basename = sprintf('%s_%s_%s%dtr_blk_%s',...
                    ph1.basename, args.level, args.regress_type, ...
                    args.shift_TRs, penalty_rest);
            elseif strcmp(args.regress_type, 'beta')
                penalty_ph2.basename = sprintf('%s_%s_%s', ...
                    ph1.basename, args.level, args.regress_type);
            end
            
            %*************** reset ph2. filename
            if args.class_selecting
                penalty_ph2.basename = sprintf('cate%s_%s', sprintf('%d',args.selected_category), penalty_ph2.basename);
            end
            
            %*************** ph3. base filename
            if args.featVox
                penalty_ph3.basename = sprintf('%s_featsel_%svox', penalty_ph2.basename, num2str(args.fsVosNum));
            else
                penalty_ph3.basename = sprintf('%s_featsel_thresh%s', penalty_ph2.basename, num2str(args.featSelThresh));
            end
            
            %*************** basename phase 4.
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
   
else
    %% *************** load mvpa results
    fname = sprintf('%s/mvpa_results_%s.mat', dirs.mvpa.group.mvpa_result{xph}, basename);
    fprintf('(+) loading group mvpa_results: %s\n', fname);
    
    load(fname);% 'mvpa_results'
    
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%============= ROC for single subject + save mvpa results in group structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if args.grp_auc{xph} 
    %*************** load mvpa_ROC
    fname = sprintf('%s/mvpa_roc_%s.mat', xoutput_dir, basename);
    if exist(fname, 'file'), load(fname); end % 'mvpa_ROC'     
    
    %% ============= SETUP 
    roc_matrix = []; perf_matrix = []; 
    roc_matrix_ex_rest = []; perf_matrix_ex_rest = []; 
    roc_matrix_each_category = []; perf_matrix_each_category = [];
    roc_matrix_each_category_ex_rest = []; perf_matrix_each_category_ex_rest = [];
    
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
        
        %*************** extract param array
%         for xrun = 1:n_runs
%             for xtrial = 1:n_trials(xrun)
%                 xcate = unique(getDATA(xmatrix', xheader, ...
%                     {'run','trial','spike'}, {xrun,xtrial,0}, ...
%                     findCol(xheader, {'category'})));
%                 xsubcate = unique(getDATA(xmatrix', xheader, ...
%                     {'run','trial','spike'}, {xrun,xtrial,0}, ...
%                     findCol(xheader, {'subcategory'})));
%                 
%                 if strcmp(args.level, 'category')
%                     xclass = xcate;
%                 elseif strcmp(args.level, 'subcategory')
%                     xclass = xsubcate + (n_subcategory * (xcate-1));
%                 end
%                 
%                 it_unit = xtrial + sum(n_trials(1:xrun-1));
%                 class_array(it_unit) = xclass;
%             end
%         end
        
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
%         
%         for xiter = 1:n_iter
%             %*************** pfmet: based on shifted regressor / (startTR:end)
%             if n_acts(xiter) ~= n_vol(xiter)
%                 xacts = horzcat(xacts, xresults.iterations(xiter).acts(:, 1:n_acts(xiter)));
%                 xdesireds = horzcat(xdesireds, xresults.iterations(xiter).perfmet.desireds(:, 1:n_acts(xiter)));
%             else
%                 xacts = horzcat(xacts, xresults.iterations(xiter).acts);
%                 xdesireds = horzcat(xdesireds, xresults.iterations(xiter).perfmet.desireds);
%             end
%         end
        
        if ~isequal(xclass_array,xdesireds)
            fprintf('xdesired number is different: %s\n', args.subject_id); 
        end
            
        n_target  = max(xdesireds);
        
        %*************** ROC inputs: perfcurve(labels, scores, posclass)
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
            
            %% ************* + rest-classification without rest feed-in
            if strcmp(args.rest, 'rest') && (xcate ~= n_target)
                xunit     = find(xdesireds~=n_target);%without rest
                xlabels   = xdesireds(xunit);
                xscores   = xacts(xcate, xunit);
                xposclass = xcate;
                
                %*************** performance
                clear t_matrix
                
                t_matrix(:, 3:4) = [xlabels', xscores'];
                t_matrix(:, 1)   = xcate;
                t_matrix(:, 2)   = xsub;
                
                perf_matrix_ex_rest = vertcat(perf_matrix_ex_rest, t_matrix);
                
                %*************** perfcurve
                clear t_matrix
                
                [x_fa, y_hit, T, AUC] = perfcurve(xlabels, xscores, xposclass);
                
                t_matrix(:, 3:5) = [x_fa, y_hit, T];
                t_matrix(:, 1)   = xcate;
                t_matrix(:, 2)   = xsub;
                
                roc_matrix_ex_rest = vertcat(roc_matrix_ex_rest, t_matrix);
                
                mvpa_ROC{xph}.roc_out.rand_auc_ex_rest{xcate}(xsub, :) = [xsub, AUC];
            end
            
            %% ************* separate category for all-subcategories classification
            if (strcmp(args.level, 'subcategory')) && ~(args.class_selecting) && (xcate <= 9)
                if strcmp(args.rest, 'rest')
                    targ_subcate = (1:n_subcategory) + (n_subcategory * fix((xcate - 1)/n_subcategory));
                    targ_subcate = [targ_subcate n_target];
                
                elseif strcmp(args.rest, 'norest') 
                    targ_subcate = (1:n_subcategory) + (n_subcategory * fix((xcate - 1)/n_subcategory));
                end
                
                xunit = [];
                for xtarg = targ_subcate 
                    xunit = horzcat(xunit, find(xdesireds==xtarg)); 
                end
                
                xunit     = sort(xunit);
                xlabels   = xdesireds(xunit);
                xscores   = xacts(xcate, xunit);
                xposclass = xcate;
                
                %*************** performance
                clear t_matrix
                
                t_matrix(:, 3:4) = [xlabels', xscores'];
                t_matrix(:, 1)   = xcate;
                t_matrix(:, 2)   = xsub;
                
                perf_matrix_each_category = vertcat(perf_matrix_each_category, t_matrix);
                
                %*************** perfcurve
                clear t_matrix
                
                [x_fa, y_hit, T, AUC] = perfcurve(xlabels, xscores, xposclass);
                
                t_matrix(:, 3:5) = [x_fa, y_hit, T];
                t_matrix(:, 1)   = xcate;
                t_matrix(:, 2)   = xsub;
                
                roc_matrix_each_category = vertcat(roc_matrix_each_category, t_matrix);
                
                mvpa_ROC{xph}.roc_out.rand_auc_each_category{xcate}(xsub, :) = [xsub, AUC];
                
                %*************** without rest from rest-classification
                if strcmp(args.rest, 'rest')
                    targ_subcate = (1:n_subcategory) + (n_subcategory * fix((xcate - 1)/n_subcategory));
                    
                    xunit = [];
                    for xtarg = targ_subcate
                        xunit = horzcat(xunit, find(xdesireds==xtarg));
                    end
                    
                    xlabels   = xdesireds(xunit);
                    xscores   = xacts(xcate, xunit);
                    xposclass = xcate;
                    
                    %*************** performance
                    clear t_matrix
                    
                    t_matrix(:, 3:4) = [xlabels', xscores'];
                    t_matrix(:, 1)   = xcate;
                    t_matrix(:, 2)   = xsub;
                    
                    perf_matrix_each_category_ex_rest = vertcat(perf_matrix_each_category_ex_rest, t_matrix);
                    
                    %*************** perfcurve
                    clear t_matrix
                    
                    [x_fa, y_hit, T, AUC] = perfcurve(xlabels, xscores, xposclass);
                    
                    t_matrix(:, 3:5) = [x_fa, y_hit, T];
                    t_matrix(:, 1)   = xcate;
                    t_matrix(:, 2)   = xsub;
                    
                    roc_matrix_each_category_ex_rest = vertcat(roc_matrix_each_category_ex_rest, t_matrix);
                    
                    mvpa_ROC{xph}.roc_out.rand_auc_each_category_ex_rest{xcate}(xsub, :) = [xsub, AUC];
                end
            end
        end
        
%         %%============= add permutation
%         bs_fname = sprintf('%s/auc_baseline_%s_%siters_%s.mat', ...
%             dirs.mvpa.group.auc_baseline{xph}, basename, num2str(args.n_iteration), args.subject_id);
%         xbase = load(bs_fname);%, 'baseline', '-v7.3')
%         mvpa_ROC{xph}.baseline{xsub} = xbase.baseline;
    end
    
    mvpa_ROC{xph}.perf.header    = {'category','subject','labels','scores'};
    mvpa_ROC{xph}.perf.matrix    = perf_matrix;
    
    mvpa_ROC{xph}.roc_out.header = {'category','subject','x_fa','y_hit','threshold'};
    mvpa_ROC{xph}.roc_out.matrix = roc_matrix;
    
    if strcmp(args.rest, 'rest')
        mvpa_ROC{xph}.roc_out.matrix_ex_rest = roc_matrix_ex_rest;
        mvpa_ROC{xph}.perf.matrix_ex_rest    = perf_matrix_ex_rest;
    end
    
    if ~(args.class_selecting) && (strcmp(args.level, 'subcategory'))
        mvpa_ROC{xph}.roc_out.matrix_each_category = roc_matrix_each_category;
        mvpa_ROC{xph}.perf.matrix_each_category    = perf_matrix_each_category;
        
        if strcmp(args.rest, 'rest')
            mvpa_ROC{xph}.roc_out.matrix_each_category_ex_rest = roc_matrix_each_category_ex_rest;
            mvpa_ROC{xph}.perf.matrix_each_category_ex_rest    = perf_matrix_each_category_ex_rest;    
        end
    end
    
    %*************** save mvpa_ROC
    fname = sprintf('%s/mvpa_roc_%s.mat', xoutput_dir, basename);
    save(fname, 'mvpa_ROC', '-v7.3')
    
else
    %% *************** load mvpa_ROC
    fname = sprintf('%s/mvpa_roc_%s.mat', xoutput_dir, basename);
    load(fname);% 'mvpa_ROC' 
    
end

%% ============= PLOTTING

if strcmp(args.cluster,'local')
    
    rmpath(genpath('/Applications/spm12'))
    
    %% ============= BASELINE HISTFIT PLOTTING
    %*************** baseline with auc score, criterion: upper 0.05
    % mvpa_ROC{xph}.baseline{xsub}{xclass}(xboots)
    
    x_tick = 0:0.2:1;
    x_lim  = [x_tick(1) x_tick(end)];
    y_tick = 0:10:80;
    y_lim  = [y_tick(1) y_tick(end)];
    fig_w  = 2400; fig_h = 1200;
    
    lable_dist = 0.4;
    n_col      = 10;
    
    n_bins     = 50;
    hist_color = [0.6 0.6 0.6];
    line_color = 'r';
    xalpha     = 0.05;
    
    for xcond = 1:n_target
        
        xfig_hist = figure;
        set(xfig_hist, 'Position', [0 0 fig_w fig_h])
        
        %***************** plotting
        for xsub = 1:args.n_sub
            if mod(args.n_sub,n_col) ~= 0
                n_row = fix(args.n_sub/n_col) + 1;
            else
                n_row = fix(args.n_sub/n_col);
            end
            
            subplot(n_row, n_col, xsub)
            plot([0 0], [0 0]); hold on
            
            set(gca,'xlim', x_lim, 'ylim', y_lim);
            set(gca,'XTick', x_tick, 'YTick', y_tick);
            
            title(sprintf('%s (s%s)', xcate_name{xcond}, num2str(sub_id(xsub))));
            xlabel('AUC');
            ylabel('frequency');
        end
        
        for xsub = args.filtered_subs
            
            xrand_aucs = mvpa_ROC{xph}.baseline{xsub}{xcond};
            
            if mod(args.n_sub,n_col) ~= 0
                n_row = fix(args.n_sub/n_col) + 1;
            else
                n_row = fix(args.n_sub/n_col);
            end
            
            subplot(n_row, n_col, xsub)
            h = histfit(xrand_aucs, n_bins, 'normal');
            h(1).FaceColor = hist_color;
            h(2).Color     = line_color;
            
            hold on
            
            set(gca,'xlim', x_lim, 'ylim', y_lim);
            set(gca,'XTick', x_tick, 'YTick', y_tick);
            
        end
        hold off;
        
        %% save fig
        fig_fname = fullfile(xoutput_dir, sprintf('plot_histfit_roc_baseline_%s_%s', basename, xcate_name{xcond}));
        
        %     savefig(xfig_hist, sprintf('%s.fig', fig_fname));
        set(gcf,'PaperUnits','inches','PaperPosition', [0 0 fig_w fig_h]/100)
        print('-djpeg', sprintf('%s.jpg',fig_fname), '-r100')
        saveas(xfig_hist, sprintf('%s.jpg',fig_fname), 'jpg')
        
%         close(xfig_hist);
        
    end
    
    %% ============= BASELINE NORMFIT PLOTTING
    %*************** baseline with auc score, criterion: upper 0.05
    % mvpa_ROC{xph}.baseline{xsub}{xcond}(xboots)
    
    x_tick = 0:0.2:1;
    x_lim  = [x_tick(1) x_tick(end)];
    if strcmp(args.level, 'category')
        y_tick = 0:2:28;
    else
        y_tick = 0:2:16;
    end
    y_lim  = [y_tick(1) y_tick(end)];
    fig_w  = 2400; fig_h = 1200;
    
    lable_dist = 0.4;
    n_col      = 10;
    
    n_bins     = 50;
    hist_color = [0.6 0.6 0.6];
    line_color = 'r';
    xalpha     = 0.05;
    
    for xcond = 1:n_target
        
        xfig_norm = figure;
        set(xfig_norm, 'Position', [0 0 fig_w fig_h])
        
        %***************** plotting
        for xsub = 1:args.n_sub
            if mod(args.n_sub,n_col) ~= 0
                n_row = fix(args.n_sub/n_col) + 1;
            else
                n_row = fix(args.n_sub/n_col);
            end
            
            subplot(n_row, n_col, xsub)
            plot([0 0], [0 0]); hold on
            
            set(gca,'xlim', x_lim, 'ylim', y_lim);
            set(gca,'XTick', x_tick, 'YTick', y_tick);
            
            title(sprintf('%s (s %s)', xcate_name{xcond}, num2str(sub_id(xsub))));
            xlabel('AUC');
            ylabel('PDF (probability density)');
        end
        
        for xsub = args.filtered_subs
            %***************** target auc
            xauc = mvpa_ROC{xph}.roc_out.rand_auc{xcond}(xsub,2);
            
            %***************** random aucs
            xrand_aucs = mvpa_ROC{xph}.baseline{xsub}{xcond};
            
            xmu = mean(xrand_aucs); xsigma = std(xrand_aucs);
            xcriterion = norminv((1-xalpha), xmu, xsigma);
            
            mvpa_ROC{xph}.baseline_criterion{xcond}(xsub,1) = xcriterion;
            
            if mod(args.n_sub,n_col) ~= 0
                n_row = fix(args.n_sub/n_col) + 1;
            else
                n_row = fix(args.n_sub/n_col);
            end
            
            subplot(n_row, n_col, xsub)
            %[h, p] = normspec([Inf, xcriterion], mean(xrand_aucs), std(xrand_aucs), 'outside');
            prob = (0.0002:0.0004:0.9998)';
            xx   = norminv(prob,xmu,xsigma);
            yy   = normpdf(xx,xmu,xsigma);
            
            plot(xx, yy, '-k'); hold on
            
            xbound = xx(xx >= xcriterion);
            ybound = yy(xx >= xcriterion);
            xfill = [xbound(1); xbound(1); xbound; xbound(end); xbound(end)];
            yfill = [        0; ybound(1); ybound; ybound(end); 0];
            fill(xfill,yfill, hist_color);
            
            plot([xcriterion xcriterion], [y_lim(1) y_lim(2)], '--k')
            plot([xauc xauc], [y_lim(1) y_lim(2)], '--r')
            
            set(gca,'xlim', x_lim, 'ylim', y_lim);
            set(gca,'XTick', x_tick, 'YTick', y_tick);
            
        end
        hold off;
        
        %% save fig
        fig_fname = fullfile(xoutput_dir, sprintf('plot_norm_roc_baseline_%s_%s', basename, xcate_name{xcond}));
        
        %     savefig(xfig_norm, sprintf('%s.fig', fig_fname));
        set(gcf,'PaperUnits','inches','PaperPosition', [0 0 fig_w fig_h]/100)
        print('-djpeg', sprintf('%s.jpg',fig_fname), '-r100')
        saveas(xfig_norm, sprintf('%s.jpg',fig_fname), 'jpg')
        
        close(xfig_norm);
        
    end
    
    fname = sprintf('%s/mvpa_roc_%s.mat', xoutput_dir, basename);
    save(fname, 'mvpa_ROC', '-v7.3')
    
    %% ============= BASELINE BAR PLOTTING
    %*************** baseline with auc score, criterion: upper 0.05
    % mvpa_ROC{xph}.baseline{xsub}{xcond}(xboots)
    for i = 2%1_baseline, 2_0.5 cutoff
        auc_cutoff_subs = [];
        
        x_tick = 0:(n_class+1);
        x_lim  = [x_tick(1) x_tick(end)];
        y_tick = 0:0.2:1;
        y_lim  = [y_tick(1) y_tick(end)];
        fig_w  = 2400; fig_h = 1200;
        
        lable_dist = 0.4;
        n_col      = 10;
        
        n_bins     = 50;
        hist_color = [0.6 0.6 0.6];
        line_color = 'r';
        xalpha     = 0.05;
        
        xfig_norm = figure;
        set(xfig_norm, 'Position', [0 0 fig_w fig_h])
        
        %***************** plotting
        for xsub = 1:args.n_sub
            if mod(args.n_sub,n_col) ~= 0
                n_row = fix(args.n_sub/n_col) + 1;
            else
                n_row = fix(args.n_sub/n_col);
            end
            
            subplot(n_row, n_col, xsub)
            plot([0 0], [0 0]); hold on
            
            set(gca,'xlim', x_lim, 'ylim', y_lim);
            set(gca,'XTick', x_tick, 'YTick', y_tick);
            
            title(sprintf('subject: %s', num2str(sub_id(xsub))));
            xlabel('operation');
            ylabel('AUC');
        end
        
        for xsub = args.filtered_subs
            
            %***************** target auc
            for xcond = 1:n_target
                aucs(xcond)       = mvpa_ROC{xph}.roc_out.rand_auc{xcond}(xsub,2);
                if i==1
                    criterions(xcond) = mvpa_ROC{xph}.baseline_criterion{xcond}(xsub,1);
                else
                    criterions(xcond) = 0.5;
                end
            end
            
            auc_cutoff{xsub} = find((aucs - criterions) < 0);
            
            if ~isempty(auc_cutoff{xsub})
                xcutoff = zeros(1,n_target);
                xcutoff(auc_cutoff{xsub}) = criterions(auc_cutoff{xsub});
                
                auc_cutoff_subs = horzcat(auc_cutoff_subs, xsub);
            end
            
            if mod(args.n_sub,n_col) ~= 0
                n_row = fix(args.n_sub/n_col) + 1;
            else
                n_row = fix(args.n_sub/n_col);
            end
            
            subplot(n_row, n_col, xsub)
            
            bar(aucs, 0.8, 'FaceColor', xcond_color{3}); hold on
            if ~isempty(auc_cutoff{xsub})
                bar(xcutoff, 0.8, 'FaceColor', hist_color);
            end
            
            for xcond = 1:n_target
                p = plot([0.6 1.4]+(xcond-1), [criterions(xcond) criterions(xcond)], '-b');
                p(1).LineWidth = 1;
            end
            
            set(gca,'xlim', x_lim, 'ylim', y_lim);
            set(gca,'XTick', x_tick, 'YTick', y_tick);
            
        end
        hold off;
        
        %*************** save fig
        if i==1
            fig_fname = fullfile(xoutput_dir, sprintf('plot_cutoff_roc_baseline_%s', basename));
        else
            fig_fname = fullfile(xoutput_dir, sprintf('plot_%1.1fcutoff_roc_baseline_%s', criterions(1), basename));
        end
        
        % savefig(xfig_norm, sprintf('%s.fig', fig_fname));
        set(gcf,'PaperUnits','inches','PaperPosition', [0 0 fig_w fig_h]/100)
        print('-djpeg', sprintf('%s.jpg',fig_fname), '-r100')
        saveas(xfig_norm, sprintf('%s.jpg',fig_fname), 'jpg')
        
        close(xfig_norm);
    end
    
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
    
    %% ============= violin plot for AUC
    xauc_data = []; xauc_baseline = []; xcolors = [];
    for xcond = 1:n_target
        xauc_data(:,xcond)     = mvpa_ROC{xph}.roc_out.rand_auc{xcond}(xsub_groups,2);
%         xauc_baseline(:,xcond) = mvpa_ROC{xph}.baseline_criterion{xcond}(xsub_groups,1);
        
        if strcmp(args.level,'category'), xclr = xcond;
        else, xclr = fix((xcond-1)/n_subcategory) + 1; end
        
        xcolors(xcond, :) = xcate_color{xclr};
    end
    
    %*************** figure
    if strcmp(args.level, 'category')
        fig_rect = [0 0 500 1000];
        y_lim    = [0.5 1.1];
    elseif strcmp(args.level, 'subcategory')
        fig_rect = [0 0 1000 1000];
        y_lim    = [0.2 1.1];
    end
    
    violin_fig = figure;
    set(violin_fig, 'Position', fig_rect);
    
    y_tick = y_lim(1):0.1:y_lim(2);
    x_tick = 1:n_target;
    
    %*************** legend
    h    = zeros(n_target, 1);
    for i=1:n_target
        h(i) = plot(NaN, NaN, 'Color', xcolors(i, :),'LineWidth',5);
        hold on;
    end

    xlegend       = xcate_name;
    if strcmp(args.level, 'category')
        lg = legend(xlegend{1}, xlegend{2}, xlegend{3},'AutoUpdate','off');
    else
        lg = legend(xlegend{1}, xlegend{2}, xlegend{3},...
            xlegend{4}, xlegend{5}, xlegend{6},...
            xlegend{7}, xlegend{8}, xlegend{9},...
            'AutoUpdate','off');
    end
    lg.Location   = 'BestOutside';
    
    %*************** violin
    violin(xauc_data,'facecolor',xcolors,'xlabel', xcate_name,...
        'edgecolor','','medc','','mc','k','plotlegend',''); hold on
    % violin(xauc_baseline,'facecolor',xbase_color,'xlabel', xcate_name,...
    %     'edgecolor','','medc','','mc','k','plotlegend','');
    
    plot([x_tick(1)-1 x_tick(end)+1], [0.5 0.5], ':k', 'LineWidth',2);
    
    title(sprintf('Localizer ACU: %s-level (N=%s)', args.level, num2str(length(xsub_groups))));
    
    set(gca,'ylim', y_lim)
    set(gca,'YTick', y_tick)
    set(gca,'XTick', x_tick, 'XTickLabel', xcate_name)
    xlabel('target category');
    ylabel('AUC');
    grid on
    
    %*************** save fig
    fig_fname = fullfile(xoutput_dir, sprintf('violin_AUC_%s_n%s', basename, num2str(length(xsub_groups))));
    
    savefig(violin_fig, sprintf('%s.fig', fig_fname));
    set(gcf,'PaperUnits','inches','PaperPosition', fig_rect/100)
    print('-djpeg', sprintf('%s.jpg',fig_fname), '-r100')
    saveas(violin_fig, sprintf('%s.jpg',fig_fname), 'jpg')
    
    close(violin_fig);
    
    %% ============= BOXPLOT for AUC
    xauc_data = []; xauc_baseline = []; xcolors = [];
    for xcond = 1:n_target
        xauc_data(:,xcond)     = mvpa_ROC{xph}.roc_out.rand_auc{xcond}(xsub_groups,2);
        xauc_baseline(:,xcond) = 0.5;%mvpa_ROC{xph}.baseline_criterion{xcond}(xsub_groups,1);
        
        if strcmp(args.level,'category'), xclr = xcond;
        else, xclr = fix((xcond-1)/n_subcategory) + 1; end
        
        xcolors(xcond, :) = xcate_color{xclr};
    end
    
    %*************** figure
    if strcmp(args.level, 'category')
        fig_rect = [0 0 500 1000];
        y_lim    = [0.5 1.0];
    elseif strcmp(args.level, 'subcategory')
        fig_rect = [0 0 1000 1000];
        y_lim    = [0.3 1.0];
    end
    
    box_fig = figure;
    set(box_fig, 'Position', fig_rect)
    
    y_tick = y_lim(1):0.1:y_lim(2);
    x_tick = 1:n_target;
    
    %*************** legend
    h    = zeros(n_target+1, 1);
    for i=1:n_target
        h(i) = plot(NaN, NaN, 'Color', xcolors(i, :),'LineWidth',5);
        hold on;
    end
    h(n_target+1) = plot(NaN, NaN, 'Color', xbase_color,'LineWidth',5);
    xlegend       = xcate_name;
    if strcmp(args.level, 'category')
        lg = legend(xlegend{1}, xlegend{2}, xlegend{3},'AutoUpdate','off');
    else
        lg = legend(xlegend{1}, xlegend{2}, xlegend{3},...
            xlegend{4}, xlegend{5}, xlegend{6},...
            xlegend{7}, xlegend{8}, xlegend{9},...
            'AutoUpdate','off');
    end
    lg.Location   = 'BestOutside';
    
    %*************** box plot
    boxplot(xauc_data, 'Notch','on','Colors',xcolors,'BoxStyle','outline','Widths',0.5,'Symbol','o')
    
    plot([x_tick(1)-1 x_tick(end)+1], [0.5 0.5], ':k', 'LineWidth',2);
    
    title(sprintf('Localizer ACU: %s-level (N=%s)', args.level, num2str(length(xsub_groups))));
    
    set(gca,'ylim', y_lim)
    set(gca,'YTick', y_tick)
    set(gca,'XTick', x_tick, 'XTickLabel', xcate_name)
    xlabel('target category');
    ylabel('AUC');
    grid on
    
    %*************** save fig
    fig_fname = fullfile(xoutput_dir, sprintf('boxplot_AUC_%s_n%s', basename, num2str(length(xsub_groups))));
    
    savefig(box_fig, sprintf('%s.fig', fig_fname));
    set(gcf,'PaperUnits','inches','PaperPosition', fig_rect/100)
    print('-djpeg', sprintf('%s.jpg',fig_fname), '-r100')
    saveas(box_fig, sprintf('%s.jpg',fig_fname), 'jpg')
    
    close(box_fig);
    
    %% ============= 2nd LEVEL SUBJECT DATA: AUC / CLASSIFIER
    %*************** random effect
    % mvpa_ROC{xph}.roc_out.rand_auc{xcate}(xsub,auc)
    n_condition = length(xcate_name);
    n_subs      = length(xsub_groups);
    xmeas_names = {'classifier_auc','classifier_accuracy','classifier_evidence'};
%     xmeas_names = {'classifier_auc','classifier_auc_baseline','classifier_accuracy','classifier_evidence'};
    
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
        
        if (xmeas == 1)% || (xmeas == 2)
            for xcate = 1:n_condition
                if xmeas == 1, xdat = mvpa_ROC{xph}.roc_out.rand_auc{xcate}(:,2); end
%                 else,          xdat = mvpa_ROC{xph}.baseline_criterion{xcate}(:,1); end
                
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
    
    %% ============= 2nd LEVEL SUBJECT DATA: EACH CATEGORY FROM ALL SUBCATEGORY CLASSIFIER
    %*************** random effect
    if (strcmp(args.level, 'subcategory')) && ~(args.class_selecting)
        clear xmatrix xnames
        
        xmatrix = args.filtered_subs';
        xnames{1} = 'subject';
        
        for xcate = 1:(n_target-1)
            xmatrix = horzcat(xmatrix, mvpa_ROC{xph}.roc_out.rand_auc_each_category{xcate}(xsub_groups,2));
            xnames{xcate + 1} = xcate_name{xcate};
        end
        
        %*************** rest classification excluding rest in AUC
        if strcmp(args.rest, 'rest')
            for xcate = 1:(n_target-1)
                xmatrix = horzcat(xmatrix, mvpa_ROC{xph}.roc_out.rand_auc_each_category_ex_rest{xcate}(:,2));
                xnames{xcate + 1 + (n_target-1)} = sprintf('%s_ex_rest', xcate_name{xcate});
            end
        end
        
        table_auc = array2table(xmatrix, 'VariableNames', xnames);
        
        %*************** write tables to csv files
        csv_name = sprintf('%s/classifier_auc_%s_each_category.csv', xoutput_dir, basename);
        writetable(table_auc, csv_name,'WriteRowNames',true)
        
    end
end
end