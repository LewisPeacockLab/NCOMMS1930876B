function[] = clearmem_mvpa_operation_05(args, dirs)
% ROC/AUC analysis
% higher than matlab2014

%% ============= UNPACK ARGS.
xph               = args.xphase;
args.regress_type = args.test_regress_type;
subject_list      = args.subject_list;

n_iteration       = args.n_iteration;

%*************** subject id num
for xsub = 1:args.n_sub
    sub_id(xsub) = str2double(subject_list(xsub).name(end-2:end));
end

%*************** filtered subject id num
for it = 1:length(args.filtered_subs)
    xsub = args.filtered_subs(it);
    sub_id_filtered(it) = str2double(subject_list(xsub).name(end-2:end));
end

%*************** output basename
basename          = args.analysis_basename;
xoutput_dir       = dirs.mvpa.group.auc{xph};

if args.four_oper_regress
    it_conds = [1 2 4 5];
    cond_names = {'Maintain','RepCategory','Suppress','Clear'};
else
    it_conds = 1:5;
    cond_names = {'Maintain','RepCategory','RepSubcate','Suppress','Clear'};
end

n_target   = length(cond_names);

%*************** plot param
for it = 1:length(it_conds)
    xcond           = it_conds(it);
    xcond_color{it} = args.cond_color{xcond};
end
xbase_color = args.base_color;
font_size   = 10;

if args.fixed_penalty{xph}
    xsub_groups = args.filtered_subs;
else
    xsub_groups = args.g_sub;
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%============= LOAD CLASSIFICATION ACCURACY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fprintf('(+) load 2nd level mvpaout\n');
% 
% %%*************** 2ND LEVEL LOAD
% fname = fullfile(dirs.mvpa.group.out{xph}, sprintf('grp_mvpaout_%s.mat', basename));
% load(fname);%'grp_mvpaout'

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%============= ROC for single subject + save mvpa results in group structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if args.grp_auc{3}
    %% ============= load group mvpa results from clearmem_mvpa_operation_04
    fname = sprintf('%s/mvpa_results_%s.mat', dirs.mvpa.group.mvpa_result{xph}, basename);
    load(fname);%, 'mvpa_results'
    
    %% ============= SETUP         
    roc_matrix = []; perf_matrix = []; 
    roc_matrix_ex_rest = []; perf_matrix_ex_rest = []; 
    
    %% ============= 1ST LEVEL
    for xsub = xsub_groups
        fprintf('%s roc: %s\n', args.phase_name{xph}, subject_list(xsub).name);
        
        %*************** reset desired
        args.subject_id = subject_list(xsub).name;
        dirs = setup_directory(dirs, args);
        
        clear index xmatrix
        load(fullfile(dirs.param, sprintf('%s_parameters_%s.mat', args.phase_name{2}, args.subject_id))); %index
        xmatrix  = index.matrix; xheader = index.header;
        xresults = mvpa_results{xsub}.results; %#ok<USENS>
        
        n_target = size(xresults.iterations(1).acts, 1);
        n_iter   = size(xresults.iterations, 2);
        for xiter = 1:n_iter
            n_vol(xiter) = size(xresults.iterations(xiter).acts, 2);
        end
        
        %***************** reset regs: shifting TR
        xregs.runs         = args.regs{xph}.selectors;
        xregs.operation_sh = zeros(1, sum(n_vol));% only operation window
        
        for it = 1:n_target
            xcond = it_conds(it);
            xunit = getDATA(xmatrix', xheader, {'condition'}, {xcond});
            
            xregs.timecourse(xunit) = xcond;
            
            xunit = find(getDATA(xmatrix', xheader, ...
                {'condition','manipulation'}, {xcond,1}));
            
            xregs.operation_sh(xunit + args.shift_TRs) = it;
        end
        
        %***************** remove spikes
        xspike = getDATA(xmatrix', xheader, {'spike'}, {1});
        xregs.operation_s(xspike) = 0;
        
        %*************** extract MVPA performance
        acts_array = [];
        for xiter = 1:n_iter
            acts_array = horzcat(acts_array, xresults.iterations(xiter).acts); %#ok<*AGROW>
        end
        
        xunit     = (xregs.operation_sh ~=0);
        xdesireds = xregs.operation_sh(xunit);
        xacts     = acts_array(:,xunit);
        
        %*************** ROC inputs: perfcurve(labels, scores, posclass)
        % scores: classifier predictions: (evidence)
        % labels: given true class labels: xdesireds
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
            
            %%============= no permutation
            mvpa_ROC{xph}.baseline{xsub}{xcate} = 1/n_target;
            
        end
        
        %%============= add permutation
%         if args.four_oper_regress
%             bs_fname = sprintf('%s/auc_baseline_four_oper_%siters_%s.mat', ...
%                 dirs.mvpa.group.auc_baseline{xph}, num2str(n_iteration), args.subject_id); %#ok<*UNRCH>
%         else
%             bs_fname = sprintf('%s/auc_baseline_%siters_%s.mat', ...
%                 dirs.mvpa.group.auc_baseline{xph}, num2str(n_iteration), args.subject_id);
%         end
%         
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
    
    %*************** save mvpa_ROC
    fname = sprintf('%s/mvpa_roc_%s.mat', xoutput_dir, basename);
    save(fname, 'mvpa_ROC', '-v7.3')
    
else
    %% *************** load mvpa_ROC
    fname = sprintf('%s/mvpa_roc_%s.mat', xoutput_dir, basename);
    load(fname);% 'mvpa_ROC'   
    
    for xsub = xsub_groups
        mvpa_ROC{xph}.baseline{xsub} = [];
        
        for xcond = 1:n_target    
            mvpa_ROC{xph}.baseline{xsub}{xcond} = 1/n_target;
        end
    end
    
end

%% ============= PLOTTING

if strcmp(args.cluster, 'local')
    %% ============= BASELINE HISTFIT PLOTTING
    %*************** baseline with auc score, criterion: upper 0.05
    % mvpa_ROC{xph}.baseline{xsub}{xcond}(xboots)
    
%     x_tick = 0:0.2:1;
%     x_lim  = [x_tick(1) x_tick(end)];
%     y_tick = 0:10:80;
%     y_lim  = [y_tick(1) y_tick(end)];
%     fig_w  = 2400; fig_h = 1200;
%     
%     lable_dist = 0.4;
%     n_col      = 10;
%     
%     n_bins     = 50;
%     hist_color = [0.6 0.6 0.6];
%     line_color = 'r';
%     xalpha     = 0.05;
%     
%     for xcond = 1:n_target
%         
%         xfig_hist = figure;
%         set(xfig_hist, 'Position', [0 0 fig_w fig_h])
%         
%         %***************** plotting
%         for xsub = 1:args.n_sub
%             if mod(args.n_sub,n_col) ~= 0
%                 n_row = fix(args.n_sub/n_col) + 1;
%             else
%                 n_row = fix(args.n_sub/n_col);
%             end
%             
%             subplot(n_row, n_col, xsub)
%             plot([0 0], [0 0]); hold on
%             
%             set(gca,'xlim', x_lim, 'ylim', y_lim);
%             set(gca,'XTick', x_tick, 'YTick', y_tick);
%             
%             title(sprintf('%s (s%s)', cond_names{xcond}, num2str(sub_id(xsub))));
%             xlabel('AUC');
%             ylabel('frequency');
%         end
%         
%         for xsub = args.g_sub
%             clear h
%             
%             xrand_aucs = mvpa_ROC{xph}.baseline{xsub}{xcond};
%             
%             if mod(args.n_sub,n_col) ~= 0
%                 n_row = fix(args.n_sub/n_col) + 1;
%             else
%                 n_row = fix(args.n_sub/n_col);
%             end
%             
%             subplot(n_row, n_col, xsub)
%             h = histfit(xrand_aucs, n_bins, 'normal');
%             h(1).FaceColor = hist_color;
%             h(2).Color     = line_color;
%             
%             hold on
%             
%             set(gca,'xlim', x_lim, 'ylim', y_lim);
%             set(gca,'XTick', x_tick, 'YTick', y_tick);
%             
%         end
%         hold off;
%         
%         %% save fig
% %         fig_fname = fullfile(xoutput_dir, sprintf('plot_histfit_roc_baseline_%s_%s', basename, cond_names{xcond}));
%         
% %         savefig(xfig_hist, sprintf('%s.fig', fig_fname));
% %         set(gcf,'PaperUnits','inches','PaperPosition', [0 0 fig_w fig_h]/100)
% %         print('-djpeg', sprintf('%s.jpg',fig_fname), '-r100')
% %         saveas(xfig_hist, sprintf('%s.jpg',fig_fname), 'jpg')
%         
%         close(xfig_hist);
%         
%     end
    
    %% ============= BASELINE NORMFIT PLOTTING
    %*************** baseline with auc score, criterion: upper 0.05
    % mvpa_ROC{xph}.baseline{xsub}{xcond}(xboots)
    
    x_tick = 0:0.2:1;
    x_lim  = [x_tick(1) x_tick(end)];
    y_tick = 0:2:20;
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
            
            title(sprintf('%s (s %s)', cond_names{xcond}, num2str(sub_id(xsub))));
            xlabel('AUC');
            ylabel('PDF (probability density)');
        end
        
        for xsub = args.g_sub
            %***************** target auc
            xauc = mvpa_ROC{xph}.roc_out.rand_auc{xcond}(xsub,2);
            
            %***************** random aucs
            xrand_aucs = sort(mvpa_ROC{xph}.baseline{xsub}{xcond});
            
            %***************** normal fitting
            xmu = mean(xrand_aucs); xsigma = std(xrand_aucs);
            xcriterion = norminv((1-xalpha), xmu, xsigma);
            
            mvpa_ROC{xph}.baseline_criterion{xcond}(xsub,1) = xcriterion;
            
            %***************** order 95%
            xcriterion_order = xrand_aucs(1000*(1-xalpha));
            mvpa_ROC{xph}.baseline_criterion_order{xcond}(xsub,1) = xcriterion_order;
            
            %***************** plotting
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
            plot([xcriterion_order xcriterion_order], [y_lim(1) y_lim(2)], '--b')
            
            plot([xauc xauc], [y_lim(1) y_lim(2)], '--r')
            
            set(gca,'xlim', x_lim, 'ylim', y_lim);
            set(gca,'XTick', x_tick, 'YTick', y_tick);
            
        end
        hold off;
        
        %% save fig
%         fig_fname = fullfile(xoutput_dir, sprintf('plot_norm_roc_baseline_%s_%s', basename, cond_names{xcond}));
        
%         savefig(xfig_norm, sprintf('%s.fig', fig_fname));
%         set(gcf,'PaperUnits','inches','PaperPosition', [0 0 fig_w fig_h]/100)
%         print('-djpeg', sprintf('%s.jpg',fig_fname), '-r100')
%         saveas(xfig_norm, sprintf('%s.jpg',fig_fname), 'jpg')
%         
        close(xfig_norm);
        
    end
    
    fname = sprintf('%s/mvpa_roc_%s.mat', xoutput_dir, basename);
    save(fname, 'mvpa_ROC', '-v7.3')
    
    %% ============= BASELINE BAR PLOTTING
    %*************** baseline with auc score, criterion: upper 0.05
    % mvpa_ROC{xph}.baseline{xsub}{xcond}(xboots)
%     x_tick = 0:(n_target+1);
%     x_lim  = [x_tick(1) x_tick(end)];
%     y_tick = 0:0.2:1;
%     y_lim  = [y_tick(1) y_tick(end)];
%     fig_w  = 2400; fig_h = 1200;
%     
%     lable_dist = 0.4;
%     n_col      = 10;
%     
%     n_bins     = 50;
%     hist_color = [0.6 0.6 0.6];
%     line_color = 'r';
%     xalpha     = 0.05;
%     
%     for i = 1:2
%         
%         auc_cutoff_subs = [];
%         
%         xfig_norm = figure;
%         set(xfig_norm, 'Position', [0 0 fig_w fig_h])
%         
%         %***************** plotting
%         for xsub = 1:args.n_sub
%             if mod(args.n_sub,n_col) ~= 0
%                 n_row = fix(args.n_sub/n_col) + 1;
%             else
%                 n_row = fix(args.n_sub/n_col);
%             end
%             
%             subplot(n_row, n_col, xsub)
%             plot([0 0], [0 0]); hold on
%             
%             set(gca,'xlim', x_lim, 'ylim', y_lim);
%             set(gca,'XTick', x_tick, 'YTick', y_tick);
%             
%             title(sprintf('subject: %s', num2str(sub_id(xsub))));
%             xlabel('operation');
%             ylabel('AUC');
%         end
%         
%         for xsub = args.g_sub
%             
%             %***************** target auc
%             for xcond = 1:n_target
%                 aucs(xcond) = mvpa_ROC{xph}.roc_out.rand_auc{xcond}(xsub,2);
%                 if i==1
%                     criterions(xcond) = mvpa_ROC{xph}.baseline_criterion{xcond}(xsub,1);
%                 else
%                     criterions(xcond) = mvpa_ROC{xph}.baseline_criterion_order{xcond}(xsub,1);
%                 end
%             end
%             
%             auc_cutoff{xsub} = find((aucs - criterions) < 0);
%             
%             if ~isempty(auc_cutoff{xsub})
%                 xcutoff = zeros(1,n_target);
%                 xcutoff(auc_cutoff{xsub}) = criterions(auc_cutoff{xsub});
%                 
%                 auc_cutoff_subs = horzcat(auc_cutoff_subs, xsub);
%             end
%             
%             if mod(args.n_sub,n_col) ~= 0
%                 n_row = fix(args.n_sub/n_col) + 1;
%             else
%                 n_row = fix(args.n_sub/n_col);
%             end
%             
%             subplot(n_row, n_col, xsub)
%             
%             bar(aucs, 0.8, 'FaceColor', xcond_color{3}); hold on
%             if ~isempty(auc_cutoff{xsub})
%                 bar(xcutoff, 0.8, 'FaceColor', hist_color);
%             end
%             
%             for xcond = 1:n_target
%                 p = plot([0.6 1.4]+(xcond-1), [criterions(xcond) criterions(xcond)], '-b');
%                 p(1).LineWidth = 1;
%             end
%             
%             set(gca,'xlim', x_lim, 'ylim', y_lim);
%             set(gca,'XTick', x_tick, 'YTick', y_tick);
%             
%         end
%         hold off;
%         
%         %*************** save fig
%         if i==1
%             fig_fname = fullfile(xoutput_dir, sprintf('plot_cutoff_roc_baseline_%s', basename));
%         else
%             fig_fname = fullfile(xoutput_dir, sprintf('plot_cutoff_roc_order_baseline_%s', basename));
%         end
% %         savefig(xfig_norm, sprintf('%s.fig', fig_fname));
%         set(gcf,'PaperUnits','inches','PaperPosition', [0 0 fig_w fig_h]/100)
%         print('-djpeg', sprintf('%s.jpg',fig_fname), '-r100')
%         saveas(xfig_norm, sprintf('%s.jpg',fig_fname), 'jpg')
%         
% %         close(xfig_norm);
%     end
    
    %% ============= BASELINE BAR PLOTTING
    %*************** baseline with auc score, criterion: upper 0.05
    % mvpa_ROC{xph}.baseline{xsub}{xcond}(xboots)
    x_tick = 0:6;
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
    
    for i = 1
        
        auc_cutoff_subs = [];
        
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
        
        for xsub = args.g_sub
            
            %***************** target auc
            for xcond = 1:n_target
                aucs(xcond) = mvpa_ROC{xph}.roc_out.rand_auc{xcond}(xsub,2);
                criterions(xcond) = 0.5;
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
            fig_fname = fullfile(xoutput_dir, sprintf('plot_0.5cutoff_roc_baseline_%s', basename));
        else
            fig_fname = fullfile(xoutput_dir, sprintf('plot_0.5cutoff_roc_order_baseline_%s', basename));
        end
        %         savefig(xfig_norm, sprintf('%s.fig', fig_fname));
        set(gcf,'PaperUnits','inches','PaperPosition', [0 0 fig_w fig_h]/100)
        print('-djpeg', sprintf('%s.jpg',fig_fname), '-r100')
        saveas(xfig_norm, sprintf('%s.jpg',fig_fname), 'jpg')
        
%         close(xfig_norm);
    end
    
    %% ============= 1st LEVEL ROC PLOTTING
    %*************** random effect
    % mvpa_ROC{xph}.roc_out.rand_auc{xcate}(xsub,auc)
    clear xmatrix xheader
    
    x_tick = 0:1;
    x_lim  = [x_tick(1) x_tick(end)];
    y_tick = x_tick; y_lim   = x_lim;
    fig_w  = 2400; fig_h = 1200;
    lable_dist = 0.4;
    
    xmatrix = mvpa_ROC{xph}.roc_out.matrix;
    xheader = mvpa_ROC{xph}.roc_out.header;
    
    for xcate = 1:n_target
        xfig = figure;
        set(xfig, 'Position', [0 0 fig_w fig_h])
        
        n_col = 10;
        if mod(args.n_sub,n_col) ~= 0, n_row = fix(args.n_sub/n_col) + 1;
        else, n_row = fix(args.n_sub/n_col); end
        
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
        
        for xsub = 1:args.n_sub
            
            xx_fa  = getDATA(xmatrix, xheader, {'category','subject'}, {xcate, xsub}, findCol(xheader,{'x_fa'}));
            yy_hit = getDATA(xmatrix, xheader, {'category','subject'}, {xcate, xsub}, findCol(xheader,{'y_hit'}));
            
            subplot(n_row, n_col, xsub)
            plot(xx_fa, yy_hit); hold on;
            plot([0 1], [0 1],'--k')
            
            set(gca,'xlim', x_lim, 'ylim', y_lim);
            set(gca,'XTick', x_tick, 'YTick', y_tick);
            
            title(sprintf('%s (subject %s)', cond_names{xcate}, num2str(sub_id(xsub))));
            
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
%         fig_fname = fullfile(xoutput_dir, sprintf('plot_mvpa_roc_%s_%s', basename, cond_names{xcate}));
        
%         savefig(xfig, sprintf('%s.fig', fig_fname));
%         set(gcf,'PaperUnits','inches','PaperPosition', [0 0 fig_w fig_h]/100)
%         print('-djpeg', sprintf('%s.jpg',fig_fname), '-r100')
%         saveas(xfig, sprintf('%s.jpg',fig_fname), 'jpg')
%         
%         close(xfig);
        
    end
    
    %% ============= AUC DISTRIBUTION
    %*************** histogram for auc: probability density
    clear auc_subs
    
    filtered_subs = args.g_sub;%(~ismember(args.g_sub, auc_cutoff_subs));
    
    for xsub = 1:length(filtered_subs)
        auc_subs(xsub) = str2double(subject_list(filtered_subs(xsub)).name(end-2:end));
    end
    
    % plot properties
    n_bin    = 10;
    x_lim    = [0.4 1];
    y_lim    = [0 0.45];
    fig_rect = [0 0 1200 200];
    
    xfig     = figure;
    set(xfig, 'Position', fig_rect)
    
    for xcate = 1:n_target
        x_auc = mvpa_ROC{xph}.roc_out.rand_auc{xcate}(filtered_subs,2);
        
        subplot(1, n_target, xcate)
        histogram(x_auc, n_bin, 'Normalization', 'probability',...
            'DisplayStyle', 'bar', 'FaceColor', xcond_color{xcate}, 'FaceAlpha', 0.8, ...
            'LineWidth', 2, 'EdgeColor', xcond_color{xcate}, 'EdgeAlpha', 1);
        hold on;
        
        set(gca,'xlim', x_lim, 'ylim', y_lim);
        set(gca,'YTick', y_lim(1):0.05:y_lim(end))
        
        title(sprintf('Density: %s', cond_names{xcate}));
        xlabel('AUC');
        ylabel('density probability');
    end
    
    %*************** SAVE FIGURE
    fig_fname = fullfile(dirs.mvpa.group.auc{xph}, ...
        sprintf('plot_operation_auc_density_bin%s_%s', num2str(n_bin), basename));
    
%     savefig(xfig, sprintf('%s.fig', fig_fname));
    set(gcf,'PaperUnits','inches','PaperPosition',fig_rect/100)
    print('-djpeg', sprintf('%s.jpg',fig_fname), '-r100')
    saveas(xfig, sprintf('%s.jpg',fig_fname), 'jpg')
    
    close(xfig);
    
    %*************** save subject_num
    fname = fullfile(dirs.mvpa.group.home{xph}, ...
        sprintf('filtered_subj_nums_%s_%s.mat', basename, num2str(length(filtered_subs))));
    save(fname, 'filtered_subs')
    
    %% ============= AUC ADD
    if args.four_oper_regress
        %*************** mvpa_ROC{xph}.roc_out.rand_auc{xcate}(:,2)
        %*************** xcate_5: clear (3:4)
        %*************** xcate_6: removal operation (2:4), xcate_7: overall (1:4)
        cond_names = {'Maintain','RepCategory','TargetSuppress','GlobalClear'};
        
        for i = 6:9
            mvpa_ROC{xph}.roc_out.rand_auc{i} = zeros(args.n_sub, 2);
        end
        
        cat_auc = [];
        for xcate = 1:n_target
            cat_auc = horzcat(cat_auc, mvpa_ROC{xph}.roc_out.rand_auc{xcate}(:,2));
        end
        
        for xcate = 5:7, mvpa_ROC{xph}.roc_out.rand_auc{xcate}(:,1) = 1:args.n_sub; end
        
        mvpa_ROC{xph}.roc_out.rand_auc{5}(:,2) = mean(cat_auc(:,3:4),2);
        mvpa_ROC{xph}.roc_out.rand_auc{6}(:,2) = mean(cat_auc(:,2:4),2);
        mvpa_ROC{xph}.roc_out.rand_auc{7}(:,2) = mean(cat_auc(:,1:4),2);
        
        total_cate = 7;
        
        cond_names{5} = 'Clears';
        cond_names{6} = 'RemovalOperations';
        cond_names{7} = 'Overall';
    else
        %*************** mvpa_ROC{xph}.roc_out.rand_auc{xcate}(:,2)
        %*************** xcate_6: replace (2:3), xcate_7: clear (4:5)
        %*************** xcate_8: removal operation (2:5), xcate_9: overall (1:5)
        cond_names = {'Maintain','RepCategory','RepSubcate','TargetSuppress','GlobalClear'};
        
        for i = 6:9
            mvpa_ROC{xph}.roc_out.rand_auc{i} = zeros(args.n_sub, 2);
        end
        
        cat_auc = [];
        for xcate = 1:n_target
            cat_auc = horzcat(cat_auc, mvpa_ROC{xph}.roc_out.rand_auc{xcate}(:,2));
        end
        
        for xcate = 6:9, mvpa_ROC{xph}.roc_out.rand_auc{xcate}(:,1) = 1:args.n_sub; end
        
        mvpa_ROC{xph}.roc_out.rand_auc{6}(:,2) = mean(cat_auc(:,2:3),2);
        mvpa_ROC{xph}.roc_out.rand_auc{7}(:,2) = mean(cat_auc(:,4:5),2);
        mvpa_ROC{xph}.roc_out.rand_auc{8}(:,2) = mean(cat_auc(:,1:4),2);
        mvpa_ROC{xph}.roc_out.rand_auc{9}(:,2) = mean(cat_auc(:,1:5),2);
        
        total_cate = 9;
        
        cond_names{6} = 'Replaces';
        cond_names{7} = 'Clears';
        cond_names{8} = 'RemovalOperations';
        cond_names{9} = 'Overall';
    end
    
    %% *************** total auc plot
    xcate         = total_cate;%total
    x_auc         = mvpa_ROC{xph}.roc_out.rand_auc{xcate}(filtered_subs,2);
    xmean         = mean(x_auc);
    xse           = std(x_auc)/sqrt(args.n_sub);
    % xse_cut       = args.se_cutoff;
    % xcutoff       = xmean - (xse_cut * xse);
    % cut_subs      = auc_subs(x_auc<=xcutoff);
    % find_cut      = find(ismember(auc_subs,cut_subs));%position of the cut_subs
    
    % se_filtered_subs = auc_subs(~ismember(auc_subs, cut_subs));
    
    xcolor_auc    = [75 172 198]/255;
    xcolor_cut    = [247 106 120]/255;
    x_lim         = [0 length(filtered_subs)+1];
    y_lim         = [0 1];
    fig_rect      = [0 0 1200 400];
    
    xfig          = figure;
    set(xfig, 'Position', fig_rect)
    
    bar(1:length(filtered_subs), x_auc, 'FaceColor', xcolor_auc, ...
        'EdgeColor',[0 0 0],'LineWidth',0.5); hold on
    
    if args.se_filtered
        cut_bar = zeros(1,args.n_sub);
        cut_bar(find_cut) = x_auc(find_cut);
        
        bar(1:args.n_sub, cut_bar,'FaceColor', xcolor_cut, ...
            'EdgeColor',[0 0 0],'LineWidth',0.5);
        
        plot([x_lim(1) x_lim(2)], [xcutoff xcutoff], '--r','LineWidth', 1);
        
        text(x_lim(1)+0.5, y_lim(2) - 0.1, ...
            sprintf('cutoff = mean - %d se, cutoff N = %s', ...
            xse_cut, num2str(length(cut_subs))), 'Color', 'r',...
            'FontSize', 15, 'FontWeight', 'bold');
    end
    
    plot([x_lim(1) x_lim(2)], [xmean xmean], '-b','LineWidth', 1);
    plot([x_lim(1) x_lim(2)], [xmean-xse xmean-xse], '--b','LineWidth', 1);
    plot([x_lim(1) x_lim(2)], [xmean+xse xmean+xse], '--b','LineWidth', 1);
    
    set(gca,'xlim', x_lim, 'ylim', y_lim);
    set(gca,'XTick', 1:args.n_sub)
    set(gca,'XTickLabel', auc_subs)
    set(gca,'YTick', y_lim(1):0.1:y_lim(end))
    
    title(sprintf('Total AUC for operation MVPA (N=%s)', num2str(length(filtered_subs))),...
        'FontSize', 15);
    xlabel('subject');
    ylabel('AUC: operation MVPA');
    
    %*************** SAVE FIGURE
    if args.se_filtered
        fig_fname = fullfile(dirs.mvpa.group.auc{xph}, ...
            sprintf('plot_operation_auc_total_cutoff%s_%s', num2str(args.se_cutoff), basename));
    else
        fig_fname = fullfile(dirs.mvpa.group.auc{xph}, ...
            sprintf('plot_operation_auc_total_%s', basename));
    end
    
    savefig(xfig, sprintf('%s.fig', fig_fname));
    set(gcf,'PaperUnits','inches','PaperPosition',fig_rect/100)
    print('-djpeg', sprintf('%s.jpg',fig_fname), '-r100')
    saveas(xfig, sprintf('%s.jpg',fig_fname), 'jpg')
    
    close(xfig);
    
    %% ============= AUC cross-correlation
    % mvpa_ROC{xph}.roc_out.rand_auc{xcond}(args.g_sub,auc)
    % cond_names: 1:5: 'Maintain','RepCategory','RepSubcate','TargetSuppress','GlobalClear'
    % 6:'Replaces',7:'Clears'
    clear heatmap_r heatmap_p it_subs
    
    n_auc    = 7;
    it_subs  = filtered_subs;
    
    x_lim    = [0.4 1]; y_lim = x_lim;
    fig_rect = [0 0 2000 2100];
    
    x_fig    = figure;
    set(x_fig, 'Position', fig_rect)
    
    for x = 1:n_auc % 1:(n_auc-1)
        for y = 1:n_auc % x:n_auc
            clear xx_auc yy_auc ax
            
            xx_auc = mvpa_ROC{xph}.roc_out.rand_auc{x}(it_subs, 2);
            yy_auc = mvpa_ROC{xph}.roc_out.rand_auc{y}(it_subs, 2);
            [xr, xpvalue] = corrcoef(xx_auc, yy_auc);
            
            heatmap_r(x, y) = xr(1,2);
            heatmap_p(x, y) = xpvalue(1,2);
            
            %*************** plot
            i = y + (n_auc * (x-1));
            
            ax = subplot(n_auc,n_auc, i); hold on;
            s  = scatter(ax, xx_auc, yy_auc,'o');
            s.MarkerEdgeColor = [0.4 0.4 0.4];
            s.SizeData  = 3;
            s.LineWidth = 0.5;
            
            L = lsline(ax);
            L.LineWidth = 1;
            L.Color = 'r';
            
            P = polyfit(xx_auc, yy_auc, 1);%linear fit
            y_fit = P(1) * xx_auc + P(2);
            
            text(x_lim(1)+0.05, y_lim(2) - 0.1, ...
                sprintf('r2 = %1.2f\n p = %1.2f', xr(1,2)^2, xpvalue(1,2)),...
                'FontSize', 8, 'FontWeight', 'bold');
            %         text(x_lim(1)+0.05, y_lim(2) - 0.2, ...
            %             sprintf('y(lin) = %1.2fx + %1.2f', P(1), P(2)),...
            %             'FontSize', font_size, 'FontWeight', 'bold');
            %
            set(gca,'xlim', x_lim, 'ylim', y_lim);
            set(gca,'YTick', y_lim(1):0.2:y_lim(end))
            
            %         title(sprintf('%s & %s (N=%s)', cond_names{x}, cond_names{y}, num2str(length(args.g_sub))));
            
            xlabel(cond_names{x});
            ylabel(cond_names{y});
        end
    end
    
    %*************** save fig
    fig_fname = fullfile(dirs.mvpa.group.auc{xph}, ...
        sprintf('plot_cross_correlation_operation_auc_%s_n%s', ...
        basename, num2str(length(it_subs))));
    
%     savefig(x_fig, sprintf('%s.fig', fig_fname));
    set(gcf,'PaperUnits','inches','PaperPosition',fig_rect/100)
    print('-djpeg', sprintf('%s.jpg',fig_fname), '-r100')
    saveas(x_fig, sprintf('%s.jpg',fig_fname), 'jpg')
    
    close(x_fig);
    
    %% --------- heatmap figure
    % heatmap_r, heatmap_p(x, y)
    for i = 1:n_auc, heat_names{i} = cond_names{i}; end
    
    fig_rect = [0 0 1000 900];
    
    heatmap_fig = figure;
    set(heatmap_fig, 'Position', fig_rect)
    
    imagesc(heatmap_r);
    colormap(flipud(pink)); %spring pink
    h = colorbar;
    set(h,'Ylim',[0.7 1])
    hold on
    
    for x = 1:n_auc
        for y = 1:n_auc
            text(0.9 + (x-1), 1 + (y-1), ...
                sprintf('%1.2f', heatmap_r(x,y)),...
                'Color','r','FontSize',10,'FontWeight','bold')
        end
    end
    
    xlabel('similarity matrix of operations','FontSize',15,'FontWeight','bold');
    
    xticks = 1:7;
    set(gca,'XTickLabel', heat_names)
    set(gca,'xaxisLocation','top')
    set(gca,'YTickLabel', heat_names)
    
    %*************** save fig
    fig_fname = fullfile(dirs.mvpa.group.auc{xph}, ...
        sprintf('plot_heatmap_operation_auc_%s_n%s', ...
        basename, num2str(length(it_subs))));
    
%     savefig(heatmap_fig, sprintf('%s.fig', fig_fname));
    set(gcf,'PaperUnits','inches','PaperPosition',fig_rect/100)
    print('-djpeg', sprintf('%s.jpg',fig_fname), '-r100')
    saveas(heatmap_fig, sprintf('%s.jpg',fig_fname), 'jpg')
    
%     close(heatmap_fig);
    
    %% ============= AUC & INDIVIDUAL DEFFERENCE
    if args.se_filtered
        it_auc_subs = se_filtered_subs;
    else
        it_auc_subs = auc_subs;
    end
    
    %*************** load difference matrix
    xprotocol_xls      = fullfile(dirs.protocols, 'individual_matrix.xlsx');
    [xmatrix, xheader] = xlsread(xprotocol_xls);
    
    %*************** filtering subjects with existing behavioral data
    behav_subs  = xmatrix(:, findCol(xheader, {'subject'}));
    xbehav_subs = it_auc_subs(ismember(it_auc_subs, behav_subs));
    xsternberg_sub = []; xstroop_sub = [];
    xtotal_personality_sub = []; zWBSIs_sub = []; zPSWQtotal_sub = []; zRRSbrood_sub = [];
    
    for xsub = 1:length(xbehav_subs)
        it_subs = xbehav_subs(xsub);
        
        %*************** sternberg
        xstern = getDATA(xmatrix, xheader, {'subject'}, {it_subs}, ...
            findCol(xheader, {'100 intrusion costs','2000 intrusion costs'}));
        
        if sum(isnan(xstern))==0
            xsternberg_sub = vertcat(xsternberg_sub, it_subs); %#ok<*NASGU>
        end
        
        %*************** stroop
        xstroop = getDATA(xmatrix, xheader, {'subject'}, {it_subs}, ...
            findCol(xheader, {'incong25 stroop effect','incong75 stroop effect'}));
        
        if sum(isnan(xstroop))==0
            xstroop_sub = vertcat(xstroop_sub, it_subs);
        end
        
        %*************** zWBSIs
        xzWBSIs = getDATA(xmatrix, xheader, {'subject'}, {it_subs}, ...
            findCol(xheader, {'zWBSIs'}));
        
        if sum(isnan(xzWBSIs))==0
            zWBSIs_sub = vertcat(zWBSIs_sub, it_subs);
        end
        
        %*************** zPSWQtotal
        zPSWQtotal = getDATA(xmatrix, xheader, {'subject'}, {it_subs}, ...
            findCol(xheader, {'zPSWQtotal'}));
        
        if sum(isnan(zPSWQtotal))==0
            zPSWQtotal_sub = vertcat(zPSWQtotal_sub, it_subs);
        end
        
        %*************** zRRSbrood
        zRRSbrood = getDATA(xmatrix, xheader, {'subject'}, {it_subs}, ...
            findCol(xheader, {'zRRSbrood'}));
        
        if sum(isnan(zRRSbrood))==0
            zRRSbrood_sub = vertcat(zRRSbrood_sub, it_subs);
        end
        
        %*************** personality total
        xpersonality = getDATA(xmatrix, xheader, {'subject'}, {it_subs}, ...
            findCol(xheader, {'personality score'}));
        
        if sum(isnan(xpersonality))==0
            xtotal_personality_sub = vertcat(xtotal_personality_sub, it_subs);
        end
    end
    
    xpersonality_sub{1} = zWBSIs_sub;
    xpersonality_sub{2} = zPSWQtotal_sub;
    xpersonality_sub{3} = zRRSbrood_sub;
    xpersonality_sub{4} = xtotal_personality_sub;
    
    %% *************** STERNBERG CORRELATION
    %***************** decreased_intrusion_costs: 100-2000 intrusion costs
    clear x_auc
    
    fig_rect  = [0 0 1200 300];
    y_lim     = [-0.4 1.2];
    
    xcolor{1} = 'r';
    xcolor{2} = 'b';
    xcolor{3} = 'k';
    
    t_sub     = xmatrix(:, findCol(xheader, {'subject'}));
    xunit     = ismember(t_sub, xsternberg_sub);
    it_sub_id = auc_subs(ismember(auc_subs, xsternberg_sub));
    it_subs    = ismember(sub_id, it_sub_id);
    
    for xcate = 1:length(cond_names)
        
        x_auc = mvpa_ROC{xph}.roc_out.rand_auc{xcate}(it_subs,2);
        
        y_removal{1} = xmatrix(xunit, findCol(xheader, {'100 intrusion costs'}));
        y_removal{2} = xmatrix(xunit, findCol(xheader, {'2000 intrusion costs'}));
        y_removal{3} = xmatrix(xunit, findCol(xheader, {'decreased intrusion costs'}));
        
        %***************** plotting
        xfig  = figure;
        set(xfig, 'Position', fig_rect)
        
        x_lim = [min(x_auc) max(x_auc)];
        
        for i = 1:3
            clear ax r pvalue L P y_fit
            [r, pvalue] = corrcoef(x_auc, y_removal{i});
            
            ax{i} = subplot(1,3, i); hold on;
            scatter(ax{i}, x_auc, y_removal{i}, xcolor{i});
            L = lsline(ax{i});
            L.LineWidth = 2;
            L.Color = xcolor{i};
            
            P = polyfit(x_auc, y_removal{i}, 1);%linear fit
            y_fit = P(1) * x_auc + P(2);
            
            text(x_lim(1)+0.05, y_lim(2) - 0.1, ...
                sprintf('r2 = %1.4f, p = %1.2f', r(1,2)^2, pvalue(1,2)),...
                'FontSize', font_size, 'FontWeight', 'bold');
            text(x_lim(1)+0.05, y_lim(2) - 0.2, ...
                sprintf('y(lin) = %1.2fx + %1.2f', P(1), P(2)),...
                'FontSize', font_size, 'FontWeight', 'bold');
            
            set(gca,'xlim', x_lim, 'ylim', y_lim);
            set(gca,'YTick', y_lim(1):0.2:y_lim(end))
            
            if args.se_filtered
                title(sprintf('%s (filtered N=%s, se: %s)', ...
                    cond_names{xcate}, num2str(length(xsternberg_sub)), num2str(args.se_cutoff)));
            else
                title(sprintf('%s (N=%s)', cond_names{xcate}, num2str(length(xsternberg_sub))));
            end
            
            xlabel('AUC');
            ylabel(xheader{i + 1});
        end
        
        %% *************** SAVE FIGURE
        if args.se_filtered
            fig_fname = fullfile(dirs.mvpa.group.auc{xph}, ...
                sprintf('plot_sternberg_operation_auc_cutoff%s_%s_%s', ...
                num2str(args.se_cutoff), cond_names{xcate}, basename));
        else
            fig_fname = fullfile(dirs.mvpa.group.auc{xph}, ...
                sprintf('plot_sternberg_operation_auc_%s_%s', cond_names{xcate}, basename));
        end
        
        savefig(xfig, sprintf('%s.fig', fig_fname));
        set(gcf,'PaperUnits','inches','PaperPosition',fig_rect/100)
        print('-djpeg', sprintf('%s.jpg',fig_fname), '-r100')
        saveas(xfig, sprintf('%s.jpg',fig_fname), 'jpg')
        
%         close(xfig);
        
    end
    
    %% *************** STROOP CORRELATION
    %***************** adapted suppress: 100-2000 intrusion costs
    clear x_auc y_removal
    fig_rect = [0 0 1200 300];
    y_lim    = [-0.2 0.8];
    
    xunit     = ismember(t_sub, xstroop_sub);
    it_sub_id = auc_subs(ismember(auc_subs, xstroop_sub));
    it_subs    = ismember(sub_id, it_sub_id);
    
    for xcate = 1:length(cond_names)
        
        x_auc = mvpa_ROC{xph}.roc_out.rand_auc{xcate}(it_subs,2);
        
        y_removal{1} = xmatrix(xunit, findCol(xheader, {'incong25 stroop effect'}));
        y_removal{2} = xmatrix(xunit, findCol(xheader, {'incong75 stroop effect'}));
        y_removal{3} = xmatrix(xunit, findCol(xheader, {'decreased stroop effect'}));
        
        %***************** plotting
        xfig  = figure;
        set(xfig, 'Position', fig_rect)
        
        x_lim = [min(x_auc) max(x_auc)];
        
        for i = 1:3
            clear r pvalue L P y_fit
            [r, pvalue] = corrcoef(x_auc, y_removal{i});
            
            ax{i} = subplot(1,3, i); hold on;
            scatter(ax{i}, x_auc, y_removal{i}, xcolor{i});
            L = lsline(ax{i});
            L.LineWidth = 2;
            L.Color = xcolor{i};
            
            P = polyfit(x_auc, y_removal{i}, 1);%linear fit
            y_fit = P(1) * x_auc + P(2);
            
            text(x_lim(1)+0.05, y_lim(2) - 0.1, ...
                sprintf('r2 = %1.4f, p = %1.2f', r(1,2)^2, pvalue(1,2)),...
                'FontSize', font_size, 'FontWeight', 'bold');
            text(x_lim(1)+0.05, y_lim(2) - 0.2, ...
                sprintf('y(lin) = %1.2fx + %1.2f', P(1), P(2)),...
                'FontSize', font_size, 'FontWeight', 'bold');
            
            set(gca,'xlim', x_lim, 'ylim', y_lim);
            set(gca,'YTick', y_lim(1):0.1:y_lim(end))
            
            if args.se_filtered
                title(sprintf('%s (filtered N=%s, se: %s)', ...
                    cond_names{xcate}, num2str(length(xstroop_sub)), num2str(args.se_cutoff)))
            else
                title(sprintf('%s (N=%s)', cond_names{xcate}, num2str(length(xstroop_sub))));
            end
            
            xlabel('AUC');
            ylabel(xheader{i + 4});
        end
        
        %% *************** SAVE FIGURE
        if args.se_filtered
            fig_fname = fullfile(dirs.mvpa.group.auc{xph}, ...
                sprintf('plot_stroopeffect_operation_auc_cutoff%s_%s_%s', ...
                num2str(args.se_cutoff), cond_names{xcate}, basename));
        else
            fig_fname = fullfile(dirs.mvpa.group.auc{xph}, ...
                sprintf('plot_stroopeffect_operation_auc_%s_%s', cond_names{xcate}, basename));
        end
        
        savefig(xfig, sprintf('%s.fig', fig_fname));
        set(gcf,'PaperUnits','inches','PaperPosition',fig_rect/100)
        print('-djpeg', sprintf('%s.jpg',fig_fname), '-r100')
        saveas(xfig, sprintf('%s.jpg',fig_fname), 'jpg')
        
        close(xfig);
    end
    
    %% *************** PERSONALITY CORRELATION
    %***************** adapted suppress: 100-2000 intrusion costs
    clear x_auc y_removal
    y_lim    = [-2.5 2.5];
    fig_rect = [0 0 1600 300];
    
    t_sub  = xmatrix(:, findCol(xheader, {'subject'}));
    
    for xcate = 1:length(cond_names)
        
        xfig  = figure;
        set(xfig, 'Position', fig_rect)
        
        for i = 1:4
            clear r pvalue L P y_fit
            y_value   = xheader{i + 7};
            xunit     = ismember(t_sub, xpersonality_sub{i});
            it_sub_id = auc_subs(ismember(auc_subs, xpersonality_sub{i}));
            it_subs    = ismember(sub_id, it_sub_id);
            
            %***************** correlation
            x_auc        = mvpa_ROC{xph}.roc_out.rand_auc{xcate}(it_subs,2);
            y_removal{i} = xmatrix(xunit, findCol(xheader, {y_value}));
            
            [r, pvalue] = corrcoef(x_auc, y_removal{i});
            
            %***************** plotting
            x_lim = [min(x_auc) max(x_auc)];
            
            ax{i} = subplot(1,4, i); hold on;
            scatter(ax{i}, x_auc, y_removal{i}, xcolor{fix(i/4)+1});
            L = lsline(ax{i});
            L.LineWidth = 2;
            L.Color = xcolor{fix(i/4)+1};
            
            P = polyfit(x_auc, y_removal{i}, 1);%linear fit
            y_fit = P(1) * x_auc + P(2);
            
            text(x_lim(1)+0.05, y_lim(2) - 0.2, ...
                sprintf('r2 = %1.4f, p = %1.4f', r(1,2)^2, pvalue(1,2)),...
                'FontSize', font_size, 'FontWeight', 'bold');
            text(x_lim(1)+0.05, y_lim(2) - 0.5, ...
                sprintf('y(lin) = %1.2fx + %1.4f', P(1), P(2)),...
                'FontSize', font_size, 'FontWeight', 'bold');
            
            set(gca,'xlim', x_lim, 'ylim', y_lim);
            set(gca,'YTick', y_lim(1):0.5:y_lim(end))
            
            if args.se_filtered
                title(sprintf('%s (filtered N=%s, se: %s)', ...
                    cond_names{xcate}, num2str(length(xpersonality_sub{i})), num2str(args.se_cutoff)))
            else
                title(sprintf('%s (N=%s)', cond_names{xcate}, num2str(length(xpersonality_sub{i}))));
            end
            
            xlabel('AUC');
            ylabel(y_value);
        end
        
        %% *************** SAVE FIGURE
        if args.se_filtered
            fig_fname = fullfile(dirs.mvpa.group.auc{xph}, ...
                sprintf('plot_personality_operation_auc_cutoff%s_%s_%s', ...
                num2str(args.se_cutoff), cond_names{xcate}, basename));
        else
            fig_fname = fullfile(dirs.mvpa.group.auc{xph}, ...
                sprintf('plot_personality_operation_auc_%s_%s', cond_names{xcate}, basename));
        end
        
        savefig(xfig, sprintf('%s.fig', fig_fname));
        set(gcf,'PaperUnits','inches','PaperPosition',fig_rect/100)
        print('-djpeg', sprintf('%s.jpg',fig_fname), '-r100')
        saveas(xfig, sprintf('%s.jpg',fig_fname), 'jpg')
        
        close(xfig);
    end
    
    %% temp for figures
     %% *************** PERSONALITY CORRELATION
    %***************** adapted suppress: 100-2000 intrusion costs
    clear x_auc y_removal
    y_lim    = [-2.5 2.5];
    fig_rect = [0 0 1500 288];
    
    t_sub  = xmatrix(:, findCol(xheader, {'subject'}));
    xfig  = figure;
    set(xfig, 'Position', fig_rect)
    i = 2;
    
    for xcate = 1:4%length(cond_names)
        
        clear r pvalue L P y_fit y_removal
        y_value   = xheader{i + 7};
        xunit     = ismember(t_sub, xpersonality_sub{i});
        it_sub_id = auc_subs(ismember(auc_subs, xpersonality_sub{i}));
        it_subs    = ismember(sub_id, it_sub_id);
        
        %***************** correlation
        x_auc     = mvpa_ROC{xph}.roc_out.rand_auc{xcate}(it_subs,2);
        y_removal = xmatrix(xunit, findCol(xheader, {y_value}));
        
        [r, pvalue] = corrcoef(x_auc, y_removal);
        
        %***************** plotting
        x_lim = [min(x_auc) max(x_auc)];
        clear ax L
        ax = subplot(1,4, xcate); hold on;
        scatter(ax, x_auc, y_removal, 'MarkerFaceColor', xcond_color{xcate}',...
            'MarkerEdgeColor', xcond_color{xcate});
        L = lsline(ax);
        L.LineWidth = 2;
        L.Color = xcond_color{xcate};
        
        P = polyfit(x_auc, y_removal, 1);%linear fit
        y_fit = P(1) * x_auc + P(2);
        
        text(x_lim(1)+0.05, y_lim(2) - 0.2, ...
            sprintf('r2 = %1.4f, p = %1.4f', r(1,2)^2, pvalue(1,2)),...
            'FontSize', font_size, 'FontWeight', 'bold');
        text(x_lim(1)+0.05, y_lim(2) - 0.5, ...
            sprintf('y(lin) = %1.2fx + %1.4f', P(1), P(2)),...
            'FontSize', font_size, 'FontWeight', 'bold');
        
        set(gca,'xlim', x_lim, 'ylim', y_lim);
        set(gca,'YTick', y_lim(1):0.5:y_lim(end))
        
        if args.se_filtered
            title(sprintf('%s (filtered N=%s, se: %s)', ...
                cond_names{xcate}, num2str(length(xpersonality_sub{i})), num2str(args.se_cutoff)))
        else
            title(sprintf('%s (N=%s)', cond_names{xcate}, num2str(length(xpersonality_sub{i}))));
        end
        
        xlabel('AUC');
        ylabel(y_value);
        
        %% *************** SAVE FIGURE
        fig_fname = fullfile(dirs.mvpa.group.auc{xph}, ...
            sprintf('tmp_plot_personality_operation_auc_%s', basename));
        
        savefig(xfig, sprintf('%s.fig', fig_fname));
        set(gcf,'PaperUnits','inches','PaperPosition',fig_rect/100)
        print('-djpeg', sprintf('%s.jpg',fig_fname), '-r100')
        saveas(xfig, sprintf('%s.jpg',fig_fname), 'jpg')
        
%         close(xfig);
        
    end
    
    %% ============= MVPA ACCURACY CONFUSION MATRIX
    cond_names = {'Maintain','RepCategory','RepSubcate','Suppress','Clear'};
    
    xguess_csv  = fullfile(dirs.mvpa.group.out{xph}, ...
        sprintf('grp_acc_confusion_matrix_%s_n%s.csv', basename, num2str(length(xsub_groups))));
    xtable      = readtable(xguess_csv);
    xacc_matrix = table2array(xtable);
    
    %% *************** confusion matrix plot
    fig_rect = [0 0 1200 450];
    
    heatmap_fig = figure;
    set(heatmap_fig, 'Position', fig_rect)
    
    subplot(1,2,1)
    imagesc(xacc_matrix); colormap(flipud(pink)); %spring pink
    h = colorbar;
    set(h,'Ylim',[0 0.6])
    hold on
    
    for x = 1:n_target
        for y = 1:n_target
            text(0.85 + (x-1), 0.98 + (y-1), ...
                sprintf('%1.2f', xacc_matrix(x,y)),'Color','r','FontSize',(45/n_target),'FontWeight','bold')
        end
    end
    
    x_lim  = [0 n_target] + 0.5;
    x_tick = 1:n_target;
    set(gca,'xlim', x_lim, 'ylim', x_lim);
    set(gca,'XTick', x_tick)
    set(gca,'YTick', x_tick)
    
    set(gca,'XTickLabel', cond_names,'FontSize',6)
    set(gca,'xaxisLocation','top')
    set(gca,'YTickLabel', cond_names,'FontSize',6)
    if strcmp(args.level, 'subcategory'), xtickangle(90); end
    
    title(sprintf('WM operation (N=%s)', num2str(length(xsub_groups))),'FontSize',10);
    xlabel('target operation','FontSize',10);
    ylabel('decode operation','FontSize',10);
    
    %*************** auc plot
    clear aucs criterions
    
    y_lim = [0.4 0.9];
    y_tick = y_lim(1):0.1:y_lim(2);
    
    t_all_aucs = []; t_all_baselines = [];
    for xcond = 1:n_target
        aucs.mean(xcond)       = mean(mvpa_ROC{xph}.roc_out.rand_auc{xcond}(xsub_groups,2));
        aucs.se(xcond)         = std(mvpa_ROC{xph}.roc_out.rand_auc{xcond}(xsub_groups,2))/sqrt(length(xsub_groups));
        criterions.mean(xcond) = mean(mvpa_ROC{xph}.baseline_criterion{xcond}(xsub_groups,1));
        criterions.se(xcond)   = std(mvpa_ROC{xph}.baseline_criterion{xcond}(xsub_groups,1))/sqrt(length(xsub_groups));
        
        t_all_aucs      = horzcat(t_all_aucs, mvpa_ROC{xph}.roc_out.rand_auc{xcond}(xsub_groups,2));
        t_all_baselines = horzcat(t_all_baselines, mvpa_ROC{xph}.baseline_criterion{xcond}(xsub_groups,1));
    end
    
    all_aucs      = mean(t_all_aucs,2);
    all_baselines = mean(t_all_baselines,2);
    
    subplot(1,2,2)
    
    for i = 1:n_target
        errorbar(i, aucs.mean(i), aucs.se(i), 'o', 'Color', xcond_color{i})
        hold on
    end
    
    for i = 1:n_target
        errorbar(i, criterions.mean(i), criterions.se(i), 'o', 'Color', xbase_color)
    end
    
    xlegend     = cond_names;
    xlegend{length(cond_names)+1} = 'baseline';
    
    lg          = legend(xlegend);
    lg.Location = 'BestOutside';
    legend(xlegend,'AutoUpdate','off')
    
    set(gca,'xlim', x_lim, 'ylim', y_lim);
    set(gca,'XTick', x_tick)
    set(gca,'YTick', y_tick)
    set(gca,'XTickLabel', cond_names,'FontSize',6)
    
    title(sprintf('WM operation ACU (N=%s)', num2str(length(xsub_groups))),'FontSize',10);
    xlabel('target category','FontSize',10);
    ylabel('AUC','FontSize',10);
    grid on
    
    %*************** mean/se
    text(1, y_lim(1) + 0.2, sprintf('auc: M=%4.4f, SE=%4.4f',...
        mean(all_aucs), std(all_aucs)/sqrt(length(xsub_groups))));
    
    text(1, y_lim(1) + 0.15, sprintf('baseline: M=%4.4f, SE=%4.4f',...
        mean(all_baselines), std(all_baselines)/sqrt(length(xsub_groups))));
    
    %*************** save fig
    fig_fname = fullfile(xoutput_dir, sprintf('plot_accuracy_roc_%s_n%s', basename, num2str(length(xsub_groups))));
    
    savefig(heatmap_fig, sprintf('%s.fig', fig_fname));
    set(gcf,'PaperUnits','inches','PaperPosition', fig_rect/100)
    print('-djpeg', sprintf('%s.jpg',fig_fname), '-r100')
    saveas(heatmap_fig, sprintf('%s.jpg',fig_fname), 'jpg')
    
    close(heatmap_fig);
    
    %% ============= violin plot for AUC
    xauc_data = []; xauc_baseline = []; xcolors = [];
    for xcond = 1:n_target
        xauc_data(:,xcond)     = mvpa_ROC{xph}.roc_out.rand_auc{xcond}(xsub_groups,2);
        xauc_baseline(:,xcond) = mvpa_ROC{xph}.baseline_criterion{xcond}(xsub_groups,1);
        
        xcolors(xcond, :) = xcond_color{xcond};
    end
    
    fig_rect = [0 0 700 1000];
    
    violin_fig = figure;
    set(violin_fig, 'Position', fig_rect)
    
    y_lim  = [0.3 1.2];
    y_tick = y_lim(1):0.1:y_lim(2);
    x_tick = 1:n_target;
    
    %*************** legend
    h    = zeros(n_target+1, 1);
    for i=1:n_target
        h(i) = plot(NaN, NaN, 'Color', xcolors(i, :),'LineWidth',5);
        hold on;
    end
    h(n_target+1) = plot(NaN, NaN, 'Color', xbase_color,'LineWidth',5);
    xlegend       = cond_names;
    xlegend{length(cond_names)+1} = 'baseline';
    lg            = legend(xlegend);
    lg.Location   = 'BestOutside';
    legend(xlegend,'AutoUpdate','off')
    
    %*************** violin
    violin(xauc_data,'facecolor',xcolors,'xlabel', cond_names,...
        'edgecolor','','medc','','mc','k','plotlegend',''); hold on
    % violin(xauc_baseline,'facecolor',xbase_color,'xlabel', cond_names,...
    %     'edgecolor','','medc','','mc','k','plotlegend','');
    
    plot([0 6], [0.5 0.5],':k','LineWidth',2);
    
    title(sprintf('Operation ACU (N=%s)', num2str(length(xsub_groups))));
    
    set(gca,'ylim', y_lim)
    set(gca,'YTick', y_tick)
    set(gca,'XTick', x_tick, 'XTickLabel', cond_names)
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
    
    %% ============= box plot for AUC
    xauc_data = []; xauc_baseline = []; xcolors = [];
    for xcond = 1:n_target
        xauc_data(:,xcond)     = mvpa_ROC{xph}.roc_out.rand_auc{xcond}(xsub_groups,2);
        xauc_baseline(:,xcond) = mvpa_ROC{xph}.baseline_criterion{xcond}(xsub_groups,1);
        
        xcolors(xcond, :) = xcond_color{xcond};
    end
    
    fig_rect = [0 0 700 1000];
    
    box_fig = figure;
    set(box_fig, 'Position', fig_rect)
    
    y_lim  = [0.5 1.0];
    y_tick = y_lim(1):0.1:y_lim(2);
    x_tick = 1:n_target;
    
    %*************** legend
    h    = zeros(n_target+1, 1);
    for i=1:n_target
        h(i) = plot(NaN, NaN, 'Color', xcolors(i, :),'LineWidth',5);
        hold on;
    end
    h(n_target+1) = plot(NaN, NaN, 'Color', xbase_color,'LineWidth',5);
    xlegend       = cond_names;
    lg            = legend(xlegend);
    lg.Location   = 'BestOutside';
    legend(xlegend,'AutoUpdate','off')
    
    %*************** box plot
    boxplot(xauc_data, 'Notch','on','Colors',xcolors,'BoxStyle','outline','Widths',0.5,'Symbol','o')
    
    plot([0 6], [0.5 0.5],':k','LineWidth',2);
    
    title(sprintf('Operation ACU (N=%s)', num2str(length(xsub_groups))));
    
    set(gca,'ylim', y_lim)
    set(gca,'YTick', y_tick)
    set(gca,'XTick', x_tick, 'XTickLabel', cond_names)
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
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%============= WRITING STATS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    xsub_groups = args.filtered_subs;
    
    %%============= 2nd LEVEL SUBJECT DATA: AUC / CLASSIFIER
    %*************** random effect
    % mvpa_ROC{xph}.roc_out.rand_auc{xcate}(xsub,auc)
    if args.four_oper_regress
        it_conds = [1 2 4 5];
        cond_names = {'Maintain','RepCategory','Suppress','Clear'};
    else
        it_conds = 1:5;
        cond_names = {'Maintain','RepCategory','RepSubcate','Suppress','Clear'};
    end
    
    n_condition = length(cond_names);
    n_subs      = length(xsub_groups);
    xmeas_names = {'classifier_auc','classifier_auc_baseline','classifier_accuracy','classifier_evidence'};
    
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
            xnames{xcate + 2} = cond_names{xcate};
        end
        
        %*************** data matrix
        xmatrix(:,1) = xsub_groups';
        xmatrix(:,2) = sub_id_filtered';
        
        if (xmeas == 1) || (xmeas == 2)
            for xcate = 1:n_condition
                if xmeas == 1, xdat = mvpa_ROC{xph}.roc_out.rand_auc{xcate}(:,2);
                else,          xdat = ones(args.n_sub, 1) * 0.5;%  mvpa_ROC{xph}.baseline_criterion{xcate}(:,1); 
                end
                
                xmatrix = horzcat(xmatrix, xdat(xsub_groups));
            end
            
            xchance = 0.5;%for auc
            
        elseif (xmeas == 3) || (xmeas == 4)
            if xmeas == 3, xcsv = fullfile(dirs.mvpa.group.out{xph}, sprintf('grp_acc_matrix_%s_n%s.csv', basename, num2str(n_subs)));
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
            for xcate = 1:(n_target-1)
                xmatrix = horzcat(xmatrix, mvpa_ROC{xph}.roc_out.rand_auc_ex_rest{xcate}(xsub_groups,2));
                xnames{xcate + 2 + n_condition} = sprintf('%s_ex_rest', cond_names{xcate});
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
        
    %% ============= COMPARING: ACCURACY/EVIDENCE
    clear xnames xconfusion
    xmeas_names = {'accuracy','evidence'};
    targ_array  = 1:n_condition;
    
    fprintf(xout_txt, '\n\n####################################################################\n');
    fprintf(xout_txt, '#### paired T test\n');
    fprintf(xout_txt, '####################################################################\n\n');

    for xmeas = 1:length(xmeas_names)%1_auc, 2_auc_baseline, 3_4_classifier accuracy, 4_classifier evidence
        clear xtable xheader xtarget xguess
        
        fprintf(xout_txt, '\n=========================================\n');
        fprintf(xout_txt, '======== %s\n\n', xmeas_names{xmeas});
        
        %*************** data matrix
        xcsv = fullfile(dirs.mvpa.group.out{xph}, ...
            sprintf('grp_%s_%s_n%s.csv', xmeas_names{xmeas}, basename, num2str(n_subs)));
        
        xtable  = readtable(xcsv);
        xheader = xtable.Properties.VariableNames;
        
        xtarget = unique(xtable.Row)';
        xguess  = xtable.Properties.VariableNames(~ismember(xtable.Properties.VariableNames, 'Row'));
        
        for xtarg = 1:length(xtarget)
            xunit = ismember(xtable.Row, xtarget{xtarg});
            xconfusion{xtarg} = xtable(xunit, :);
        end
        
        for xtarg = 1:length(xtarget)
            
            fprintf(xout_txt, '#### target: %s\n', cond_names{xtarg});
            
            for gg = 1:length(xtarget)
                xx = xconfusion{xtarg}(:, findCol(xheader, xguess{gg}));
                fprintf(xout_txt, 'guess %s: M=%4.4f, SE=%4.4f\n', cond_names{gg}, ...
                    mean(xx.Variables), std(xx.Variables)/sqrt(length(xsub_groups)));
            end
            
            fprintf(xout_txt,'\n');
            
            clear xpvalue xstats
            it_targs = targ_array(~ismember(targ_array, xtarg));
            for it_y = 1:length(it_targs)
                ytarg = it_targs(it_y);
                
                %*************** ttest
                xx = xconfusion{xtarg}(:, findCol(xheader, xguess{xtarg}));
                yy = xconfusion{xtarg}(:, findCol(xheader, xguess{ytarg}));
                
                [~, xpvalue(it_y), ~, xstats] = ttest(xx.Variables, yy.Variables, 'Alpha', 0.05);
               
                fprintf(xout_txt, '%s vs. %s: T(%s)=%4.4f, P=%4.4f\n', ...
                    cond_names{xtarg}, cond_names{ytarg}, num2str(xstats.df), xstats.tstat, xpvalue(it_y));
                
            end
            
            fprintf(xout_txt, '\n\n********* FDR p-value\n');
            
            [~, ~, ~, xadj_p] = fdr_bh(xpvalue, 0.05, 'pdep', 'no');
            
            for it_y = 1:length(it_targs)
                ytarg = it_targs(it_y);
                fprintf(xout_txt, '%s vs. %s: P=%4.4f\n', ...
                    cond_names{xtarg}, cond_names{ytarg}, xadj_p(it_y));
            end
            
            fprintf(xout_txt, '\n');
            
        end
    end
    
    %% ============= COMPARING: REPLACES: ACCURACY/EVIDENCE
    targ_names = {'target','nontarget'};
    if ~(args.four_oper_regress)
        clear xconfusion
        xmeas_names = {'accuracy','evidence'};
        targ_array  = 1:n_condition;
        
        fprintf(xout_txt, '\n\n####################################################################\n');
        fprintf(xout_txt, '#### Replaces: paired T test\n');
        fprintf(xout_txt, '####################################################################\n\n');
        
        for xmeas = 1:length(xmeas_names)%1_auc, 2_auc_baseline, 3_4_classifier accuracy, 4_classifier evidence
            clear xtable xheader xtarget xguess
            
            fprintf(xout_txt, '\n=========================================\n');
            fprintf(xout_txt, '======== %s\n\n', xmeas_names{xmeas});
            
            %*************** data matrix
            xcsv = fullfile(dirs.mvpa.group.out{xph}, ...
                sprintf('grp_%s_%s_n%s.csv', xmeas_names{xmeas}, basename, num2str(n_subs)));
            
            xtable  = readtable(xcsv);
            xheader = xtable.Properties.VariableNames;
            
            xtarget = unique(xtable.Row)';
            xguess  = xtable.Properties.VariableNames(~ismember(xtable.Properties.VariableNames, 'Row'));
            
            for xtarg = 1:length(xtarget)
                xunit = ismember(xtable.Row, xtarget{xtarg});
                xconfusion{xtarg} = xtable(xunit, :);
            end
            
            %*************** 1_target (actual), 2_guessed (predicted)
            xx1 = xconfusion{2}(:, findCol(xheader, {xguess{2}, xguess{3}}));
            xx2 = xconfusion{3}(:, findCol(xheader, {xguess{3}, xguess{2}}));
            
            fprintf(xout_txt, '#### replaces: target vs. nontarg\n');
            
            xx = vertcat(xx1, xx2);
            
            for xtarg = 1:2
                yy = xx.Variables;
                
                fprintf(xout_txt, '%s: M=%4.4f, SE=%4.4f\n', targ_names{xtarg}, ...
                    mean(yy(:,xtarg)), std(yy(:,xtarg))/sqrt(length(xsub_groups)));
            end
            
            fprintf(xout_txt,'\n');
            
            clear xpvalue xstats
            
            [~, xpvalue, ~, xstats] = ttest(yy(:,1), yy(:,2), 'Alpha', 0.05);
            
            fprintf(xout_txt, 'targ vs. nontarg: T(%s)=%4.4f, P=%4.4f\n', ...
                num2str(xstats.df), xstats.tstat, xpvalue);
            
        end
    else
        %% ============= COMPARING CONDITIONS: suppress(3) vs. clear(4)
        clear xconfusion
        xmeas_names = {'accuracy','evidence'};
        targ_array  = 1:n_condition;
        
        fprintf(xout_txt, '\n\n####################################################################\n');
        fprintf(xout_txt, '#### Suppress vs. Clear: paired T test\n');
        fprintf(xout_txt, '####################################################################\n\n');
        
        for xmeas = 1:length(xmeas_names)%1_auc, 2_auc_baseline, 3_4_classifier accuracy, 4_classifier evidence
            clear xtable xheader xtarget xguess
            
            fprintf(xout_txt, '\n=========================================\n');
            fprintf(xout_txt, '======== %s\n\n', xmeas_names{xmeas});
            
            %*************** data matrix
            xcsv = fullfile(dirs.mvpa.group.out{xph}, ...
                sprintf('grp_%s_%s_n%s.csv', xmeas_names{xmeas}, basename, num2str(n_subs)));
            
            xtable  = readtable(xcsv);
            xheader = xtable.Properties.VariableNames;
            
            xtarget = unique(xtable.Row)';
            xguess  = xtable.Properties.VariableNames(~ismember(xtable.Properties.VariableNames, 'Row'));
            
            for xtarg = 1:length(xtarget)
                xunit = ismember(xtable.Row, xtarget{xtarg});
                xconfusion{xtarg} = xtable(xunit, :);
            end
            
            %*************** 1_target (actual), 2_guessed (predicted)
            % suppress(3) vs. clear(4)
            t_xx1 = xconfusion{3}(:, findCol(xheader, {xguess{3}, xguess{4}}));
            t_xx2 = xconfusion{4}(:, findCol(xheader, {xguess{4}, xguess{3}}));
            
            xx1 = t_xx1.Variables;
            xx2 = t_xx2.Variables;
            
            fprintf(xout_txt, '#### suppress vs. clear: target vs. nontarg\n');
            
            xx = vertcat(xx1, xx2);
            
            for xtarg = 1:2
                fprintf(xout_txt, '%s: M=%4.4f, SE=%4.4f\n', targ_names{xtarg}, ...
                    mean(xx(:,xtarg)), std(xx(:,xtarg))/sqrt(length(xsub_groups)));
            end
            
            fprintf(xout_txt,'\n');
            
            clear xpvalue xstats
            
            [~, xpvalue, ~, xstats] = ttest(xx(:,1), xx(:,2), 'Alpha', 0.05);
            
            fprintf(xout_txt, 'targ vs. nontarg: T(%s)=%4.4f, P=%4.4f\n', ...
                num2str(xstats.df), xstats.tstat, xpvalue);
            
        end
        
        %% ============= COMPARING CONDITIONS: replace(2) vs. two_removals (3,4)
        clear xconfusion
        xmeas_names = {'accuracy','evidence'};
        targ_array  = 1:n_condition;
        
        fprintf(xout_txt, '\n\n####################################################################\n');
        fprintf(xout_txt, '#### Replace vs. Removals: paired T test\n');
        fprintf(xout_txt, '####################################################################\n\n');
        
        for xmeas = 1:length(xmeas_names)%1_auc, 2_auc_baseline, 3_4_classifier accuracy, 4_classifier evidence
            clear xtable xheader xtarget xguess xx1 xx2
            
            fprintf(xout_txt, '\n=========================================\n');
            fprintf(xout_txt, '======== %s\n\n', xmeas_names{xmeas});
            
            %*************** data matrix
            xcsv = fullfile(dirs.mvpa.group.out{xph}, ...
                sprintf('grp_%s_%s_n%s.csv', xmeas_names{xmeas}, basename, num2str(n_subs)));
            
            xtable  = readtable(xcsv);
            xheader = xtable.Properties.VariableNames;
            
            xtarget = unique(xtable.Row)';
            xguess  = xtable.Properties.VariableNames(~ismember(xtable.Properties.VariableNames, 'Row'));
            
            for xtarg = 1:length(xtarget)
                xunit = ismember(xtable.Row, xtarget{xtarg});
                xconfusion{xtarg} = xtable(xunit, :);
            end
            
            %*************** 1_target (actual), 2_guessed (predicted)
            % replace(2) vs. two_removals (3,4)
            t_xx1 = xconfusion{2}(:, findCol(xheader, {xguess{2}, xguess{3}, xguess{4}}));
            t_xx2 = xconfusion{3}(:, findCol(xheader, {xguess{3}, xguess{2}}));
            t_xx3 = xconfusion{4}(:, findCol(xheader, {xguess{4}, xguess{2}}));
            
            xx1 = t_xx1.Variables; xx1 = [xx1(:,1) mean(xx1(:,2:3),2)];
            xx2 = t_xx2.Variables;
            xx3 = t_xx2.Variables;
            
            xx = vertcat(xx1, xx2, xx3);
            
            fprintf(xout_txt, '#### suppress vs. clear: target vs. nontarg\n');
            
            for xtarg = 1:2
                
                fprintf(xout_txt, '%s: M=%4.4f, SE=%4.4f\n', targ_names{xtarg}, ...
                    mean(xx(:,xtarg)), std(xx(:,xtarg))/sqrt(length(xsub_groups)));
            end
            
            fprintf(xout_txt,'\n');
            
            clear xpvalue xstats
            
            [~, xpvalue, ~, xstats] = ttest(xx(:,1), xx(:,2), 'Alpha', 0.05);
            
            fprintf(xout_txt, 'targ vs. nontarg: T(%s)=%4.4f, P=%4.4f\n', ...
                num2str(xstats.df), xstats.tstat, xpvalue);
            
        end
        
    end
    
%     fclose(xout_txt);
    
end

end