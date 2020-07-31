function[] = grp_analysis_timecourse_operation_zscored(args, grp_mvpaout, dirs)

%---------------------------------------------------------------------------
%*************** group: timecourse of operation during study
%---------------------------------------------------------------------------

xph = args.xphase;
fprintf('\n(+) timecourse of working memory: %s\n', args.phase_name{xph});

xsub_groups = args.filtered_subs;%args.g_sub

%% ============= UNPACK PARAMETERS
xindex           = args.index{xph};% param index from study

xparam           = xindex.param;
xcond_name       = xparam.conds_names;
n_condition      = length(xcond_name);
dur_stim         = xparam.dur_stim;
dur_manipulation = xparam.dur_manipulation;% 6 TR
it_trs           = 1:xparam.n_tc_trs;%timecourse trs

xalpha           = 0.05;
condition_names  = {'maintain','replace category',...
    'replace subcategory','target suppress','global clear'};

%*************** BASELINE CORRECTION
baseline_trs    = args.baseline_trs;

%*************** output basename
base_name       = sprintf('zscored_%s_n%s', args.analysis_basename, num2str(length(xsub_groups)));
xout_dir        = fullfile(dirs.mvpa.group.out{xph},'zscored');
if ~isdir(xout_dir), mkdir(xout_dir); end

%% ============= PLOT PARAMS
xcolor{1}  = [238, 20, 91]/255;% targ
xcolor{2}  = [144, 144, 144]/255;% non_target

x_tick     = it_trs;
x_lim      = [x_tick(1) x_tick(end)];% + [-0.5 0.5];
y_lim      = [-2 3];
n_ticks    = length(x_tick);
line_w     = 2;
line_w_err = 0.1;

fig_rect   = [0 0 1200 800]; 

%*************** separated operations
xcond_color  = args.cond_color;
% xcate_color  = args.cate_color;
% 
% xbase_color  = args.base_color;
% xonset_color = args.onset_color;

%% ============= ZSCORING: MEAN & SD OF EVIDENCE WITHIN CONDITION
n_targ = 2;%targ, nontarg

for xcond = 1:n_condition
    xevidence = [];
    for xtr = it_trs
        for xtarg = 1:n_targ%1_targ, 2_nontarg
            
            %*************** fixed effect
            for xsub = xsub_groups
                xevidence = horzcat(xevidence, ...
                    mean(grp_mvpaout{xsub}.decode.timecourse.operation{xcond}.target{xtarg}.evidence{xtr})); %#ok<*AGROW>
            end
        end
    end
    
    mean_evidence{xcond} = mean(xevidence);
    sd_evidence{xcond}   = std(xevidence);
end

%% ============= TIMECOURSE STRUCTURE
%*************** working memory operation
%*************** timecourse_fixed.mean{xcond}(1,:): target
%*************** timecourse_fixed.mean{xcond}(2,:): nontarget baseline
n_targ = 2;%targ, nontarg

for xcond = 1:n_condition
    for xtr = it_trs
        for xtarg = 1:n_targ%1_targ, 2_nontarg
            xevidence_random = [];
            
            %*************** random effect
            for xsub = xsub_groups
                xevidence_random = horzcat(xevidence_random, ...
                    mean(grp_mvpaout{xsub}.decode.timecourse.operation{xcond}.target{xtarg}.evidence{xtr}));
            end
            
            z_evidence = (xevidence_random - mean_evidence{xcond})/sd_evidence{xcond};
            
            timecourse_random.evidence{xcond}{xtarg}{xtr} = z_evidence;
            timecourse_random.mean{xcond}(xtarg, xtr)     = mean(z_evidence);
            timecourse_random.se{xcond}(xtarg, xtr)       = std(z_evidence)/sqrt(length(xsub_groups));
        end
    end
end

%% ============= STATS
for xcond = 1:n_condition
    for xtr = it_trs

        timecourse_random.t_test{xcond}{xtr}.pvalue = nan(n_targ, n_targ);
        
        for xcol = 1:(n_targ-1)
            for xrow = (xcol+1):n_targ
                xtarg_col = timecourse_random.evidence{xcond}{xcol}{xtr};%array_1
                xtarg_row = timecourse_random.evidence{xcond}{xrow}{xtr};%array_2
                
                %*************** ttest
                [xhypothesis, xpvalue, ~, xstats] = ttest(xtarg_row, xtarg_col, 'Alpha', 0.05);
                
                timecourse_random.t_test{xcond}{xtr}.pvalue(xrow, xcol)     = xpvalue;
                timecourse_random.t_test{xcond}{xtr}.hypothesis(xrow, xcol) = xhypothesis;
                timecourse_random.t_test{xcond}{xtr}.stats{xrow, xcol}      = xstats;
            end
        end
    end
end

%% ============= ZSCORING: MEAN & SD OF EVIDENCE WITHIN CONDITION
clear mean_evidence sd_evidence
n_targ = n_condition;

for xcond = 1:n_condition
    xevidence = [];
    for xtr = it_trs
        for xtarg = 1:n_targ%1_targ, 2_nontarg
            
            %*************** fixed effect
            for xsub = xsub_groups
                xevidence = horzcat(xevidence, ...
                    mean(grp_mvpaout{xsub}.decode.timecourse.operation{xcond}.decoded_operation{xtarg}.evidence{xtr}));
            end
        end
    end
    
    mean_evidence{xcond} = mean(xevidence);
    sd_evidence{xcond}   = std(xevidence);
end

%% ============= TIMECOURSE STRUCTURE: ALL TARGET OPERATIONS
%*************** separated working memory operation
%*************** decode.timecourse.operation{xcond}.decoded_operation{xoper}.evidence{xtr}
n_targ = n_condition;% all operations

for xcond = 1:n_condition
    for xtr = it_trs
        for xtarg = 1:n_targ
            xevidence_random = [];
            
            %*************** random effect
            for xsub = xsub_groups
                xevidence_random = horzcat(xevidence_random, ...
                    mean(grp_mvpaout{xsub}.decode.timecourse.operation{xcond}.decoded_operation{xtarg}.evidence{xtr}));
            end
                    
            z_evidence = (xevidence_random - mean_evidence{xcond})/sd_evidence{xcond};

            all_timecourse_random.evidence{xcond}{xtarg}{xtr} = z_evidence;
            all_timecourse_random.mean{xcond}(xtarg, xtr)     = mean(z_evidence);
            all_timecourse_random.se{xcond}(xtarg, xtr)       = std(z_evidence)/sqrt(length(xsub_groups));
        end
    end
end

%% ============= STATS
for xcond = 1:n_condition
    for xtr = it_trs

        all_timecourse_random.t_test{xcond}{xtr}.pvalue = nan(n_targ-1, n_targ-1);
                
        for xcol = 1:(n_targ-1)
            for xrow = (xcol+1):n_targ
                xtarg_col = all_timecourse_random.evidence{xcond}{xcol}{xtr};%array_1
                xtarg_row = all_timecourse_random.evidence{xcond}{xrow}{xtr};%array_2
                
                %*************** ttest
                [xhypothesis, xpvalue, ~, xstats] = ttest(xtarg_row, xtarg_col, 'Alpha', 0.05);
                
                all_timecourse_random.t_test{xcond}{xtr}.pvalue(xrow, xcol)     = xpvalue;
                all_timecourse_random.t_test{xcond}{xtr}.hypothesis(xrow, xcol) = xhypothesis;
                all_timecourse_random.t_test{xcond}{xtr}.stats{xrow, xcol}      = xstats;
            end
        end
    end
end

%% ============= TIMECOURSE STRUCTURE: ALL TARGET OPERATIONS
%*************** separated working memory operation
%*************** decode.timecourse.operation{xcond}.decoded_operation{xoper}.evidence{xtr}

for xcond = 1:n_condition
    
    for xtr = it_trs
        xevidence_random = all_timecourse_random.evidence{xcond}{xcond}{xtr};
        
        target_timecourse.evidence{xcond}{xtr} = xevidence_random;
        target_timecourse.mean{xcond}(xtr)     = mean(xevidence_random);
        target_timecourse.se{xcond}(xtr)       = std(xevidence_random)/sqrt(length(xsub_groups));
    end
end

%% ============= STATS
for xtr = it_trs
    
    target_timecourse.t_test{xtr}.pvalue = nan(n_condition, n_condition);
    
    for xcol = 1:(n_condition-1)
        for xrow = (xcol+1):n_condition
            xtarg_col = target_timecourse.evidence{xcol}{xtr};%array_1
            xtarg_row = target_timecourse.evidence{xrow}{xtr};%array_2
            
            %*************** ttest
            [xhypothesis, xpvalue, ~, xstats] = ttest(xtarg_row, xtarg_col, 'Alpha', 0.05);
            
            target_timecourse.t_test{xtr}.pvalue(xrow, xcol)     = xpvalue;
            target_timecourse.t_test{xtr}.hypothesis(xrow, xcol) = xhypothesis;
            target_timecourse.t_test{xtr}.stats{xrow, xcol}      = xstats;
        end
    end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%============= BASELINE CORRECTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ============= TIMECOURSE STRUCTURE: ALL OPERATIONS
%*************** separated working memory operation
%*************** decode.timecourse.operation{xcond}.decoded_operation{xoper}.evidence{xtr}
n_targ = n_condition;% all operations

for xcond = 1:n_condition
    
    %*************** baseline_trs: 6TRs
    xbase_evidence = zeros(n_targ, baseline_trs);
    for xtarg = 1:n_targ
        for xtr = 1:baseline_trs
            it_tr = xtr + dur_stim;
            xbase_evidence(xtarg, xtr) = mean(all_timecourse_random.evidence{xcond}{xtarg}{it_tr});
        end
        
        mean_baseline(xtarg) = mean(xbase_evidence(xtarg,:));
    end
    
    %*************** corrected evidence
    for xtr = it_trs
        for xtarg = 1:n_targ
            xcorrect_evidence = ...
                all_timecourse_random.evidence{xcond}{xtarg}{xtr} - mean_baseline(xtarg);
            
            all_timecourse_random.baseline.evidence{xcond}{xtarg}{xtr} = xcorrect_evidence;
            all_timecourse_random.baseline.mean{xcond}(xtarg, xtr)     = mean(xcorrect_evidence);
            all_timecourse_random.baseline.se{xcond}(xtarg, xtr)       = std(xcorrect_evidence)/sqrt(length(xsub_groups));
        end
    end
end

%% ============= STATS
for xcond = 1:n_condition
    for xtr = it_trs

        all_timecourse_random.baseline.t_test{xcond}{xtr}.pvalue = nan(n_targ-1, n_targ-1);
                
        for xcol = 1:(n_targ-1)
            for xrow = (xcol+1):n_targ
                xtarg_col = all_timecourse_random.baseline.evidence{xcond}{xcol}{xtr};%array_1
                xtarg_row = all_timecourse_random.baseline.evidence{xcond}{xrow}{xtr};%array_2
                
                %*************** ttest
                [xhypothesis, xpvalue, ~, xstats] = ttest(xtarg_row, xtarg_col, 'Alpha', 0.05);
                
                all_timecourse_random.baseline.t_test{xcond}{xtr}.pvalue(xrow, xcol)     = xpvalue;
                all_timecourse_random.baseline.t_test{xcond}{xtr}.hypothesis(xrow, xcol) = xhypothesis;
                all_timecourse_random.baseline.t_test{xcond}{xtr}.stats{xrow, xcol}      = xstats;
            end
        end
    end
end

%% ============= TIMECOURSE STRUCTURE: ALL TARGET OPERATIONS
%*************** separated working memory operation
%*************** decode.timecourse.operation{xcond}.decoded_operation{xoper}.evidence{xtr}

%*************** baseline_trs: 6TRs
xbase_evidence = zeros(n_condition, baseline_trs);
for xcond = 1:n_condition 
    for xtr = 1:baseline_trs
        it_tr = xtr + dur_stim;
        xbase_evidence(xcond, xtr) = mean(target_timecourse.evidence{xcond}{it_tr});
    end
    
    mean_baseline(xcond) = mean(xbase_evidence(xcond,:));
end

%*************** corrected evidence
for xcond = 1:n_condition
    for xtr = it_trs
        xcorrect_evidence = ...
            target_timecourse.evidence{xcond}{xtr} - mean_baseline(xcond);
        
        target_timecourse.baseline.evidence{xcond}{xtr} = xcorrect_evidence;
        target_timecourse.baseline.mean{xcond}(xtr)     = mean(xcorrect_evidence);
        target_timecourse.baseline.se{xcond}(xtr)       = std(xcorrect_evidence)/sqrt(length(xsub_groups));
    end
end

%% ============= STATS
for xtr = it_trs
    
    target_timecourse.baseline.t_test{xtr}.pvalue = nan(n_condition, n_condition);
    
    for xcol = 1:(n_condition-1)
        for xrow = (xcol+1):n_condition
            xtarg_col = target_timecourse.baseline.evidence{xcol}{xtr};%array_1
            xtarg_row = target_timecourse.baseline.evidence{xrow}{xtr};%array_2
            
            %*************** ttest
            [xhypothesis, xpvalue, ~, xstats] = ttest(xtarg_row, xtarg_col, 'Alpha', 0.05);
            
            target_timecourse.baseline.t_test{xtr}.pvalue(xrow, xcol)     = xpvalue;
            target_timecourse.baseline.t_test{xtr}.hypothesis(xrow, xcol) = xhypothesis;
            target_timecourse.baseline.t_test{xtr}.stats{xrow, xcol}      = xstats;
        end
    end
end

%% ============= SAVE
z_grp_timecourse.random     = timecourse_random; 
z_grp_timecourse.all.random = all_timecourse_random;
z_grp_timecourse.target     = target_timecourse;%#ok<*STRNU> %random effect

fname = fullfile(xout_dir, sprintf('grp_operation_timecourse_%s.mat', base_name));
save(fname, 'z_grp_timecourse');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%============= PLOTTING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%============= CONTINUOUS PLOT: RANDOM EFFECT
n_targ = 2;
xgap   = 0.1;

for xcond = 1:n_condition
    clear xmean xstd fity
    
    LogRegs_fig = figure;
    set(LogRegs_fig, 'Position', fig_rect)
    
    %*************** mean line plots
    fitx = linspace(1, n_ticks, n_ticks*10);
    
    for xtarg = 1:n_targ
        xmean{xtarg} = timecourse_random.mean{xcond}(xtarg, x_tick);
        xstd{xtarg}  = timecourse_random.se{xcond}(xtarg, x_tick); %#ok<*AGROW>
        
        fity{xtarg}  = interp1(x_tick, xmean{xtarg}, fitx,'spline');
        
        plot(fitx, fity{xtarg}, '-','Color', xcolor{xtarg}, 'LineWidth', line_w); hold on;
    end
    
    %*************** legend
    xlegend = {'target','baseline'};
    lg          = legend(xlegend);
    lg.Location = 'SouthEast';
    legend(xlegend,'AutoUpdate','off')
    grid on
    
    %*************** stats stars
    for xcol = 1:(n_targ-1)
        for xrow = (xcol+1):n_targ
            xheight = y_lim(2) - xgap;
            
            for xtr = x_tick
                xpvalue = timecourse_random.t_test{xcond}{xtr}.pvalue(xrow, xcol);
                
                if xpvalue <= xalpha
                    text(xtr-0.1, xheight, '*', 'FontSize', 20);
                    text(xtr-0.1, xheight-xgap, sprintf('%1.2f', xpvalue), 'FontSize', 10);
                end
                
                plot([xtr xtr], y_lim, ':', 'Color', [0.75 0.75 0.75])
                
            end
        end
    end
    
    %*************** std error-bar filling
    for xtarg = 1:n_targ
        clear xerr
        xerr(1,:) = xmean{xtarg} - xstd{xtarg};
        xerr(2,:) = xmean{xtarg} + xstd{xtarg};
        
        for i = 1:2, fit_err(:,i) = interp1(x_tick, xerr(i,:), fitx,'spline'); end
        
        in_between = [fit_err(:,1)', fliplr(fit_err(:,2)')];
        h = fill([fitx fliplr(fitx)], in_between, xcolor{xtarg});
        set(h,'facealpha', .2)
        
        for i=1:2, plot(fitx, fit_err(:,i), '-', 'Color', xcolor{xtarg},'LineWidth',line_w_err); end
    end
     
    set(gca,'xlim', x_lim, 'ylim', y_lim);
    set(gca,'XTick', x_tick, 'XTickLabel', x_tick-1)
    set(gca,'YTick', y_lim(1):0.5:y_lim(end))
    
    title(sprintf('Zscored Timecourse %s (N=%s)', condition_names{xcond}, num2str(length(xsub_groups))));
    xlabel('Volume (tr)');
    ylabel('classifier evidence');
    
    %*************** onset lines
    plot([dur_stim dur_stim] + 1, y_lim,'-', 'Color', 'r')
    plot([dur_stim+dur_manipulation, dur_stim+dur_manipulation] + 1, y_lim,'-', 'Color', 'r')
    
    text(1, y_lim(1) + xgap, 'stim onset', 'Color','r');
    text(dur_stim + 1, y_lim(1) + xgap, 'operation onset','Color','r');
    text(dur_stim + dur_manipulation + 1, y_lim(1) + xgap, 'fixation onset','Color','r');
    
    %% ============= SAVE FIGURE
    fig_fname = fullfile(xout_dir, sprintf('grp_plot_continous_operation_timecourse_random_%s_%s', ...
        xcond_name{xcond}, base_name));
    
    savefig(LogRegs_fig, sprintf('%s.fig', fig_fname));
    set(gcf,'PaperUnits','inches','PaperPosition',fig_rect/100)
    print('-djpeg', sprintf('%s.jpg',fig_fname), '-r100')
    saveas(LogRegs_fig, sprintf('%s.jpg',fig_fname), 'jpg')
    
    close(LogRegs_fig);
    
end

%% ============= ALL OPERATIONS: CONTINUOUS PLOT: RANDOM EFFECT
xgap  = 0.1;
y_lim = [-2 3];

for xfill = 0:1
    std_fill   = xfill;
    n_targ     = n_condition;%all operations
    xsel_conds = 1:5;
    
    for xcond = 1:n_condition
        clear xmean xstd fity
        
        LogRegs_fig = figure;
        set(LogRegs_fig, 'Position', fig_rect)
        
        %*************** mean line plots
        fitx = linspace(1, n_ticks, n_ticks*10);
        
        for xtarg = 1:n_targ
            xmean{xtarg} = all_timecourse_random.mean{xcond}(xtarg, x_tick);
            xstd{xtarg}  = all_timecourse_random.se{xcond}(xtarg, x_tick); %#ok<*AGROW>
            
            fity{xtarg}  = interp1(x_tick, xmean{xtarg}, fitx,'spline');
            
            if xtarg==xcond
                plot(fitx, fity{xtarg}, '-', 'Color', xcond_color{xtarg}, 'LineWidth', line_w + 2); hold on;
            else
                plot(fitx, fity{xtarg}, '-', 'Color', xcond_color{xtarg}, 'LineWidth', line_w); hold on;
            end
        end
        
        %*************** legend
        xlegend     = condition_names;
        lg          = legend(xlegend);
        lg.Location = 'SouthEast';
        legend(xlegend,'AutoUpdate','off')
        grid on
        
        %*************** stats stars
        n = 0;
        
        for i = 1:(n_targ-1)
            xcol = xsel_conds(i);
            for j = (i+1):n_targ
                xrow = xsel_conds(j);
                xheight = (y_lim(2) - xgap) - (n * xgap);
                text(-2, xheight + xgap, sprintf('%s vs. %s', num2str(xcol), num2str(xrow)));
                plot(x_lim, [xheight xheight], ':', 'Color', [0.75 0.75 0.75])
                
                for xtr = x_tick
                    xpvalue = all_timecourse_random.t_test{xcond}{xtr}.pvalue(xrow, xcol);
                    
                    if xpvalue <= xalpha
                        text(xtr-0.1, xheight, '*', 'FontSize', 20);
                    end
                    
                    plot([xtr xtr], y_lim, ':', 'Color', [0.75 0.75 0.75])
                    
                end
                n = n + 1;
            end
        end
        
        %*************** std error-bar filling
        if std_fill
            for xtarg = 1:n_targ
                clear xerr
                xerr(1,:) = xmean{xtarg} - xstd{xtarg};
                xerr(2,:) = xmean{xtarg} + xstd{xtarg};
                
                for i = 1:2, fit_err(:,i) = interp1(x_tick, xerr(i,:), fitx,'spline'); end
                
                in_between = [fit_err(:,1)', fliplr(fit_err(:,2)')];
                h = fill([fitx fliplr(fitx)], in_between, xcond_color{xtarg});
                set(h,'facealpha', .2)
                
                for i=1:2, plot(fitx, fit_err(:,i), '-', 'Color', xcond_color{xtarg},'LineWidth',line_w_err); end
            end
        end
        
        set(gca,'xlim', x_lim, 'ylim', y_lim);
        set(gca,'XTick', x_tick, 'XTickLabel', x_tick-1)
        set(gca,'YTick', y_lim(1):0.5:y_lim(end))
        
        title(sprintf('Timecourse %s (N=%s)', condition_names{xcond}, num2str(length(xsub_groups))));
        xlabel('Volume (tr)');
        ylabel('classifier evidence');
        
        %*************** onset lines
        plot([dur_stim dur_stim] + 1, y_lim,'-', 'Color', 'r')
        plot([dur_stim+dur_manipulation, dur_stim+dur_manipulation] + 1, y_lim,'-', 'Color', 'r')
        
        text(1, y_lim(1) + xgap, 'stim onset', 'Color','r');
        text(dur_stim + 1, y_lim(1) + xgap, 'operation onset','Color','r');
        text(dur_stim + dur_manipulation + 1, y_lim(1) + xgap, 'fixation onset','Color','r');
        
        %% ============= SAVE FIGURE
        if std_fill
            fig_fname = fullfile(xout_dir, ...
                sprintf('grp_plot_continous_operation_all_timecourse_random_%s_%s', ...
                xcond_name{xcond}, base_name));
        else
            fig_fname = fullfile(xout_dir, ...
                sprintf('grp_plot_continous_operation_all_timecourse_line_random_%s_%s', ...
                xcond_name{xcond}, base_name));
        end
        
        savefig(LogRegs_fig, sprintf('%s.fig', fig_fname));
        set(gcf,'PaperUnits','inches','PaperPosition',fig_rect/100)
        print('-djpeg', sprintf('%s.jpg',fig_fname), '-r100')
        saveas(LogRegs_fig, sprintf('%s.jpg',fig_fname), 'jpg')
        
        close(LogRegs_fig);
        
    end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%============= BASELINE CORRECTED
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ============= ALL OPERATIONS: CONTINUOUS PLOT: RANDOM EFFECT
for xfill = 0:1
    std_fill   = xfill;
    y_lim_corr = [-2.5 3];%y_lim;% - 0.5;
    
    n_targ     = n_condition;%all operations    
    xsel_conds = 1:5;
    
    for xcond = 1:n_condition
        clear xmean xstd fity
        
        LogRegs_fig = figure;
        set(LogRegs_fig, 'Position', fig_rect)
        
        %*************** mean line plots
        fitx = linspace(1, n_ticks, n_ticks*10);
        
        for xtarg = 1:n_targ
            xmean{xtarg} = all_timecourse_random.baseline.mean{xcond}(xtarg, x_tick);
            xstd{xtarg}  = all_timecourse_random.baseline.se{xcond}(xtarg, x_tick); %#ok<*AGROW>
            
            fity{xtarg}  = interp1(x_tick, xmean{xtarg}, fitx,'spline');
            
            if xtarg==xcond
                plot(fitx, fity{xtarg}, '-', 'Color', xcond_color{xtarg}, 'LineWidth', line_w + 2); hold on;
            else
                plot(fitx, fity{xtarg}, '-', 'Color', xcond_color{xtarg}, 'LineWidth', line_w); hold on;
            end
        end
        
        %*************** legend
        xlegend     = condition_names;
        lg          = legend(xlegend);
        lg.Location = 'SouthEast';
        legend(xlegend,'AutoUpdate','off')
        grid on
        
        %*************** stats stars
        n = 0;
        
        for i = 1:(n_targ-1)
            xcol = xsel_conds(i);
            for j = (i+1):n_targ
                xrow = xsel_conds(j);
                xheight = (y_lim_corr(2)-xgap) - (n * xgap);
                text(-2, xheight + 0.02, sprintf('%s vs. %s', num2str(xcol), num2str(xrow)));
                plot(x_lim, [xheight xheight], ':', 'Color', [0.75 0.75 0.75])
                
                for xtr = x_tick
                    xpvalue = all_timecourse_random.baseline.t_test{xcond}{xtr}.pvalue(xrow, xcol);
                    
                    if xpvalue <= xalpha
                        text(xtr-0.1, xheight, '*', 'FontSize', 20);
                    end
                    
                    plot([xtr xtr], y_lim_corr, ':', 'Color', [0.75 0.75 0.75])
                    
                end
                n = n + 1;
            end
        end
        
        %*************** std error-bar filling
        if std_fill
            for xtarg = 1:n_targ
                clear xerr
                xerr(1,:) = xmean{xtarg} - xstd{xtarg};
                xerr(2,:) = xmean{xtarg} + xstd{xtarg};
                
                for i = 1:2, fit_err(:,i) = interp1(x_tick, xerr(i,:), fitx,'spline'); end
                
                in_between = [fit_err(:,1)', fliplr(fit_err(:,2)')];
                h = fill([fitx fliplr(fitx)], in_between, xcond_color{xtarg});
                set(h,'facealpha', .2)
                
                for i=1:2, plot(fitx, fit_err(:,i), '-', 'Color', xcond_color{xtarg},'LineWidth',line_w_err); end
            end
        end
        
        set(gca,'xlim', x_lim, 'ylim', y_lim_corr);
        set(gca,'XTick', x_tick, 'XTickLabel', x_tick-1)
        set(gca,'YTick', y_lim_corr(1):0.5:y_lim_corr(end))
        
        title(sprintf('Zscored Timecourse %s (N=%s)', condition_names{xcond}, num2str(length(xsub_groups))));
        xlabel('Volume (tr)');
        ylabel('classifier evidence');
        
        %*************** onset lines
        plot(x_lim, [0 0], ':', 'Color', 'r', 'LineWidth', 1.5)%baseline
        plot([dur_stim dur_stim] + 1, y_lim_corr,'-', 'Color', 'r')
        plot([dur_stim+dur_manipulation, dur_stim+dur_manipulation] + 1, y_lim_corr,'-', 'Color', 'r')
        
        text(1, y_lim_corr(1) + xgap, 'stim onset', 'Color','r');
        text(dur_stim + 1, y_lim_corr(1) + xgap, 'operation onset','Color','r');
        text(dur_stim + dur_manipulation + 1, y_lim_corr(1) + xgap, 'fixation onset','Color','r');
        
        %% ============= SAVE FIGURE
        if std_fill
            fig_fname = fullfile(xout_dir, ...
                sprintf('grp_plot_corrected_continous_operation_all_timecourse_random_%s_%s', ...
                xcond_name{xcond}, base_name));
        else
            fig_fname = fullfile(xout_dir, ...
                sprintf('grp_plot_corrected_continous_operation_all_timecourse_line_random_%s_%s', ...
                xcond_name{xcond}, base_name));
        end
        
        savefig(LogRegs_fig, sprintf('%s.fig', fig_fname));
        set(gcf,'PaperUnits','inches','PaperPosition',fig_rect/100)
        print('-djpeg', sprintf('%s.jpg',fig_fname), '-r100')
        saveas(LogRegs_fig, sprintf('%s.jpg',fig_fname), 'jpg')
        
        close(LogRegs_fig);
        
    end
end

%% ============= TARGET OPERATIONS: CONTINUOUS PLOT: RANDOM EFFECT
for xfill = 0:1
    std_fill   = xfill;
%     y_lim_corr = y_lim - 0.5;
    xsel_conds = 1:5;
    
    LogRegs_fig = figure;
    set(LogRegs_fig, 'Position', fig_rect)
    
    %*************** mean line plots
    fitx = linspace(1, n_ticks, n_ticks*10);
    
    clear xmean xstd
    for xtarg = 1:n_condition
        xmean{xtarg} = target_timecourse.baseline.mean{xtarg}(x_tick);
        xstd{xtarg}  = target_timecourse.baseline.se{xtarg}(x_tick); %#ok<*AGROW>
        
        fity{xtarg}  = interp1(x_tick, xmean{xtarg}, fitx,'spline');
        
        plot(fitx, fity{xtarg}, '-', 'Color', xcond_color{xtarg}, 'LineWidth', line_w); hold on;
    end
    
    %*************** legend
    xlegend     = condition_names;
    lg          = legend(xlegend);
    lg.Location = 'SouthEast';
    legend(xlegend,'AutoUpdate','off')
    grid on
    
    %*************** stats stars
    n = 0;
    
    for i = 1:(n_condition-1)
        xcol = xsel_conds(i);
        for j = (i+1):n_condition
            xrow = xsel_conds(j);
            xheight = (y_lim_corr(2) - xgap) - (n * xgap);
            text(-2, xheight + 0.02, sprintf('%s vs. %s', num2str(xcol), num2str(xrow)));
            plot(x_lim, [xheight xheight], ':', 'Color', [0.75 0.75 0.75])
            
            for xtr = x_tick
                xpvalue = target_timecourse.baseline.t_test{xtr}.pvalue(xrow, xcol);
                
                if xpvalue <= xalpha
                    text(xtr-0.1, xheight, '*', 'FontSize', 20);
                end
                
                plot([xtr xtr], y_lim_corr, ':', 'Color', [0.75 0.75 0.75])
                
            end
            n = n + 1;
        end
    end
    
    %*************** std error-bar filling
    if std_fill
        for xtarg = 1:n_condition
            clear xerr
            xerr(1,:) = xmean{xtarg} - xstd{xtarg};
            xerr(2,:) = xmean{xtarg} + xstd{xtarg};
            
            for i = 1:2, fit_err(:,i) = interp1(x_tick, xerr(i,:), fitx,'spline'); end
            
            in_between = [fit_err(:,1)', fliplr(fit_err(:,2)')];
            h = fill([fitx fliplr(fitx)], in_between, xcond_color{xtarg});
            set(h,'facealpha', .2)
            
            for i=1:2, plot(fitx, fit_err(:,i), '-', 'Color', xcond_color{xtarg},'LineWidth',line_w_err); end
        end
    end
    
    set(gca,'xlim', x_lim, 'ylim', y_lim_corr);
    set(gca,'XTick', x_tick, 'XTickLabel', x_tick-1)
    set(gca,'YTick', y_lim_corr(1):0.5:y_lim_corr(end))
    
    title(sprintf('Zscored Timecourse: target-only (N=%s)', num2str(length(xsub_groups))));
    xlabel('Volume (tr)');
    ylabel('classifier evidence');
    
    %*************** onset lines
    plot(x_lim, [0 0], ':', 'Color', 'r', 'LineWidth', 1.5)%baseline
    plot([dur_stim dur_stim] + 1, y_lim_corr,'-', 'Color', 'r')
    plot([dur_stim+dur_manipulation, dur_stim+dur_manipulation] + 1, y_lim_corr,'-', 'Color', 'r')
    
    text(1, y_lim_corr(1) + xgap, 'stim onset', 'Color','r');
    text(dur_stim + 1, y_lim_corr(1) + xgap, 'operation onset','Color','r');
    text(dur_stim + dur_manipulation + 1, y_lim_corr(1) + xgap, 'fixation onset','Color','r');
    
    %% ============= SAVE FIGURE
    if std_fill
        fig_fname = fullfile(xout_dir, ...
            sprintf('grp_plot_corrected_continous_operation_target_timecourse_random_%s_%s', ...
            xcond_name{xcond}, base_name));
    else
        fig_fname = fullfile(xout_dir, ...
            sprintf('grp_plot_line_corrected_continous_operation_target_timecourse_random_%s_%s', ...
            xcond_name{xcond}, base_name));
    end
    
    savefig(LogRegs_fig, sprintf('%s.fig', fig_fname));
    set(gcf,'PaperUnits','inches','PaperPosition',fig_rect/100)
    print('-djpeg', sprintf('%s.jpg',fig_fname), '-r100')
    saveas(LogRegs_fig, sprintf('%s.jpg',fig_fname), 'jpg')
    
    close(LogRegs_fig);

end
end