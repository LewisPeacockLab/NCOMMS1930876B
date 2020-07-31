function[] = grp_analysis_timecourse_workingmemory_sorted_cate(args, grp_mvpaout, dirs)

%---------------------------------------------------------------------------
%*************** group: timecourse of the working memory (category|subcategory)
%% ---------------------------------------------------------------------------

xph = args.xphase;
fprintf('\n(+) timecourse of working memory: %s\n', args.phase_name{xph});

xsub_groups     = args.filtered_subs;%args.g_sub;

for select_category = 1:3
    cate_array      = 1:3;
    new_category    = cate_array(~ismember(cate_array, select_category));
    
    %% ============= UNPACK PARAMETERS
    xindex          = args.index{xph};% param index from study
    
    xparam          = xindex.param;
    condition_names = {'Maintain','RepCat','RepSubcat','Suppress','Clear'};
    dur_stim        = xparam.dur_stim;
    dur_manip       = xparam.dur_manipulation;
    n_tc_trs        = args.tc_tr_disp;
    ttest_trs       = args.ttest_trs;
    anova_trs       = args.anova_trs;
    line_w          = 1;
    
    n_tr_blks       = args.n_tr_blks;%trs in one stats block: 5=2.3s/3=1.38s
    n_stat_blks     = (length(anova_trs)/n_tr_blks);
    xalpha          = args.alpha;
    
    %*************** BASELINE CORRECTION
    baseline_trs    = args.baseline_trs;
    
    %*************** output basename
    base_name       = sprintf('%s_n%s', args.analysis_basename, num2str(length(xsub_groups)));
    
    %% ============= PLOT PARAMS
    xcond_color  = args.cond_color;
    xonset_color = args.onset_color;
    
    %% ============= PARSE TIMECOURSE
    % fname = fullfile(dirs.mvpa.group.out{xph}, sprintf('grp_wmem_sorted_timecourse_blk%d_%s.mat', n_tr_blks, base_name));
    %% ============= TIMECOURSE: RANDOM EFFECT
    %*************** working memory contents
    % condition: 1_maintain,2_replace_category,4_target_suppress,5_global_clear
    % xtr: 12 presentation + 5:9 fixation = 21 max
    
    clear tc tc_diff xevidence
    
    xtarg     = 1;
    it_conds  = [2 4 5];
    all_conds = [1 it_conds];
    
    for xcond = all_conds
        %*************** baseline
        xmean_base = zeros(1, length(baseline_trs));
        for xtr = 1:baseline_trs
            xbase = zeros(1, length(xsub_groups));
            for it = 1:length(xsub_groups)
                xsub = xsub_groups(it);
                if xcond==2
                    tbase = [];
                    for xnew = new_category
                        tbase = horzcat(tbase, ...
                            mean(grp_mvpaout{xsub}.decode.cate_timecourse.condition{xcond}.category{select_category}.newcategory{xnew}.target{xtarg}.evidence{xtr}));
                    end
                    xbase(it) = mean(tbase);
                else
                    xbase(it) = ...
                        mean(grp_mvpaout{xsub}.decode.cate_timecourse.condition{xcond}.category{select_category}.target{xtarg}.evidence{xtr});
                end
            end
            
            xmean_base(xtr) = mean(xbase);
        end
        
        for xtr = 1:n_tc_trs
            %*************** baseline corrected
            tevidence = zeros(1, length(xsub_groups));
            for it = 1:length(xsub_groups)
                xsub = xsub_groups(it);
                if xcond==2
                    ttevidence = [];
                    for xnew = new_category
                        ttevidence = horzcat(ttevidence, ...
                            mean(grp_mvpaout{xsub}.decode.cate_timecourse.condition{xcond}.category{select_category}.newcategory{xnew}.target{xtarg}.evidence{xtr}));
                    end
                    tevidence(it) = mean(ttevidence);
                else
                    tevidence(it) = ...
                        mean(grp_mvpaout{xsub}.decode.cate_timecourse.condition{xcond}.category{select_category}.target{xtarg}.evidence{xtr});
                end
            end
            
            xevidence = tevidence - mean(xmean_base);
            tc.evidence{xcond}.tr(:, xtr) = xevidence;
            tc.mean(xcond, xtr)           = mean(xevidence);
            tc.se(xcond, xtr)             = std(xevidence)/sqrt(length(xsub_groups));
        end
    end
    
    %% *************** ANOVA
    for xblk = 1:n_stat_blks
        clear xevidence xdata xheader xtable
        
        %*************** collecting evidence
        xevidence = zeros(length(xsub_groups),length(all_conds));
        xheader   = cell(1, 1+length(all_conds));
        
        it_trs = (anova_trs(1)-1) + (1:n_tr_blks) + (n_tr_blks * (xblk-1));
        
        for it = 1:length(all_conds)
            xcond = all_conds(it);
            
            t_evidence = zeros(length(xsub_groups),length(it_trs));
            for xtr = it_trs
                t_evidence(:, xtr) = tc.evidence{xcond}.tr(:, xtr);
            end
            
            xevidence(:, it) = mean(t_evidence, 2);
            xheader{it+1}    = condition_names{xcond};
        end
        
        %*************** data table for ANOVA
        xdata(:,1) = xsub_groups';
        xdata(:, 2:size(xevidence, 2)+1) = xevidence;
        xheader{1} = 'SN';
        xtable     = array2table(xdata, 'VariableNames', xheader);
        
        xanova     = {'target_evidence'};
        xmeasures  = table((1:length(all_conds))','VariableNames', xanova);
        xrepmeas   = fitrm(xtable,'Maintain-Clear~1','WithinDesign',xmeasures);
        xanova_out = ranova(xrepmeas);
        
        xtx_out = sprintf('F(%s, %s)=%4.4f, p=%1.4f', ...
            num2str(xanova_out.DF(1)), num2str(xanova_out.DF(2)),...
            xanova_out.F(1), xanova_out.pValue(1));
        
        tc.anova.blk{xblk}.pvalue = xanova_out.pValue(1);
        tc.anova.blk{xblk}.stats  = xtx_out;
        
    end
    
    %% ============= TIMECOURSE DIFFERENCE
    clear xevidence
    
    for xcond = it_conds
        for xtr = 1:n_tc_trs
            xevidence = tc.evidence{xcond}.tr(:, xtr) - tc.evidence{1}.tr(:, xtr);
            
            tc_diff.evidence{xcond}.tr(:, xtr) = xevidence;
            tc_diff.mean(xcond, xtr)           = mean(xevidence);
            tc_diff.se(xcond, xtr)             = std(xevidence)/sqrt(length(xsub_groups));
        end
    end
    
    %*************** one-sample ttest with mean=0
    for it = 1:length(it_conds)
        xcond = it_conds(it);
        
        %*************** collecting evidence
        for xblk = 1:n_stat_blks
            clear t_evidence xevidence
            it_trs = (anova_trs(1)-1) + (1:n_tr_blks) + (n_tr_blks * (xblk-1));
            
            t_evidence = zeros(length(xsub_groups),length(it_trs));
            for i = 1:length(it_trs)
                xtr = it_trs(i);
                t_evidence(:, i) = tc_diff.evidence{xcond}.tr(:, xtr);
            end
            
            xevidence = mean(t_evidence, 2);
            
            [~, xpvalue, ~, xstats] = ttest(xevidence, 0, 'Alpha', args.alpha);
            
            tc_diff.ttest.cond{xcond}.pvalue.blk(xblk) = xpvalue;
            tc_diff.ttest.cond{xcond}.stats.blk{xblk}  = ...
                sprintf('T(%s)=%4.4f, p=%1.4f', num2str(xstats.df), xstats.tstat, xpvalue);
        end
    end
    
    %*************** TR-BASIS one-sample ttest with mean=0
    for xcond = it_conds
        for xtr = ttest_trs
            xevidnece = tc_diff.evidence{xcond}.tr(:, xtr);
            
            [~, xpvalue, ~, xstats] = ttest(xevidnece, 0, 'Alpha', args.alpha);
            
            tc_diff.ttest.cond{xcond}.pvalue.tr(xtr) = xpvalue;
            tc_diff.ttest.cond{xcond}.stats.tr{xtr}  = ...
                sprintf('T(%s)=%4.4f, p=%1.4f', num2str(xstats.df), xstats.tstat, xpvalue);
        end
    end
    
    %*************** pairwise ttest
    for xblk = 1:n_stat_blks
        
        tc_diff.pwttest.blk{xblk}.pvalue = nan(max(it_conds), max(it_conds));
        tc_diff.pwttest.blk{xblk}.stats  = cell(max(it_conds), max(it_conds));
        
        it_trs = (anova_trs(1)-1) + (1:n_tr_blks) + (n_tr_blks * (xblk-1));
        
        for xrow = it_conds
            for xcol = it_conds
                xtarg_col = []; xtarg_row = [];
                
                for xtr = it_trs
                    xtarg_row = horzcat(xtarg_row, tc_diff.evidence{xrow}.tr(:, xtr));
                    xtarg_col = horzcat(xtarg_col, tc_diff.evidence{xcol}.tr(:, xtr));
                end
                
                [~, xpvalue, ~, xstats] = ttest(mean(xtarg_row, 2), mean(xtarg_col, 2), 'Alpha', args.alpha);
                
                tc_diff.pwttest.blk{xblk}.pvalue(xrow, xcol) = xpvalue;
                tc_diff.pwttest.blk{xblk}.stats{xrow, xcol}  = ...
                    sprintf('T(%s)=%4.4f, p=%1.4f', num2str(xstats.df), xstats.tstat, xpvalue);
            end
        end
    end
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % TIMECOURSE PLOTS: BASELINE CORRECTED
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    sel_conds{1} = [2 4];
    sel_conds{2} = [2 5];
    sel_conds{3} = [4 5];
    
    x_tick       = 1:n_tc_trs;
    x_ticklable  = 1:n_tc_trs;
    x_lim        = [x_tick(1)-1 x_tick(end)+1];% + [-0.5 0.5];
    fig_rect     = [0 0 1200 1800];
    
    xfig = figure;
    set(xfig, 'Position', fig_rect)
    
    for xsubplt = [1 3:5]
        
        %*************** subplot
        if xsubplt==1, subplot(5,1, 1:2);
        else,          subplot(5,1, xsubplt); end
        
        %%%%%%%%%%%%%%%% timecourse
        %*************** mean line plots
        clear xmean xstd fity tt pp xlegend
        
        fitx    = linspace(1, n_tc_trs, n_tc_trs*10);
        
        if xsubplt==1
            xconds = all_conds;
            xmean  = tc.mean;
            xse    = tc.se;
            if strcmp(args.level, 'category')
                y_lim  = [-0.2 0.5];
            elseif strcmp(args.level, 'subcategory')
                y_lim  = [-0.2 0.3];
            end
        else
            xconds = sel_conds{xsubplt-2};
            xmean  = tc_diff.mean;
            xse    = tc_diff.se;
            y_lim  = [-0.3 0.2];
        end
        
        xlegend = cell(1, length(xconds));
        
        for it = 1:length(xconds)
            xcond       = xconds(it);
            xlegend{it} = condition_names{xcond};
            
            fity = interp1(x_tick, xmean(xcond, x_tick), fitx,'spline');
            
            plot(fitx, fity, '-','Color', xcond_color{xcond}, 'LineWidth', line_w); hold on;
        end
        
        %*************** legend
        if xsubplt==1
            lg = legend(xlegend{1}, xlegend{2}, xlegend{3}, xlegend{4});
            legend(xlegend{1}, xlegend{2}, xlegend{3}, xlegend{4},'AutoUpdate','off')
        else
            lg = legend(xlegend{1}, xlegend{2});
            legend(xlegend{1}, xlegend{2},'AutoUpdate','off')
        end
        lg.Location = 'SouthEast';
        lg.FontSize = 10;
        
        %*************** std error-bar filling
        for xcond = xconds
            clear xerr fit_err
            xerr(1,:) = xmean(xcond, x_tick) - xse(xcond, x_tick);
            xerr(2,:) = xmean(xcond, x_tick) + xse(xcond, x_tick);
            
            for i = 1:2, fit_err(:,i) = interp1(x_tick, xerr(i,:), fitx,'spline'); end %#ok<*AGROW>
            
            in_between = [fit_err(:,1)', fliplr(fit_err(:,2)')];
            h = fill([fitx fliplr(fitx)], in_between, xcond_color{xcond});
            h.FaceAlpha = .2;
            h.EdgeAlpha = .2;
            h.EdgeColor = xcond_color{xcond};
        end
        
        %*************** stats
        clear pvalue_matrix
        for xblk = 1:n_stat_blks
            if xsubplt==1
                %*************** ANOVA
                xpvalue = tc.anova.blk{xblk}.pvalue;
            else
                %*************** ttest
                xpvalue = tc_diff.pwttest.blk{xblk}.pvalue(xconds(1), xconds(2));
            end
            
            pvalue_matrix(xblk) = xpvalue;
        end
        
        %*************** FDR multiple comparison correction across timepoints (blocks)
        % method: 'pdep' the original Bejnamini & Hochberg FDR procedure is used,
        %          which is guaranteed to be accurate if the individual tests are independent or
        %          positively dependent (e.g., Gaussian variables that are positively correlated or independent)
        
        [xh, xcrit_p, xadj_ci_cvrg, xadj_p] = fdr_bh(pvalue_matrix, args.alpha, 'pdep', 'no');
        
        multi_out = sprintf('critical p=%1.4f, CI coverage=%1.4f', xcrit_p, xadj_ci_cvrg);
        
        if sum(xh)
            xp_max = max(pvalue_matrix(xh));
            xwhich = pvalue_matrix==xp_max;
            
            %*************** summary
            if xsubplt==1
                text(0.5, y_lim(2) - 0.1, ...
                    sprintf('%s, Max: %s', multi_out, tc.anova.blk{xwhich}.stats));
            else
                text(0.5, y_lim(2) - 0.1, ...
                    sprintf('%s, Min: %s', multi_out, tc_diff.pwttest.blk{xwhich}.stats{xconds(1), xconds(2)}));
            end
        end
        
        for xblk = 1:n_stat_blks
            %*************** block boundary
            xx = (anova_trs(1)-1) + 1 + (n_tr_blks * (xblk-1));
            plot([xx xx]-0.5, y_lim, '--', 'Color', [0.75 0.75 0.75])
            
            if xblk==n_stat_blks
                xx_r = anova_trs(1) + (n_tr_blks * (xblk));
                plot([xx_r xx_r]-0.5, y_lim, '--', 'Color', [0.75 0.75 0.75])
            end
            
            if xh(xblk)
                if xsubplt==1
                    %*************** ANOVA
                    xstats = tc.anova.blk{xblk}.stats;
                else
                    %*************** ttest
                    xstats = tc_diff.pwttest.blk{xblk}.stats{xconds(1), xconds(2)};
                end
                
                text(xx, y_lim(2) - 0.05, '*', 'FontSize', 20);
                h = text(xx, y_lim(1) + 0.05, sprintf('%s\nadj: p=%s', xstats, num2str(xadj_p(xblk))), 'FontSize', 8);
                set(h, 'Rotation', 90);
            else
                h = text(xx, y_lim(1) + 0.05, sprintf('adj: p=%s', num2str(xadj_p(xblk))), 'FontSize', 8);
                set(h, 'Rotation', 90);
            end
        end
        
        %*************** baseline 0
        plot(x_lim, [0 0], '--', 'Color', 'k', 'LineWidth', 1.5)%baseline
        
        set(gca,'xlim', x_lim, 'ylim', y_lim);
        set(gca,'XTick', x_tick,'YTick', y_lim(1):0.05:y_lim(end))
        set(gca,'XTickLabel', x_ticklable)
        
        title(sprintf('%s classifier: %s (N=%s, FDR alpha=%1.2f, TR/blk=%s)', ...
            args.level, category_names{select_category}, num2str(length(xsub_groups)), xalpha, num2str(n_tr_blks)));
        xlabel('Volume (tr)');
        ylabel('classifier evidence (target)');
        
        if xsubplt==1
            %*************** real onset lines
            plot([dur_stim dur_stim]+1, y_lim, '--','Color', xonset_color)
            plot([(dur_stim + dur_manip) (dur_stim + dur_manip)]+1, y_lim, '--', 'Color', xonset_color)
            
            %*************** real onset
            h1 = text(1.5, y_lim(1)+0.01, '(1) stim onset', 'Color', xonset_color);
            h2 = text(dur_stim + 1.5, y_lim(1)+0.01, '(2) operation onset', 'Color', xonset_color);
            h3 = text(dur_stim + dur_manip + 1.5, y_lim(1)+0.01, '(3) fixation onset', 'Color', xonset_color);
            set(h1,'Rotation', 90); set(h2,'Rotation', 90); set(h3,'Rotation', 90);
        end
    end
    
    %*************** save fig
    fig_fname = fullfile(dirs.mvpa.group.out{xph}, ...
        sprintf('tc_cate%d_targets_blk%d_%s', select_category, n_tr_blks, base_name));
    savefig(xfig, sprintf('%s.fig', fig_fname));
    set(gcf,'PaperUnits','inches','PaperPosition', fig_rect/100)
    print('-dpng', '-r300', sprintf('%s.png',fig_fname))
    saveas(xfig, sprintf('%s.png', fig_fname), 'png')
    
    close(xfig);
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % TIMECOURSE PLOTS: ONE-SAMPLE TTEST
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    x_tick       = 1:n_tc_trs;
    x_ticklable  = 1:n_tc_trs;
    x_lim        = [x_tick(1)-1 x_tick(end)+1];% + [-0.5 0.5];
    y_lim        = [-0.3 0.2];
    fig_rect     = [0 0 1200 1000];
    
    xfig = figure;
    set(xfig, 'Position', fig_rect)
    
    for xsubplt = 1:length(it_conds)
        xcond = it_conds(xsubplt);
        
        subplot(length(it_conds), 1, xsubplt);
        
        %%%%%%%%%%%%%%%% timecourse
        %*************** mean line plots
        clear xmean xstd fity tt pp xlegend
        
        xmean = tc_diff.mean(xcond, x_tick);
        xse   = tc_diff.se(xcond, x_tick);
        fitx  = linspace(1, n_tc_trs, n_tc_trs*10);
        
        fity  = interp1(x_tick, xmean, fitx,'spline');
        
        plot(fitx, fity, '-','Color', xcond_color{xcond}, 'LineWidth', line_w); hold on;
        
        %*************** legend
        xlegend     = condition_names{xcond};
        lg          = legend(xlegend);
        lg.Location = 'SouthEast';
        lg.FontSize = 10;
        legend(xlegend,'AutoUpdate','off')
        
        %*************** std error-bar filling
        clear xerr fit_err
        xerr(1,:) = xmean - xse;
        xerr(2,:) = xmean + xse;
        
        for i = 1:2, fit_err(:,i) = interp1(x_tick, xerr(i,:), fitx,'spline'); end %#ok<*AGROW>
        
        in_between = [fit_err(:,1)', fliplr(fit_err(:,2)')];
        h = fill([fitx fliplr(fitx)], in_between, xcond_color{xcond});
        h.FaceAlpha = .2;
        h.EdgeAlpha = .2;
        h.EdgeColor = xcond_color{xcond};
        
        %*************** stats: one-sample ttest
        blk_array     = 1:n_stat_blks;
        pvalue_matrix = tc_diff.ttest.cond{xcond}.pvalue.blk(blk_array);
        
        %*************** FDR multiple comparison correction across timepoints (blocks)
        % method: 'pdep' the original Bejnamini & Hochberg FDR procedure is used,
        %          which is guaranteed to be accurate if the individual tests are independent or
        %          positively dependent (e.g., Gaussian variables that are positively correlated or independent)
        
        [xh, xcrit_p, xadj_ci_cvrg, xadj_p] = fdr_bh(pvalue_matrix, args.alpha, 'pdep', 'no');
        
        multi_out = sprintf('critical p=%1.4f, CI coverage=%1.4f', xcrit_p, xadj_ci_cvrg);
        
        if sum(xh)
            xp_max = max(pvalue_matrix(xh));
            xwhich = pvalue_matrix==xp_max;
            
            %*************** summary
            text(0.5, y_lim(2) - 0.1, ...
                sprintf('%s, Min: %s', multi_out, tc_diff.ttest.cond{xcond}.stats.blk{blk_array(xwhich)}));
        end
        
        for xblk = blk_array
            it_trs = (anova_trs(1)-1) + (1:n_tr_blks) + (n_tr_blks * (xblk-1));
            xtr    = it_trs(1);
            
            %*************** block boundary
            plot([xtr xtr]-0.5, y_lim, '--', 'Color', [0.75 0.75 0.75], 'LineWidth', 1)
            if xblk==length(blk_array)
                plot([xtr xtr] + n_tr_blks - 0.5, y_lim, '--', 'Color', [0.75 0.75 0.75], 'LineWidth', 1)
            end
            
            if xh(xblk), text(it_trs(2), y_lim(2) - 0.05, '*', 'FontSize', 20); end
            
            xstats = tc_diff.ttest.cond{xcond}.stats.blk{xblk};
            h = text(xtr, y_lim(1) + 0.05, sprintf('%s\nadj: p=%s', xstats, num2str(xadj_p(xblk))), 'FontSize', 8);
            set(h, 'Rotation', 90);
        end
        
        %*************** baseline 0
        plot(x_lim, [0 0], '--', 'Color', 'k', 'LineWidth', 1.5)%baseline
        
        set(gca,'xlim', x_lim, 'ylim', y_lim);
        set(gca,'XTick', x_tick,'YTick', y_lim(1):0.1:y_lim(end))
        set(gca,'XTickLabel', x_ticklable)
        
        title(sprintf('%s classifier: %s, %s (N=%s, FDR alpha=%1.2f, TR/blk=%s)', ...
            args.level, category_names{select_category}, xlegend, num2str(length(xsub_groups)), xalpha, num2str(n_tr_blks)));
        xlabel('Volume (tr)');
        ylabel('classifier evidence (target)');
        
        %*************** real onset lines
        plot([dur_stim dur_stim]+1, y_lim, '--','Color', xonset_color)
        plot([(dur_stim + dur_manip) (dur_stim + dur_manip)]+1, y_lim, '--', 'Color', xonset_color)
        
    end
    
    %*************** save fig
    fig_fname = fullfile(dirs.mvpa.group.out{xph}, ...
        sprintf('tc_cate%d_cond_blk%d_%s', select_category, n_tr_blks, base_name));
    savefig(xfig, sprintf('%s.fig', fig_fname));
    set(gcf,'PaperUnits','inches','PaperPosition', fig_rect/100)
    print('-dpng', '-r300', sprintf('%s.png',fig_fname))
    saveas(xfig, sprintf('%s.png', fig_fname), 'png')
    
    close(xfig);
    
end
end