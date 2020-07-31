function[] = grp_analysis_timecourse_workingmemory_add(args, grp_mvpaout_add, dirs)

%---------------------------------------------------------------------------
%*************** group: timecourse of the working memory (category|subcategory)
%---------------------------------------------------------------------------

xph = args.xphase;
fprintf('\n(+) timecourse of working memory: %s\n', args.phase_name{xph});

xsub_groups = args.filtered_subs;%args.g_sub;

%% ============= UNPACK PARAMETERS
xindex          = args.index{xph};% param index from study

xparam          = xindex.param;
condition_names = {'Maintain','RepCat','RepSubcat','Suppress','Clear'};
n_condition     = length(condition_names);
dur_stim        = xparam.dur_stim;
dur_manip       = xparam.dur_manipulation;
n_tc_trs        = args.tc_tr_disp;
line_w          = 2; 
dur_sync        = args.dur_sync_disp;

n_tr_blks       = args.n_tr_blks;%trs in one stats block: 5=2.3s/3=1.38s
n_stat_blks     = (n_tc_trs/n_tr_blks);
xalpha          = args.alpha;
xcorr_alpha     = xalpha/n_stat_blks;

%*************** BASELINE CORRECTION
baseline_trs    = args.baseline_trs;

%*************** output basename
base_name       = sprintf('%s_n%s', args.analysis_basename, num2str(length(xsub_groups)));

%*************** same/different label
if strcmp(args.level, 'category')
    it_sames   = 2; 
    same_names = {'same','different'};
    it_subplot = [1 2];
elseif strcmp(args.level, 'subcategory')
    it_sames   = 3;
    same_names = {'same','different','related'};
    it_subplot = [1 3 2];
end

%% ============= LOAD PARSED TIMECOURSE
fname     = fullfile(dirs.mvpa.group.out{xph}, sprintf('grp_wmem_timecourse_blk%d_%s.mat', n_tr_blks, base_name));
load(fname);% grp_timecourse

fname_add = fullfile(dirs.mvpa.group.out{xph}, sprintf('grp_wmem_timecourse_add_blk%d_%s.mat', n_tr_blks, base_name));

%% ============= TIMECOURSE: SAME SEQUENCE N=N+1: PER CONDITION
%*************** working memory contents
% same category/subcategory for N+1 trial: collecting only target evidence
% condition: 1_maintain,2_replace_category,3_replace_item,4_target_suppress,5_global_clear
% decode.timecourse.sameseq{xsame}.condition{xcond}.evidence{xtr}
% decode.timecourse.sameseq{xsame}.sync.condition{xcond}.evidence{xtr}
% xsame: category: 1_same, 2_differ /subcategory: 1_same, 2_differ, 3_related

clear timecourse timecourse_bcorr

for xcond = 1:n_condition
    for xtr = 1:n_tc_trs
        
        for xsame = 1:it_sames
            %*************** random effect
            xevidence = []; xevidence_sync = [];
            
            for xsub = xsub_groups
                xevidence = horzcat(xevidence, ...
                    mean(grp_mvpaout_add{xsub}.decode.timecourse.sameseq{xsame}.condition{xcond}.evidence{xtr}));
                
                if xtr <= dur_sync
                    xevidence_sync = horzcat(xevidence_sync, ...
                        mean(grp_mvpaout_add{xsub}.decode.timecourse.sameseq{xsame}.sync.condition{xcond}.evidence{xtr}));
                end
            end
            
            timecourse.sameseq.cond{xcond}.same{xsame}.tr{xtr} = xevidence;
            timecourse.sameseq.cond{xcond}.mean(xsame, xtr)    = mean(xevidence);
            timecourse.sameseq.cond{xcond}.se(xsame, xtr)      = std(xevidence)/sqrt(length(xsub_groups));
            
            %*************** sync
            if xtr <= dur_sync
                timecourse.sameseq.sync.cond{xcond}.same{xsame}.tr{xtr} = xevidence_sync;
                timecourse.sameseq.sync.cond{xcond}.mean(xsame, xtr)    = mean(xevidence_sync);
                timecourse.sameseq.sync.cond{xcond}.se(xsame, xtr)      = std(xevidence_sync)/sqrt(length(xsub_groups));
            end
        end
    end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%============= BASELINE CORRECTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for xcond = 1:n_condition
    %*************** baseline_trs: 3TRs
    xbase_evidence = zeros(1, baseline_trs);
    
    for xtr = 1:baseline_trs
        xbase_evidence(xtr) = ...
            mean(grp_timecourse.tc.evidence{xcond}.targ{1}.tr{xtr}); %#ok<*NODEF>
    end
    
    mean_baseline = mean(xbase_evidence);
    
    for xsame = 1:it_sames
        for xtr = 1:n_tc_trs
            
            xevidence = ...
                timecourse.sameseq.cond{xcond}.same{xsame}.tr{xtr} - mean_baseline;
            
            timecourse_bcorr.sameseq.cond{xcond}.same{xsame}.tr{xtr} = xevidence;
            timecourse_bcorr.sameseq.cond{xcond}.mean(xsame, xtr)    = mean(xevidence);
            timecourse_bcorr.sameseq.cond{xcond}.se(xsame, xtr)      = std(xevidence)/sqrt(length(xsub_groups));
            
            %*************** sync
            if xtr <= dur_sync
                
                xevidence_sync = ...
                    timecourse.sameseq.sync.cond{xcond}.same{xsame}.tr{xtr} - mean_baseline;
                
                timecourse_bcorr.sameseq.sync.cond{xcond}.same{xsame}.tr{xtr} = xevidence_sync;
                timecourse_bcorr.sameseq.sync.cond{xcond}.mean(xsame, xtr)    = mean(xevidence_sync);
                timecourse_bcorr.sameseq.sync.cond{xcond}.se(xsame, xtr)      = std(xevidence_sync)/sqrt(length(xsub_groups));
                
            end
        end
    end
end

%% ============= STATS: B/W SAME/DIFFERENT/RELATED
for xcond = 1:n_condition
    for xblk = 1:n_stat_blks
        
        timecourse_bcorr.sameseq.cond{xcond}.t_test.blk{xblk}.pvalue = ...
            nan(it_sames-1, it_sames-1);
        
        it_trs = (1:n_tr_blks) + (n_tr_blks * (xblk-1));
        
        for xcol = 1:it_sames-1
            for xrow = xcol+1:it_sames
                
                xtarg_col = []; xtarg_row = [];
                
                for xtr = it_trs
                    xtarg_col = vertcat(xtarg_col, timecourse_bcorr.sameseq.cond{xcond}.same{xcol}.tr{xtr});
                    xtarg_row = vertcat(xtarg_row, timecourse_bcorr.sameseq.cond{xcond}.same{xrow}.tr{xtr});
                end
                
                %*************** ttest
                [~, xpvalue, ~, xstats] = ttest(mean(xtarg_col), mean(xtarg_row), 'Alpha', xcorr_alpha);
                
                timecourse_bcorr.sameseq.cond{xcond}.t_test.blk{xblk}.pvalue(xrow, xcol) = xpvalue;
                timecourse_bcorr.sameseq.cond{xcond}.t_test.blk{xblk}.stats{xrow, xcol}  = xstats;
            end
        end
    end
    
    %*************** sync
    for xblk = 1:(dur_sync/n_tr_blks)
        
        timecourse_bcorr.sameseq.sync.cond{xcond}.t_test.blk{xblk}.pvalue = ...
            nan(it_sames-1, it_sames-1);
        
        it_trs = (1:n_tr_blks) + (n_tr_blks * (xblk-1));
        
        for xcol = 1:it_sames-1
            for xrow = xcol+1:it_sames
                
                xtarg_col = []; xtarg_row = [];
                
                for xtr = it_trs
                    xtarg_col = vertcat(xtarg_col, timecourse_bcorr.sameseq.sync.cond{xcond}.same{xcol}.tr{xtr});
                    xtarg_row = vertcat(xtarg_row, timecourse_bcorr.sameseq.sync.cond{xcond}.same{xrow}.tr{xtr});
                end
                
                %*************** ttest
                [~, xpvalue, ~, xstats] = ttest(mean(xtarg_col), mean(xtarg_row), 'Alpha', xcorr_alpha);
                
                timecourse_bcorr.sameseq.sync.cond{xcond}.t_test.blk{xblk}.pvalue(xrow, xcol) = xpvalue;
                timecourse_bcorr.sameseq.sync.cond{xcond}.t_test.blk{xblk}.stats{xrow, xcol}  = xstats;
            end
        end
    end
end

%% ============= STATS: B/W TARGETS
for xsame = 1:it_sames
    for xblk = 1:n_stat_blks
        
        timecourse_bcorr.sameseq.same{xsame}.t_test.blk{xblk}.pvalue = ...
            nan(n_condition-1, n_condition-1);
        
        it_trs = (1:n_tr_blks) + (n_tr_blks * (xblk-1));
        
        for xcol = 1:n_condition-1
            for xrow = xcol+1:n_condition
                
                xtarg_col = []; xtarg_row = [];
                
                for xtr = it_trs
                    xtarg_col = vertcat(xtarg_col, timecourse_bcorr.sameseq.cond{xcol}.same{xsame}.tr{xtr});
                    xtarg_row = vertcat(xtarg_row, timecourse_bcorr.sameseq.cond{xrow}.same{xsame}.tr{xtr});
                end
                
                %*************** ttest
                [~, xpvalue, ~, xstats] = ttest(mean(xtarg_col), mean(xtarg_row), 'Alpha', xcorr_alpha);
                
                timecourse_bcorr.sameseq.same{xsame}.t_test.blk{xblk}.pvalue(xrow, xcol) = xpvalue;
                timecourse_bcorr.sameseq.same{xsame}.t_test.blk{xblk}.stats{xrow, xcol}  = xstats;
            end
        end
    end
    
    %*************** sync
    for xblk = 1:(dur_sync/n_tr_blks)
        
        timecourse_bcorr.sameseq.sync.same{xsame}.t_test.blk{xblk}.pvalue = ...
            nan(n_condition-1, n_condition-1);
        
        it_trs = (1:n_tr_blks) + (n_tr_blks * (xblk-1));
        
        for xcol = 1:n_condition-1
            for xrow = xcol+1:n_condition
                
                xtarg_col = []; xtarg_row = [];
                
                for xtr = it_trs
                    xtarg_col = vertcat(xtarg_col, timecourse_bcorr.sameseq.sync.cond{xcol}.same{xsame}.tr{xtr});
                    xtarg_row = vertcat(xtarg_row, timecourse_bcorr.sameseq.sync.cond{xrow}.same{xsame}.tr{xtr});
                end
                
                %*************** ttest
                [~, xpvalue, ~, xstats] = ttest(mean(xtarg_col), mean(xtarg_row), 'Alpha', xcorr_alpha);
                
                timecourse_bcorr.sameseq.sync.same{xsame}.t_test.blk{xblk}.pvalue(xrow, xcol) = xpvalue;
                timecourse_bcorr.sameseq.sync.same{xsame}.t_test.blk{xblk}.stats{xrow, xcol}  = xstats;
            end
        end
    end 
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%============= REPLACE CONTROL: TIMECOURSE
if strcmp(args.level, 'category')
    it_conds = [4, 5];
    for xit = 1:length(it_conds)
        for xtr = 1:n_tc_trs
            xevidence = [];
            for xsub = xsub_groups
                xevidence = horzcat(xevidence, ...
                    mean(grp_mvpaout_add{xsub}.decode.control.condition{xit}.new_evidence{xtr}));
            end
            timecourse.new_control.cond{xit}.tr{xtr} = xevidence;
        end
    end
    
    %%============= BASELINE CORRECTION
    for xit = 1:length(it_conds)
        %*************** baseline_trs: 3TRs
        xbase_evidence = zeros(1, baseline_trs);
        
        for xtr = 1:baseline_trs
            xbase_evidence(xtr) = ...
                mean(timecourse.new_control.cond{xit}.tr{xtr});
        end
        
        mean_baseline = mean(xbase_evidence);
        
        for xtr = 1:n_tc_trs
            timecourse_bcorr.new_control.cond{xit}.tr{xtr} = ...
                timecourse.new_control.cond{xit}.tr{xtr} - mean_baseline;
        end
    end
    
    %%============= STATS: B/W TARGETS
    % grp_timecourse: 2. replace category: {1} target {2} new category
    % targ: 1: item 2: new, 3: rand_supp, 4: rand_clear
    n_tr_blks   = args.n_tr_blks;%trs in one stats block: 5=2.3s/3=1.38s
    anova_trs   = args.anova_trs;
    n_stat_blks = (length(anova_trs)/n_tr_blks);

    n_targ = 5;
    timecourse_bcorr.new_control.targ{1}.tr = grp_timecourse.tc_bcorr.evidence{2}.targ{1}.tr;
    timecourse_bcorr.new_control.targ{2}.tr = grp_timecourse.tc_bcorr.evidence{2}.targ{2}.tr;
    timecourse_bcorr.new_control.targ{3}.tr = timecourse_bcorr.new_control.cond{1}.tr;
    timecourse_bcorr.new_control.targ{4}.tr = timecourse_bcorr.new_control.cond{2}.tr;
    
    for xtr = 1:n_tc_trs
        timecourse_bcorr.new_control.targ{5}.tr{xtr} = ...
            mean([timecourse_bcorr.new_control.targ{3}.tr{xtr}; timecourse_bcorr.new_control.targ{4}.tr{xtr}]);
    end

    for xtarg = 1:n_targ
        for xtr = 1:n_tc_trs
            xevidence = timecourse_bcorr.new_control.targ{xtarg}.tr{xtr};
            timecourse_bcorr.new_control.targ{xtarg}.mean(xtr) = mean(xevidence);
            timecourse_bcorr.new_control.targ{xtarg}.se(xtr) = std(xevidence)/sqrt(length(xsub_groups));
        end
    end
    
    for xblk = 1:n_stat_blks
        
        it_trs = (1:n_tr_blks) + (n_tr_blks * (xblk-1));
        
        for xtarg = 1:n_targ, xtarg_evi{xtarg} = []; end
        
        for xtr = it_trs
            for xtarg = 1:n_targ
                xtarg_evi{xtarg} = vertcat(xtarg_evi{xtarg}, timecourse_bcorr.new_control.targ{xtarg}.tr{xtr});
            end
        end
        
        for xcol = 1:(n_targ-1)
            for xrow = (xcol + 1):n_targ
                %*************** ttest
                % targ: 1: item 2: new, 3: rand_supp, 4: rand_clear
                [~, xpvalue, ~, xstats] = ttest(mean(xtarg_evi{xcol}), mean(xtarg_evi{xrow}), 'Alpha', xcorr_alpha);
                
                timecourse_bcorr.new_control.t_test.blk{xblk}.pvalue(xrow, xcol) = xpvalue;
                timecourse_bcorr.new_control.t_test.blk{xblk}.stats{xrow, xcol}  = xstats;
            end
        end
    end
end

%% ============= SAVE
grp_timecourse_add.tc       = timecourse;
grp_timecourse_add.tc_bcorr = timecourse_bcorr; %#ok<*STRNU>

save(fname_add, 'grp_timecourse_add');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BASELINE CORRECTED PLOTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ============= PLOT PARAMS
xcond_color  = args.cond_color;
xonset_color = args.onset_color;
xcolor{1}    = [232 14 138]/255;%same
xcolor{2}    = [101 47 142]/255;%differ
xcolor{3}    = [0 166 156]/255;%related

x_tick       = 1:n_tc_trs;
x_ticklable  = 1:n_tc_trs;
x_lim        = [x_tick(1) x_tick(end)];% + [-0.5 0.5];
y_lim        = [-0.3 0.5];

col_w = 800; row_h = 400;

%% ============= REPLACE CONTROL: TIMECOURSE
if strcmp(args.level, 'category')
    % targ: 1: item 2: new, 3: rand_supp, 4: rand_clear
    n_targ = 4;
    for xtarg = 1:n_targ
        targ_color{xtarg} = xcond_color{xtarg + 1};%target
    end
    
    fig_rect = [0 0 1200 900];
    xfig     = figure;
    set(xfig, 'Position', fig_rect)
        
    %*************** legend
    for xtarg = 1:n_targ
        plot([0 0], [0 0], '-','Color', targ_color{xtarg}); hold on;
    end
    xlegend     = {'replace-target','replace-new','suppress-rand','clear-rand'};
    lg          = legend(xlegend);
    lg.Location = 'SouthEast';
    legend(xlegend,'AutoUpdate','off')
    
    %*************** mean line plots
    fitx = linspace(1, n_tc_trs, n_tc_trs*10);

    for xtarg = 1:n_targ
        clear xmean xse xerr fit_err fity
        xmean = timecourse_bcorr.new_control.targ{xtarg}.mean;
        xse   = timecourse_bcorr.new_control.targ{xtarg}.se;
        
        fity  = interp1(x_tick, xmean, fitx,'spline');
        plot(fitx, fity, '-','Color', targ_color{xtarg}, 'LineWidth', line_w); hold on;
        
        %*************** std error-bar filling
        xerr(1,:) = xmean - xse;
        xerr(2,:) = xmean + xse;
        
        for i = 1:2, fit_err(:,i) = interp1(x_tick, xerr(i,:), fitx,'spline'); end
        
        in_between = [fit_err(:,1)', fliplr(fit_err(:,2)')];
        h = fill([fitx fliplr(fitx)], in_between, targ_color{xtarg});
        set(h,'facealpha', .4)
        
        plot(fitx, fit_err(:,1), '-', 'Color', targ_color{xtarg})
        plot(fitx, fit_err(:,2), '-', 'Color', targ_color{xtarg})
    end
    
    %*************** stats stars
    n = 0;
    for xcol = 1:(n_targ-1)
        for xrow = (xcol + 1):n_targ
            xheight = (y_lim(2)-0.03) - (n * 0.03);
            text(-3, xheight, sprintf('%s vs. %s', xlegend{xcol}, xlegend{xrow}));
            plot(x_lim, [xheight xheight], ':', 'Color', [0.75 0.75 0.75])
            
            t_pvalue = [];
            it_blk = 3:7;
            for xblk = it_blk
                t_pvalue = horzcat(t_pvalue, timecourse_bcorr.new_control.t_test.blk{xblk}.pvalue(xrow, xcol));
            end
            
            [~, ~, ~, xadj_p] = fdr_bh(t_pvalue, args.alpha, 'pdep', 'no');
            
            for xit = 1:length(it_blk)
                xblk = it_blk(xit);
                
                xtr = (xblk-1) * n_tr_blks;
                plot([xtr xtr]+0.5, y_lim, '--', 'Color', [0.75 0.75 0.75])
                
                xpvalue = xadj_p(xit);
                
                if xpvalue <= xcorr_alpha
                    text(xtr+(n_tr_blks/2)+0.5, xheight, '*', 'FontSize', 20);
                elseif xpvalue <= (xcorr_alpha/5)%0.01
                    text(xtr+(n_tr_blks/2)+0.4, xheight, '**', 'FontSize', 20);
                elseif xpvalue <= (xcorr_alpha/50)%0.001
                    text(xtr+(n_tr_blks/2)+0.3, xheight, '***', 'FontSize', 20);
                end
                
                if (xcol==(n_targ-1)) && (xrow==n_targ)
                    xstats = timecourse_bcorr.new_control.t_test.blk{xblk}.stats{xrow, xcol};
                    xstat_txt = sprintf('T(%s)=%4.4f, p=%1.4f', num2str(xstats.df), xstats.tstat, xpvalue);
                    h = text(xtr+(n_tr_blks/2)+0.5, xheight - 0.05, xstat_txt, 'FontSize', 10);
                    set(h, 'Rotation', -90);
                end
                
            end
            n = n + 1;
            
        end
    end
    
    xtr = xblk * n_tr_blks;
    plot([xtr xtr]+0.5, y_lim, '--', 'Color', [0.75 0.75 0.75])
            
    set(gca,'xlim', x_lim, 'ylim', y_lim);
    set(gca,'XTick', x_tick, 'YTick', y_lim(1):0.1:y_lim(end))
    set(gca,'XTickLabel', x_ticklable)
    
    title(sprintf('%s classifier: replace control (N=%s, p<%1.4f)', ...
        args.level, num2str(length(xsub_groups)), xcorr_alpha));
    xlabel('Volume (tr)');
    ylabel('classifier evidence');
    
    %%*************** save fig
    fig_fname = fullfile(dirs.mvpa.group.out{xph}, ...
        sprintf('grp_plot_replace_control1_wn_tc_blk%d_%s', ...
        n_tr_blks, base_name));
    
    savefig(xfig, sprintf('%s.fig', fig_fname));
    set(gcf,'PaperUnits','inches','PaperPosition', fig_rect/100)
    print('-dpng', '-r300', sprintf('%s.png',fig_fname))
    saveas(xfig, sprintf('%s.png',fig_fname), 'png')
    
%     close(xfig);
    
    %% ============= MEAN SUPPRESS & CLEAR
    % targ: 1: item 2: new, 3: rand_supp, 4: rand_clear
    fig_rect = [0 0 1200 900];
    xfig     = figure;
    set(xfig, 'Position', fig_rect)
        
    targ_color{5} = [100 100 100]/255;
    it_targs = [1, 2, 5];
    %*************** legend
    for xtarg = it_targs
        plot([0 0], [0 0], '-','Color', targ_color{xtarg}); hold on;
    end
    xlegend     = {'replace-target','replace-new','removal-rand'};
    lg          = legend(xlegend);
    lg.Location = 'SouthEast';
    legend(xlegend,'AutoUpdate','off')
    
    %*************** mean line plots
    fitx = linspace(1, n_tc_trs, n_tc_trs*10);

    for xtarg = it_targs
        clear xmean xse xerr fit_err fity
        xmean = timecourse_bcorr.new_control.targ{xtarg}.mean;
        xse   = timecourse_bcorr.new_control.targ{xtarg}.se;
        
        fity  = interp1(x_tick, xmean, fitx,'spline');
        plot(fitx, fity, '-','Color', targ_color{xtarg}, 'LineWidth', line_w); hold on;
        
        %*************** std error-bar filling
        xerr(1,:) = xmean - xse;
        xerr(2,:) = xmean + xse;
        
        for i = 1:2, fit_err(:,i) = interp1(x_tick, xerr(i,:), fitx,'spline'); end
        
        in_between = [fit_err(:,1)', fliplr(fit_err(:,2)')];
        h = fill([fitx fliplr(fitx)], in_between, targ_color{xtarg});
        set(h,'facealpha', .4)
        
        plot(fitx, fit_err(:,1), '-', 'Color', targ_color{xtarg})
        plot(fitx, fit_err(:,2), '-', 'Color', targ_color{xtarg})
    end
    
    %*************** stats stars
    n = 0;
    for it_col = 1:(length(it_targs)-1)
        xcol = it_targs(it_col);
        for it_row = (it_col + 1):length(it_targs)
            xrow = it_targs(it_row);
            
            xheight = (y_lim(2)-0.03) - (n * 0.03);
            text(-3, xheight, sprintf('%s vs. %s', xlegend{it_col}, xlegend{it_row}));
            plot(x_lim, [xheight xheight], ':', 'Color', [0.75 0.75 0.75])
            
            t_pvalue = [];
            it_blk = 3:7;
            for xblk = it_blk
                t_pvalue = horzcat(t_pvalue, timecourse_bcorr.new_control.t_test.blk{xblk}.pvalue(xrow, xcol));
            end
            
            [~, ~, ~, xadj_p] = fdr_bh(t_pvalue, args.alpha, 'pdep', 'no');
            
            for xit = 1:length(it_blk)
                xblk = it_blk(xit);
                
                xtr = (xblk-1) * n_tr_blks;
                plot([xtr xtr]+0.5, y_lim, '--', 'Color', [0.75 0.75 0.75])
                
                xpvalue = xadj_p(xit);
                
                if xpvalue <= xcorr_alpha
                    text(xtr+(n_tr_blks/2)+0.5, xheight, '*', 'FontSize', 20);
                elseif xpvalue <= (xcorr_alpha/5)%0.01
                    text(xtr+(n_tr_blks/2)+0.4, xheight, '**', 'FontSize', 20);
                elseif xpvalue <= (xcorr_alpha/50)%0.001
                    text(xtr+(n_tr_blks/2)+0.3, xheight, '***', 'FontSize', 20);
                end
                
                if (it_col==(length(it_targs)-1)) && (it_row==length(it_targs))
                    xstats = timecourse_bcorr.new_control.t_test.blk{xblk}.stats{xrow, xcol};
                    xstat_txt = sprintf('T(%s)=%4.4f, p=%1.4f', num2str(xstats.df), xstats.tstat, xpvalue);
                    h = text(xtr+(n_tr_blks/2)+0.5, xheight - 0.05, xstat_txt, 'FontSize', 10);
                    set(h, 'Rotation', -90);
                end
            end
            n = n + 1;
            
            xtr = xblk * n_tr_blks;
            plot([xtr xtr]+0.5, y_lim, '--', 'Color', [0.75 0.75 0.75])
        end
    end
    
    set(gca,'xlim', x_lim, 'ylim', y_lim);
    set(gca,'XTick', x_tick, 'YTick', y_lim(1):0.1:y_lim(end))
    set(gca,'XTickLabel', x_ticklable)
    
    title(sprintf('%s classifier: replace control (N=%s, p<%1.4f)', ...
        args.level, num2str(length(xsub_groups)), xcorr_alpha));
    xlabel('Volume (tr)');
    ylabel('classifier evidence');
    
    %%*************** save fig
    fig_fname = fullfile(dirs.mvpa.group.out{xph}, ...
        sprintf('grp_plot_replace_control2_wn_tc_blk%d_%s', ...
        n_tr_blks, base_name));
    
    savefig(xfig, sprintf('%s.fig', fig_fname));
    set(gcf,'PaperUnits','inches','PaperPosition', fig_rect/100)
    print('-dpng', '-r300', sprintf('%s.png',fig_fname))
    saveas(xfig, sprintf('%s.png',fig_fname), 'png')
    
%     close(xfig);
    
end


%% ============= WITHIN CONDITION
% same vs. different vs. related

fig_rect = [0 0 col_w row_h];

for xcond = 1:n_condition
    
    xfig     = figure;
    set(xfig, 'Position', fig_rect)

    clear xmean xstd fity
    
    %*************** mean line plots
    fitx = linspace(1, n_tc_trs, n_tc_trs*10);
    
    for xsame = 1:it_sames
        xmean{xsame} = timecourse_bcorr.sameseq.cond{xcond}.mean(xsame, x_tick);
        xstd{xsame}  = timecourse_bcorr.sameseq.cond{xcond}.se(xsame, x_tick); %#ok<*AGROW>
        
        fity{xsame}  = interp1(x_tick, xmean{xsame}, fitx,'spline');
        
        plot(fitx, fity{xsame}, '-','Color', xcolor{xsame}, 'LineWidth', line_w); hold on;
    end
    
    %*************** legend
    xlegend        = same_names;
    lg             = legend(same_names);
    lg.Location    = 'SouthEast';
    legend(xlegend,'AutoUpdate','off')
    
    %*************** stats stars
    n = 0;
    for xcol = 1:(it_sames-1)
        for xrow = (xcol+1):it_sames
            xheight = (y_lim(2)-0.05) - (n * 0.04);
            text(-3, xheight, sprintf('%s vs. %s', ...
                same_names{xcol}, same_names{xrow}));
            plot(x_lim, [xheight xheight], ':', 'Color', [0.75 0.75 0.75])
            
            for xblk = 1:n_stat_blks
                
                xtr = (xblk-1) * n_tr_blks;
                plot([xtr xtr]+0.5, y_lim, '--', 'Color', [0.75 0.75 0.75])
                
                xpvalue = timecourse_bcorr.sameseq.cond{xcond}.t_test.blk{xblk}.pvalue(xrow, xcol);
                
                if xpvalue <= xcorr_alpha
                    text(xtr+(n_tr_blks/2)+0.5, xheight, '*', 'FontSize', 20);
                elseif xpvalue <= (xcorr_alpha/5)%0.01
                    text(xtr+(n_tr_blks/2)+0.4, xheight, '**', 'FontSize', 20);
                elseif xpvalue <= (xcorr_alpha/50)%0.001
                    text(xtr+(n_tr_blks/2)+0.3, xheight, '***', 'FontSize', 20);
                end
            end
            n = n + 1;
        end
    end
    
    %*************** std error-bar filling
    for xsame = 1:it_sames
        clear xerr fit_err
        xerr(1,:) = xmean{xsame} - xstd{xsame};
        xerr(2,:) = xmean{xsame} + xstd{xsame};
        
        for i = 1:2, fit_err(:,i) = interp1(x_tick, xerr(i,:), fitx,'spline'); end
        
        in_between = [fit_err(:,1)', fliplr(fit_err(:,2)')];
        h = fill([fitx fliplr(fitx)], in_between, xcolor{xsame});
        set(h,'facealpha', .4)
        
        plot(fitx, fit_err(:,1), '-', 'Color', xcolor{xsame})
        plot(fitx, fit_err(:,2), '-', 'Color', xcolor{xsame})
    end
    
    %*************** real onset lines
    plot([dur_stim dur_stim]+1, y_lim, '-', 'Color', xonset_color)
    plot([(dur_stim + dur_manip) (dur_stim + dur_manip)]+1, y_lim, '-', 'Color', xonset_color)
    
    set(gca,'xlim', x_lim, 'ylim', y_lim);
    set(gca,'XTick', x_tick, 'YTick', y_lim(1):0.1:y_lim(end))
    set(gca,'XTickLabel', x_ticklable)
    
    title(sprintf('%s classifier: %s (N=%s, p<%1.4f)', ...
        args.level, condition_names{xcond}, num2str(length(xsub_groups)), xcorr_alpha));
    xlabel('Volume (tr)');
    ylabel('classifier evidence');
    
    h1 = text(1.5, y_lim(1)+0.05, '(1) stim onset', 'Color', xonset_color);
    h2 = text(dur_stim + 1.5, y_lim(1)+0.05, '(2) operation onset', 'Color', xonset_color);
    h3 = text(dur_stim + dur_manip + 1.5, y_lim(1)+0.05, '(3) fixation onset', 'Color', xonset_color);
    set(h1,'Rotation', 90); set(h2,'Rotation', 90); set(h3,'Rotation', 90);
    
    %% *************** save fig
    fig_fname = fullfile(dirs.mvpa.group.out{xph}, ...
        sprintf('grp_plot_sameseq_wn_tc_%s_blk%d_%s', ...
        condition_names{xcond}, n_tr_blks, base_name));
    
    savefig(xfig, sprintf('%s.fig', fig_fname));
    set(gcf,'PaperUnits','inches','PaperPosition', fig_rect/100)
    print('-dpng', '-r300', sprintf('%s.png',fig_fname))
    saveas(xfig, sprintf('%s.png',fig_fname), 'png')
    
    close(xfig);
    
end

%% ============= BETWEEN CONDITION
fig_rect     = [0 0 col_w (row_h * it_sames)];
it_fill      = 1;

sel_conds{1} = 1:n_condition;
sel_conds{2} = [1 4];
sel_conds{3} = [2 3];
sel_conds{4} = [1 2];
sel_conds{5} = [4 5];
sel_conds{6} = [2 4];
sel_conds{7} = [2 5];
    
for xfill = it_fill
    for xsel = 1:length(sel_conds)
        xsel_conds = sel_conds{xsel};
        
        xfig = figure;
        set(xfig, 'Position', fig_rect)
        
        %*************** mean line plots
        fitx = linspace(1, n_tc_trs, n_tc_trs*10);
        
        for xsame = 1:it_sames
            clear xmean xstd fity
        
            subplot(it_sames, 1, it_subplot(xsame))
            
            for xcond = xsel_conds
            
                xmean{xcond} = timecourse_bcorr.sameseq.cond{xcond}.mean(xsame, x_tick);
                xstd{xcond}  = timecourse_bcorr.sameseq.cond{xcond}.se(xsame, x_tick); %#ok<*AGROW>
                
                fity{xcond}  = interp1(x_tick, xmean{xcond}, fitx,'spline');
                
                plot(fitx, fity{xcond}, '-','Color', xcond_color{xcond}, 'LineWidth', line_w); hold on;
            end
            
            %*************** legend
            clear xlegend
            for i=1:length(xsel_conds); xlegend{i} = condition_names{xsel_conds(i)}; end
            lg             = legend(xlegend);
            lg.Location    = 'BestOutside';
            legend(xlegend,'AutoUpdate','off')
            
            %*************** stats stars
            n = 0;
            for i = 1:(length(xsel_conds)-1)
                xcol = xsel_conds(i);
                for j = i+1:length(xsel_conds)
                    xrow    = xsel_conds(j);
                    xheight = (y_lim(2)-0.05) - (n * 0.04);
                    text(-3, xheight, sprintf('%d vs. %d', xrow, xcol));
                    plot(x_lim, [xheight xheight], ':', 'Color', [0.75 0.75 0.75])
                    
                    for xblk = 1:n_stat_blks
                        
                        xtr = (xblk-1) * n_tr_blks;
                        plot([xtr xtr]+0.5, y_lim, '--', 'Color', [0.75 0.75 0.75])
                        
                        xpvalue = timecourse_bcorr.sameseq.same{xsame}.t_test.blk{xblk}.pvalue(xrow, xcol);
                        
                        if xpvalue <= xcorr_alpha
                            text(xtr+(n_tr_blks/2)+0.5, xheight, '*', 'FontSize', 20);
                        elseif xpvalue <= (xcorr_alpha/5)%0.01
                            text(xtr+(n_tr_blks/2)+0.4, xheight, '**', 'FontSize', 20);
                        elseif xpvalue <= (xcorr_alpha/50)%0.001
                            text(xtr+(n_tr_blks/2)+0.3, xheight, '***', 'FontSize', 20);
                        end
                    end
                    n = n + 1;
                end
            end
            
            %*************** std error-bar filling
            if xfill
                for xcond = xsel_conds
                    clear xerr fit_err
                    xerr(1,:) = xmean{xcond} - xstd{xcond};
                    xerr(2,:) = xmean{xcond} + xstd{xcond};
                    
                    for i = 1:2, fit_err(:,i) = interp1(x_tick, xerr(i,:), fitx,'spline'); end
                    
                    in_between = [fit_err(:,1)', fliplr(fit_err(:,2)')];
                    h = fill([fitx fliplr(fitx)], in_between, xcond_color{xcond});
                    set(h,'facealpha', .2)
                    
                    plot(fitx, fit_err(:,1), '-', 'Color', xcond_color{xcond}, 'LineWidth', line_w)
                    plot(fitx, fit_err(:,2), '-', 'Color', xcond_color{xcond}, 'LineWidth', line_w)
                end
            end
            
            %*************** baseline 0
            plot(x_lim, [0 0], ':', 'Color', 'r', 'LineWidth', 1.5)%baseline
            
            %*************** real onset lines
            plot([dur_stim dur_stim]+1, y_lim, '-','Color', xonset_color)
            plot([(dur_stim + dur_manip) (dur_stim + dur_manip)]+1, y_lim, '-', 'Color', xonset_color)
            
            set(gca,'xlim', x_lim, 'ylim', y_lim);
            set(gca,'XTick', x_tick,'YTick', y_lim(1):0.1:y_lim(end))
            set(gca,'XTickLabel', x_ticklable)
            
            title(sprintf('%s classifier: %s (N=%s, p<%1.4f)', ...
                args.level, same_names{xsame}, num2str(length(xsub_groups)), xcorr_alpha));
            xlabel('Volume (tr)');
            ylabel('classifier evidence (target)');
            
            %*************** real onset
            h1 = text(1.5, y_lim(1)+0.01, '(1) stim onset', 'Color', xonset_color);
            h2 = text(dur_stim + 1.5, y_lim(1)+0.01, '(2) operation onset', 'Color', xonset_color);
            h3 = text(dur_stim + dur_manip + 1.5, y_lim(1)+0.01, '(3) fixation onset', 'Color', xonset_color);
            set(h1,'Rotation', 90); set(h2,'Rotation', 90); set(h3,'Rotation', 90);
        end
        
        %% *************** save fig
        if xfill
            fig_fname = fullfile(dirs.mvpa.group.out{xph}, ...
                sprintf('grp_plot_sameseq_wn_targets_sel%d_blk%d_%s', xsel, n_tr_blks, base_name));
        else
            fig_fname = fullfile(dirs.mvpa.group.out{xph}, ...
                sprintf('grp_plot_sameseq_wn_targets_line_sel%d_blk%d_%s', xsel, n_tr_blks, base_name));
        end
        
        savefig(xfig, sprintf('%s.fig', fig_fname));
        set(gcf,'PaperUnits','inches','PaperPosition', fig_rect/100)
        print('-dpng', '-r300', sprintf('%s.png',fig_fname))
        saveas(xfig, sprintf('%s.png',fig_fname), 'png')
        
        close(xfig);
        
    end
end

%% ============= SYNCH: BETWEEN CONDITION
fig_rect    = [0 0 col_w/2 (row_h * it_sames)];
x_lim_sync  = [1 dur_sync];
x_tick_sync = x_lim_sync(1):x_lim_sync(end);

for xfill = it_fill
    for xsel = 1:length(sel_conds)
        xsel_conds = sel_conds{xsel};
        
        xfig = figure;
        set(xfig, 'Position', fig_rect)
        
        %*************** mean line plots
        fitx = linspace(1, n_tc_trs, n_tc_trs*10);
        
        for xsame = 1:it_sames
            clear xmean xstd fity
        
            subplot(it_sames, 1, it_subplot(xsame))
            
            for xcond = xsel_conds
            
                xmean{xcond} = timecourse_bcorr.sameseq.sync.cond{xcond}.mean(xsame, x_tick_sync);
                xstd{xcond}  = timecourse_bcorr.sameseq.sync.cond{xcond}.se(xsame, x_tick_sync); %#ok<*AGROW>
                
                fity{xcond}  = interp1(x_tick_sync, xmean{xcond}, fitx,'spline');
                
                plot(fitx, fity{xcond}, '-','Color', xcond_color{xcond}, 'LineWidth', line_w); hold on;
            end
            
            %*************** legend
            clear xlegend
            for i=1:length(xsel_conds); xlegend{i} = condition_names{xsel_conds(i)}; end
            lg             = legend(xlegend);
            lg.Location    = 'BestOutside';
            legend(xlegend,'AutoUpdate','off')
            
            %*************** stats stars
            n = 0;
            for i = 1:(length(xsel_conds)-1)
                xcol = xsel_conds(i);
                for j = i+1:length(xsel_conds)
                    xrow    = xsel_conds(j);
                    xheight = (y_lim(2)-0.05) - (n * 0.04);
                    text(-3, xheight, sprintf('%d vs. %d', xrow, xcol));
                    plot(x_lim_sync, [xheight xheight], ':', 'Color', [0.75 0.75 0.75])
                    
                    for xblk = 1:(dur_sync/n_tr_blks)
                        
                        xtr = (xblk-1) * n_tr_blks;
                        plot([xtr xtr]+0.5, y_lim, '--', 'Color', [0.75 0.75 0.75])
                        
                        xpvalue = timecourse_bcorr.sameseq.sync.same{xsame}.t_test.blk{xblk}.pvalue(xrow, xcol);
                        
                        if xpvalue <= xcorr_alpha
                            text(xtr+(n_tr_blks/2)+0.5, xheight, '*', 'FontSize', 20);
                        elseif xpvalue <= (xcorr_alpha/5)%0.01
                            text(xtr+(n_tr_blks/2)+0.4, xheight, '**', 'FontSize', 20);
                        elseif xpvalue <= (xcorr_alpha/50)%0.001
                            text(xtr+(n_tr_blks/2)+0.3, xheight, '***', 'FontSize', 20);
                        end
                    end
                    n = n + 1;
                end
            end
            
            %*************** std error-bar filling
            if xfill
                for xcond = xsel_conds
                    clear xerr fit_err
                    xerr(1,:) = xmean{xcond} - xstd{xcond};
                    xerr(2,:) = xmean{xcond} + xstd{xcond};
                    
                    for i = 1:2, fit_err(:,i) = interp1(x_tick_sync, xerr(i,:), fitx,'spline'); end
                    
                    in_between = [fit_err(:,1)', fliplr(fit_err(:,2)')];
                    h = fill([fitx fliplr(fitx)], in_between, xcond_color{xcond});
                    set(h,'facealpha', .2)
                    
                    plot(fitx, fit_err(:,1), '-', 'Color', xcond_color{xcond}, 'LineWidth', line_w)
                    plot(fitx, fit_err(:,2), '-', 'Color', xcond_color{xcond}, 'LineWidth', line_w)
                end
            end
            
            %*************** baseline 0
            plot(x_lim_sync, [0 0], ':', 'Color', 'r', 'LineWidth', 1.5)%baseline
            
            %*************** real onset lines
            plot([dur_stim dur_stim]+1, y_lim, '-','Color', xonset_color)
            
            set(gca,'xlim', x_lim_sync, 'ylim', y_lim);
            set(gca,'XTick', x_tick_sync,'YTick', y_lim(1):0.1:y_lim(end))
            set(gca,'XTickLabel', x_ticklable)
            
            title(sprintf('%s classifier: %s (N=%s, p<%1.4f)', ...
                args.level, same_names{xsame}, num2str(length(xsub_groups)), xcorr_alpha));
            xlabel('Volume (tr)');
            ylabel('classifier evidence (target)');
            
            %*************** real onset
            h1 = text(1.5, y_lim(1)+0.01, '(1) stim onset', 'Color', xonset_color);
            h2 = text(dur_stim + 1.5, y_lim(1)+0.01, '(2) operation onset', 'Color', xonset_color);
            set(h1,'Rotation', 90); set(h2,'Rotation', 90); 
        end
        
        %% *************** save fig
        if xfill
            fig_fname = fullfile(dirs.mvpa.group.out{xph}, ...
                sprintf('grp_plot_sameseq_wn_sync_targets_sel%d_blk%d_%s', xsel, n_tr_blks, base_name));
        else
            fig_fname = fullfile(dirs.mvpa.group.out{xph}, ...
                sprintf('grp_plot_sameseq_wn_sync_targets_line_sel%d_blk%d_%s', xsel, n_tr_blks, base_name));
        end
        
        savefig(xfig, sprintf('%s.fig', fig_fname));
        set(gcf,'PaperUnits','inches','PaperPosition', fig_rect/100)
        print('-dpng', '-r300', sprintf('%s.png',fig_fname))
        saveas(xfig, sprintf('%s.png',fig_fname), 'png')
        
        close(xfig);
        
    end
end

end