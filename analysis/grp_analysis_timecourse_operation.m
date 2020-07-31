function[] = grp_analysis_timecourse_operation(args, grp_mvpaout, dirs)

%---------------------------------------------------------------------------
%*************** group: timecourse of operation during study
%---------------------------------------------------------------------------

xph = args.xphase;
fprintf('\n(+) timecourse of working memory: %s\n', args.phase_name{xph});

xsub_groups = args.filtered_subs;%args.g_sub
% xsub_groups = args.filtered_subs(~ismember(args.filtered_subs, [1 3 20 22 27 44 40 49 51 59]));


%% ============= UNPACK PARAMETERS
xindex           = args.index{xph};% param index from study

xparam           = xindex.param;
xcond_name       = xparam.conds_names;
n_condition      = length(xcond_name);
dur_stim         = xparam.dur_stim;
dur_manipulation = xparam.dur_manipulation;% 6 TR
tc_tr_disp       = args.tc_tr_disp;%for stats/display
% it_trs           = 1:args.tc_tr;%timecourse trs
it_trs           = 1:tc_tr_disp;%timecourse trs

n_tr_blks        = args.n_tr_blks;
n_stat_blks      = (tc_tr_disp/n_tr_blks);
xalpha           = args.alpha/n_stat_blks;

condition_names  = {'maintain','replace category',...
    'replace subcategory','target suppress','global clear'};

%*************** BASELINE CORRECTION
baseline_trs    = args.baseline_trs;

%*************** output basename
base_name       = sprintf('%s_n%s', args.analysis_basename, num2str(length(xsub_groups)));

%% ============= PLOT PARAMS
% xcolor{1}  = [238, 20, 91]/255;% targ
% xcolor{2}  = [144, 144, 144]/255;% non_target

x_tick     = 1:args.tc_tr_disp;
x_lim      = [x_tick(1)-1 x_tick(end)+1];% + [-0.5 0.5];
n_ticks    = args.tc_tr_disp;
line_w     = 2;
line_w_err = 0.1;
rect_w     = 550;
rect_h     = 500;
font_size  = 5;

%*************** separated operations
xcond_color  = args.cond_color;
% xcate_color  = args.cate_color;

%% ============= TIMECOURSE: ALL OPERATIONS
%*************** separated working memory operation
%*************** decode.timecourse.operation{xcond}.decoded_operation{xoper}.evidence{zz}.tr{xtr} 
clear tc_all

n_targ = n_condition;% all operations

for zz = 1:2%1_raw, 2_zscored
    for xcond = 1:n_condition
        clear mean_baseline
        %*************** random effect
        for xtr = it_trs
            for xtarg = 1:n_targ
                xevidence = [];
                
                for xsub = xsub_groups
                    xevidence = horzcat(xevidence, ...
                        mean(grp_mvpaout{xsub}.decode.timecourse.operation{xcond}.decoded_operation{xtarg}.evidence{zz}.tr{xtr}));
                end
                
                tc_all.cond{xcond}.evidence{zz}.d_cond{xtarg}.tr{xtr} = xevidence;
                tc_all.cond{xcond}.evidence{zz}.mean(xtarg, xtr)      = mean(xevidence);
                tc_all.cond{xcond}.evidence{zz}.se(xtarg, xtr)        = std(xevidence)/sqrt(length(xsub_groups));
            end
        end
        
        %*************** baseline correction
        xbase_evidence = zeros(n_targ, baseline_trs);
        for xtarg = 1:n_targ
            for xtr = 1:baseline_trs
                it_tr = xtr + dur_stim;
                xbase_evidence(xtarg, xtr) = mean(tc_all.cond{xcond}.evidence{zz}.d_cond{xtarg}.tr{it_tr});
            end
            
            mean_baseline(xtarg) = mean(xbase_evidence(xtarg,:));
        end
        
        %*************** corrected evidence
        for xtr = it_trs
            for xtarg = 1:n_targ
                xcorrect_evidence = ...
                    tc_all.cond{xcond}.evidence{zz}.d_cond{xtarg}.tr{xtr} - mean_baseline(xtarg);
                
                tc_all.bs_corr.cond{xcond}.evidence{zz}.d_cond{xtarg}.tr{xtr} = xcorrect_evidence;
                tc_all.bs_corr.cond{xcond}.evidence{zz}.mean(xtarg, xtr)      = mean(xcorrect_evidence);
                tc_all.bs_corr.cond{xcond}.evidence{zz}.se(xtarg, xtr)        = std(xcorrect_evidence)/sqrt(length(xsub_groups));
                
            end
        end
    end
    
    %% ============= STATS: PAIRWISE TTEST
    for xcond = 1:n_condition
        for xblk = 1:n_stat_blks
            xblk_tr = (1:n_tr_blks) + (n_tr_blks * (xblk-1));

            tc_all.bs_corr.cond{xcond}.evidence{zz}.ttest.blk{xblk}.pvalue = nan(n_targ-1, n_targ-1);
            
            for xcol = 1:(n_targ-1)
                for xrow = (xcol+1):n_targ
                    
                    xpat = tc_all.bs_corr.cond{xcond}.evidence{zz};
                                        
                    xtarg_col = []; xtarg_row = [];
                    
                    for xtr = xblk_tr
                        xtarg_col = horzcat(xtarg_col, xpat.d_cond{xcol}.tr{xtr});
                        xtarg_row = horzcat(xtarg_row, xpat.d_cond{xrow}.tr{xtr});
                    end
                    
                    %*************** ttest
                    [~, xpvalue, ~, xstats] = ttest(xtarg_row, xtarg_col, 'Alpha', 0.05);
                    
                    tc_all.bs_corr.cond{xcond}.evidence{zz}.ttest.blk{xblk}.pvalue(xrow, xcol) = xpvalue;
                    tc_all.bs_corr.cond{xcond}.evidence{zz}.ttest.blk{xblk}.stats{xrow, xcol}  = xstats;
                end
            end
        end
    end
end

%% ============= TIMECOURSE: TARGET-ONLY: ALL OPERATIONS
%*************** target timecourses
clear tc_targs

for zz = 1:2%1_raw, 2_zscored, 3_locked zscore
    clear mean_baseline xevidence xcorrect_evidence
    
    for xcond = 1:n_condition
        for xtr = it_trs
            xevidence = tc_all.cond{xcond}.evidence{zz}.d_cond{xcond}.tr{xtr};
            
            tc_targs.evidence{zz}.cond{xcond}.tr{xtr}   = xevidence;
            tc_targs.evidence{zz}.cond{xcond}.mean(xtr) = mean(xevidence);
            tc_targs.evidence{zz}.cond{xcond}.se(xtr)   = std(xevidence)/sqrt(length(xsub_groups));
        end
    end
    
    %*************** baseline correction
    xbase_evidence = zeros(n_condition, baseline_trs);
    for xcond = 1:n_condition
        for xtr = 1:baseline_trs
            it_tr = xtr + dur_stim;
            xbase_evidence(xcond, xtr) = mean(tc_targs.evidence{zz}.cond{xcond}.tr{it_tr});
        end
        
        mean_baseline(xcond) = mean(xbase_evidence(xcond,:));
    end
    
    %*************** corrected evidence
    for xcond = 1:n_condition
        for xtr = it_trs
            xcorrect_evidence = ...
                tc_targs.evidence{zz}.cond{xcond}.tr{xtr} - mean_baseline(xcond);
            
            tc_targs.bs_corr.evidence{zz}.cond{xcond}.tr{xtr}   = xcorrect_evidence;
            tc_targs.bs_corr.evidence{zz}.cond{xcond}.mean(xtr) = mean(xcorrect_evidence);
            tc_targs.bs_corr.evidence{zz}.cond{xcond}.se(xtr)   = std(xcorrect_evidence)/sqrt(length(xsub_groups));
        end
    end
end

%% ============= ZSCORING TIMECOURSE: LOCKED TIMECOURSE
%*************** target-only: zscored >> baseline corrected
zz    = 1;%evidence
it_zz = 3;%timecourse locked zscore
within_sub = 1;
clear xevidence xmean xsd mean_baseline

%*************** zscore
if within_sub
    for xcond = 1:n_condition
        xevidence = [];
        for xtr = it_trs
            xevidence = vertcat(xevidence, tc_targs.evidence{zz}.cond{xcond}.tr{xtr});
        end
        
        for xsub = 1:length(xsub_groups)
            z_evidence(:,xsub) = zscore(xevidence(:, xsub));
        end
        
        for xtr = it_trs
            tc_targs.evidence{it_zz}.cond{xcond}.tr{xtr} = z_evidence(xtr,:);
        end
    end
else
    for xcond = 1:n_condition
        xevidence = [];
        for xtr = it_trs
            xevidence = horzcat(xevidence, tc_targs.evidence{zz}.cond{xcond}.tr{xtr});
        end
        
        z_evidence = reshape(zscore(xevidence), [length(xsub_groups), length(it_trs)])';
        
        for xtr = it_trs
            tc_targs.evidence{it_zz}.cond{xcond}.tr{xtr} = z_evidence(xtr,:);
        end
    end
end

%*************** baseline correction
xbase_evidence = zeros(n_condition, baseline_trs);
for xcond = 1:n_condition
    for xtr = 1:baseline_trs
        it_tr = xtr + dur_stim;
        
        xbase_evidence(xcond, xtr) = mean(tc_targs.evidence{it_zz}.cond{xcond}.tr{it_tr});
    end
    
    mean_baseline(xcond) = mean(xbase_evidence(xcond,:));
end

%*************** zscore > corrected evidence
for xcond = 1:n_condition
    for xtr = it_trs
        
        xcorrect_evidence = tc_targs.evidence{it_zz}.cond{xcond}.tr{xtr} -  mean_baseline(xcond);
        
        tc_targs.bs_corr.evidence{it_zz}.cond{xcond}.tr{xtr}   = xcorrect_evidence;
        tc_targs.bs_corr.evidence{it_zz}.cond{xcond}.mean(xtr) = mean(xcorrect_evidence);
        tc_targs.bs_corr.evidence{it_zz}.cond{xcond}.se(xtr)   = std(xcorrect_evidence)/sqrt(length(xsub_groups));
        
    end
end

%% ============= STATS: BASELINE_CORRECTED
for zz = 1:3%1_raw, 2_zscored, 3_locked zscore
    %% ============= STATS: PAIRWISE TTEST
    n_targ = n_condition;
    
    for xblk = 1:n_stat_blks
        xblk_tr = (1:n_tr_blks) + (n_tr_blks * (xblk-1));
        
        tc_targs.bs_corr.evidence{zz}.ttest.blk{xblk}.pvalue = nan(n_targ, n_targ);
        
        for xcol = 1:(n_targ-1)
            for xrow = (xcol + 1):n_targ
                
                xpat = tc_targs.bs_corr.evidence{zz};
                
                xtarg_col = []; xtarg_row = [];
                
                for xtr = xblk_tr
                    xtarg_col = horzcat(xtarg_col, xpat.cond{xcol}.tr{xtr});
                    xtarg_row = horzcat(xtarg_row, xpat.cond{xrow}.tr{xtr});
                end
                
                %*************** ttest
                [~, xpvalue, ~, xstats] = ttest(xtarg_row, xtarg_col, 'Alpha', 0.05);
                
                tc_targs.bs_corr.evidence{zz}.ttest.blk{xblk}.pvalue(xrow, xcol) = xpvalue;
                tc_targs.bs_corr.evidence{zz}.ttest.blk{xblk}.stats{xrow, xcol}  = xstats;
            end
        end
    end

    %% ============= STATS: ANOVA: 5/4/3 operations
    %*************** anova{ff}: 1_5, 2_4, 1:3 operations (category-level)
    for ff = 1:3
        for xblk = 1:n_stat_blks
            xblk_tr = (1:n_tr_blks) + (n_tr_blks * (xblk-1));
            
            clear xevidence xheader xdata
            
            if ff==1
                it_conds = 1:n_condition;
                xheader  = {'SN','Maintain','RepCat','RepSubcat','Suppress','Clear'};
            elseif ff==2
                it_conds = [1 2 4 5]; 
                xheader  = {'SN','Maintain','RepCat','Suppress','Clear'};
            elseif ff==3
                it_conds = [1 2 3]; 
                xheader  = {'SN','Maintain','RepCat','RepSubcat'};
            end
            
            xevidence  = [];
            for it = 1:length(it_conds)
                xcond = it_conds(it);
                
                tevi = [];
                for xtr = xblk_tr
                    tevi = horzcat(tevi, tc_targs.bs_corr.evidence{zz}.cond{xcond}.tr{xtr}');
                end
                
                xevidence(:, it) = mean(tevi, 2);
            end
            
            xdata(:,1) = xsub_groups';
            xdata(:, 2:size(xevidence, 2)+1) = xevidence;
            
            xtable      = array2table(xdata, 'VariableNames', xheader);
            
            xanova      = {'target_evidence'};
            xmeasures   = table((1:length(it_conds))','VariableNames', xanova);
            
            xrepmeas    = fitrm(xtable, sprintf('%s-%s~1', xheader{2}, xheader{end}),'WithinDesign',xmeasures);
            xanova_out  = ranova(xrepmeas);
            
            var_names   = [xanova xanova_out.Properties.VariableNames];
            row_names   = xanova_out.Properties.RowNames;
            
            tmp_val     = [row_names num2cell(xanova_out.Variables)];
            xanva_table = cell2table(tmp_val, 'VariableNames', var_names);
            
            xout = sprintf('F(%s, %s)=%4.4f, p=%1.4f', ...
                num2str(xanova_out.DF(1)), num2str(xanova_out.DF(2)),...
                xanova_out.F(1), xanova_out.pValue(1));
            
            fprintf('%s\n', xout);
            
            tc_targs.bs_corr.evidence{zz}.anova{ff}.blk{xblk}.pvalue = xanova_out.pValue(1);
            tc_targs.bs_corr.evidence{zz}.anova{ff}.blk{xblk}.table  = xanva_table;
            tc_targs.bs_corr.evidence{zz}.anova{ff}.blk{xblk}.stats  = xout;
            
        end
    end
end

%% ============= SAVE
grp_timecourse.all   = tc_all;
grp_timecourse.targs = tc_targs; %#ok<*STRNU>

fname = fullfile(dirs.mvpa.group.out{xph}, sprintf('grp_operation_timecourse_%s.mat', base_name));
save(fname, 'grp_timecourse');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%============= PLOTTING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ============= PLOT: ALL OPERATIONS: BASELINE CORRECTED
%*************** column: 1_evidence, 2_zscored_evidence
clear y_lim
n_targ     = n_condition;
disp_stats = 1; %ttest
n_zz       = 2;

y_lim{1}   = [-0.4 0.6];
y_lim{2}   = [-2 1.5];%zscore
y_ticks{1} = y_lim{1}(1):0.1:y_lim{1}(end);
y_ticks{2} = y_lim{2}(1):0.5:y_lim{2}(end);
fig_rect   = [0 0 (rect_w * n_zz) (rect_h * n_condition)]; 

xfig = figure;
set(xfig, 'Position', fig_rect)

for std_fill = 1%0:1
    for zz = 1:n_zz
        for xcond = 1:n_condition
            clear xmean xstd fity
            
            subplot(n_condition, n_zz, zz + n_zz * (xcond-1));
%             xcond + (xcond-1) + (zz-1));
            
            %*************** mean line plots
            fitx = linspace(1, n_ticks, n_ticks*10);
            
            for xtarg = 1:n_targ
                xmean{xtarg} = tc_all.bs_corr.cond{xcond}.evidence{zz}.mean(xtarg, x_tick);
                xstd{xtarg}  = tc_all.bs_corr.cond{xcond}.evidence{zz}.se(xtarg, x_tick); %#ok<*AGROW>
                
                fity{xtarg}  = interp1(x_tick, xmean{xtarg}, fitx,'spline');
                
                if xtarg==xcond
                    plot(fitx, fity{xtarg}, '-', 'Color', xcond_color{xtarg}, 'LineWidth', line_w + 2); hold on;
                else
                    plot(fitx, fity{xtarg}, '-', 'Color', xcond_color{xtarg}, 'LineWidth', line_w); hold on;
                end
            end
            
            %*************** legend
            if zz==n_zz
                xlegend     = condition_names;
                lg          = legend(xlegend);
                lg.Location = 'SouthWest';
                lg.FontSize = font_size;
                
                legend(xlegend,'AutoUpdate','off')
            end
            grid on
                
            %*************** ttest
            if disp_stats
                n = 1;
                
                for xcol = 1:(n_targ-1)
                    for xrow = (xcol+1):n_targ
                        y_off   = (y_lim{zz}(2)-y_lim{zz}(1))/20;
                        yy_sig = y_lim{zz}(2) - (n * y_off);
                        text(-2, yy_sig + 0.02, sprintf('%s vs. %s', num2str(xcol), num2str(xrow)), 'FontSize', font_size);
                        plot(x_lim, [yy_sig yy_sig], ':', 'Color', [0.75 0.75 0.75])
                        
                        for xblk = 1:n_stat_blks
                            xp = tc_all.bs_corr.cond{xcond}.evidence{zz}.ttest.blk{xblk}.pvalue(xrow, xcol);
                            
                            if xp <= 0.1
                                if (xp <= (0.1/n_stat_blks)) && (xp > xalpha)
                                    xsig = '+';
                                elseif (xp <= xalpha) && (xp > (xalpha/5))
                                    xsig = '*';
                                elseif (xp <= (xalpha/5)) && (xp > (xalpha/50))
                                    xsig = '**';
                                elseif xp <= (xalpha/50)%0.001
                                    xsig = '***';
                                end
                                
                                text((xblk * n_tr_blks)-0.1, yy_sig, xsig, 'FontSize', font_size*2);
                                %                             text(xblk-0.1, yy_sig-y_off, sprintf('p=%4.4f', xp), 'FontSize', 5);
                            end
                        end
                        n = n + 1;
                    end
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
                    h.FaceAlpha = .2;
                    h.EdgeAlpha = .2;
                    h.EdgeColor = xcond_color{xtarg};
                    
                    for i=1:2, plot(fitx, fit_err(:,i), '-', 'Color', xcond_color{xtarg},'LineWidth',line_w_err); end
                end
            end
            
            set(gca,'xlim', x_lim, 'ylim', y_lim{zz});
            set(gca,'XTick', x_tick, 'YTick', y_ticks{zz}, 'XTickLabel', x_tick, 'FontSize', font_size)
            
            if zz==1
                title(sprintf('%s_z (N=%s)', condition_names{xcond}, num2str(length(xsub_groups))));
            else
                title(sprintf('%s (N=%s)', condition_names{xcond}, num2str(length(xsub_groups))));
            end
            xlabel('Volume (tr)');
            ylabel('classifier evidence');
            
            %*************** onset lines
            plot([dur_stim dur_stim] + 1, y_lim{zz},'-', 'Color', 'r')
            plot([dur_stim+dur_manipulation, dur_stim+dur_manipulation] + 1, y_lim{zz},'-', 'Color', 'r')
            
            text(1, y_lim{zz}(1)+0.05, 'stim', 'Color','r', 'FontSize', font_size);
            text(dur_stim + 1, y_lim{zz}(1)+0.05, 'operation','Color','r', 'FontSize', font_size);
            text(dur_stim + dur_manipulation + 1, y_lim{zz}(1)+0.05, 'fixation','Color','r', 'FontSize', font_size);
            
        end
    end
    
    %% ============= SAVE FIGURE
    if std_fill
        fig_fname = fullfile(dirs.mvpa.group.out{xph}, ...
            sprintf('grp_tc_all_%s', base_name));
    else
        fig_fname = fullfile(dirs.mvpa.group.out{xph}, ...
            sprintf('grp_tc_all_line_%s', base_name));
    end
    
    savefig(xfig, sprintf('%s.fig', fig_fname));
    set(gcf,'PaperUnits','inches','PaperPosition',fig_rect/100)
    print('-dpng', sprintf('%s.png',fig_fname), '-r300')
    saveas(xfig, sprintf('%s.png',fig_fname), 'png')
    
    close(xfig);
    
end

%% ============= PLOT: ALL OPERATIONS: BASELINE CORRECTED
%*************** column: 1_evidence, 2_zscored_evidence
clear y_lim 

sel_conds{1}  = 1:n_condition;
sel_conds{2}  = [1 2 4 5];
sel_conds{3}  = [1 2 3];
sel_conds{4}  = [1 4];
sel_conds{5}  = [1 2];
sel_conds{6}  = [2 3];
sel_conds{7}  = [2 4];
sel_conds{8}  = [2 5];
sel_conds{9}  = [4 5];

it_zz         = [1 3];%1:3
n_zz          = length(it_zz);
y_lim{1}      = [-0.2 0.5];
y_lim{2}      = [-0.5 1];%zscore
y_lim{3}      = [-2 3];%zscore
for zz = it_zz
    xgrid(zz)   = (y_lim{zz}(end)-y_lim{zz}(1))/10;
    y_ticks{zz} = y_lim{zz}(1):xgrid(zz):y_lim{zz}(end); 
end
font_size     = 10;
fig_rect      = [0 0 (rect_w * 3) rect_h]; 

for std_fill = 1%0:1
    for xsel = 1:length(sel_conds)    
        
        xfig = figure;
        set(xfig, 'Position', fig_rect)
        
        for ss = 1:length(it_zz)
            zz = it_zz(ss);
            clear xmean xstd fity
            
            xsel_conds = sel_conds{xsel};
            n_targ     = length(xsel_conds);
            
            subplot(1, n_zz, ss)
            
            %*************** mean line plots
            fitx = linspace(1, n_ticks, n_ticks*10);
            
            for xtarg = 1:n_targ
                xcond = xsel_conds(xtarg);
                
                xmean{xtarg} = tc_targs.bs_corr.evidence{zz}.cond{xcond}.mean(x_tick);
                xstd{xtarg}  = tc_targs.bs_corr.evidence{zz}.cond{xcond}.se(x_tick);
                
                fity{xtarg}  = interp1(x_tick, xmean{xtarg}, fitx,'spline');
                
                plot(fitx, fity{xtarg}, '-', 'Color', xcond_color{xcond}, 'LineWidth', line_w); hold on;
            end
            
            %*************** legend
            if zz==1
                clear xlegend
                for xtarg = 1:n_targ
                    xcond = xsel_conds(xtarg);
                    xlegend{xtarg} = condition_names{xcond};
                end
                lg             = legend(xlegend);
                lg.Location    = 'SouthWest';%'bestoutside';'best';%
                lg.FontSize    = font_size;
                legend(xlegend,'AutoUpdate','off')
            end
            grid on
           
            %*************** std error-bar filling
            if std_fill
                for xtarg = 1:n_targ
                    xcond = xsel_conds(xtarg);
                    
                    clear xerr
                    xerr(1,:) = xmean{xtarg} - xstd{xtarg};
                    xerr(2,:) = xmean{xtarg} + xstd{xtarg};
                    
                    for i = 1:2, fit_err(:,i) = interp1(x_tick, xerr(i,:), fitx,'spline'); end
                    
                    in_between = [fit_err(:,1)', fliplr(fit_err(:,2)')];
                    h = fill([fitx fliplr(fitx)], in_between, xcond_color{xcond});
                    h.FaceAlpha = 0.2;
                    h.EdgeColor = xcond_color{xcond};
                    h.EdgeAlpha = 0.2;
                end
            end
            
            %*************** ttest
            n = 1;
            for i = 1:(n_targ-1)
                xcol = xsel_conds(i);
                for j = (i+1):n_targ
                    xrow = xsel_conds(j);
                    
                    for ss = 1:2, tt{ss} = []; pp{ss} = []; end
                    
                    y_off  = xgrid(zz)/2;
                    yy_sig = y_lim{zz}(2) - (n * y_off);
                    
                    if n_targ > 2
                        text(-5, yy_sig, sprintf('%s vs. %s', ...
                            num2str(xcol), num2str(xrow)));
                    end
                    
                    plot(x_lim, [yy_sig yy_sig], ':', 'Color', [0.75 0.75 0.75])
                    
                    for xblk = 1:n_stat_blks
                        xp     = tc_targs.bs_corr.evidence{zz}.ttest.blk{xblk}.pvalue(xrow, xcol);
                        xstats = tc_targs.bs_corr.evidence{zz}.ttest.blk{xblk}.stats{xrow, xcol};
                        
                        if xp <= (0.1/n_stat_blks)
                            if (xp <= (0.1/n_stat_blks)) && (xp > xalpha)
                                xsig = '+';
                            elseif (xp <= xalpha) && (xp > (xalpha/5))
                                xsig = '*';
                            elseif (xp <= (xalpha/5)) && (xp > (xalpha/50))
                                xsig = '**';
                            elseif xp <= (xalpha/50)%0.001
                                xsig = '***';
                            end
                            
                            text((xblk * n_tr_blks)-1, yy_sig, xsig, 'FontSize', font_size*1.5);
                        
                            if (xsel~=1) && (xsel~=2)                                
                                if length(xsel_conds) < 3
                                    h = text((xblk * n_tr_blks)-1, y_lim{zz}(end) - (xgrid(zz)*4), sprintf('T=%4.4f, P=%4.4f', ...
                                        xstats.tstat,  xp), 'FontSize', 8);
                                    set(h,'Rotation',90);
                                    
                                    if xstats.tstat>0, ss = 1; else, ss = 2; end %#ok<*FXSET>
                                    tt{ss} = horzcat(tt{ss}, xstats.tstat);
                                    pp{ss} = horzcat(pp{ss}, xp);
                                end
                            end
                        end
                        
                        xx = n_tr_blks + (n_tr_blks * (xblk-1)) + 0.5;
                        plot([xx xx], y_lim{zz}, '--', 'Color', 'b')
                    end
                    n = n + 1;
                    
                    if length(xsel_conds) < 3
                        if ~isempty(tt{1})
                            text(1, 0.2 - (n * 0.05), ...
                                sprintf('Ts(%s)>%4.4f, Ps<%4.4f', num2str(xstats.df), min(tt{1}), max(pp{1})), 'FontSize', 10);
                        end
                        if  ~isempty(tt{2})
                            text(15, 0.2 - (n * 0.05), ...
                                sprintf('Ts(%s)<%4.4f, Ps<%4.4f', num2str(xstats.df), max(tt{2}), max(pp{2})), 'FontSize', 10);
                        end
                    end
                end
            end
            
            %*************** anova
            if xsel < 4
                
                yy_sig = y_lim{zz}(1) + (y_off*2);
                
                for xblk = 1:n_stat_blks
                    xp     = tc_targs.bs_corr.evidence{zz}.anova{xsel}.blk{xblk}.pvalue;
                    xstats = tc_targs.bs_corr.evidence{zz}.anova{xsel}.blk{xblk}.stats;
                    if xp <= (0.1/n_stat_blks)
                        if (xp <= 0.1) && (xp > xalpha)
                            xsig = '+';
                        elseif (xp <= xalpha) && (xp > (xalpha/5))
                            xsig = '*';
                        elseif (xp <= (xalpha/5)) && (xp > (xalpha/50))
                            xsig = '**';
                        elseif xp <= (xalpha/50)%0.001
                            xsig = '***';
                        end
                        
                        text((xblk * n_tr_blks)-1, yy_sig, xsig, 'FontSize', font_size*1.5);
                        h = text((xblk * n_tr_blks)-1, yy_sig+y_off, xstats, 'FontSize', 7);
                        set(h,'Rotation',90);
                    end
                end
            end
            
            set(gca,'xlim', x_lim, 'ylim', y_lim{zz});
            set(gca,'XTick', x_tick, 'YTick', y_ticks{zz}, 'XTickLabel', x_tick, 'FontSize', font_size)
            
            if zz==1
                title(sprintf('targets_z (N=%s)', num2str(length(xsub_groups))));
            else
                title(sprintf('targets (N=%s)', num2str(length(xsub_groups))));
            end
            xlabel('Volume (tr)');
            ylabel('classifier evidence');
            
            %*************** onset lines
            plot([dur_stim dur_stim] + 1, y_lim{zz},'-', 'Color', 'r')
            plot([dur_stim+dur_manipulation, dur_stim+dur_manipulation] + 1, y_lim{zz},'-', 'Color', 'r')
            
            t{1} = text(1.5, y_lim{zz}(1)+0.05, 'stim', 'Color','r', 'FontSize', font_size);
            t{2} = text(dur_stim + 1.5, y_lim{zz}(1)+0.05, 'operation','Color','r','FontSize', font_size);
            t{3} = text(dur_stim + dur_manipulation + 1.5, y_lim{zz}(1)+0.05, 'fixation','Color','r','FontSize', font_size);
            for i=1:3, set(t{i},'Rotation', 90); end
        end
        
        %% ============= SAVE FIGURE
        if std_fill
            fig_fname = fullfile(dirs.mvpa.group.out{xph}, ...
                sprintf('grp_sel%s_tc_targets_%s', num2str(xsel), base_name));
        else
            fig_fname = fullfile(dirs.mvpa.group.out{xph}, ...
                sprintf('grp_sel%s_tc_targets_line_%s', num2str(xsel), base_name));
        end
        
        savefig(xfig, sprintf('%s.fig', fig_fname));
        set(gcf,'PaperUnits','inches','PaperPosition',fig_rect/100)
        print('-dpng', sprintf('%s.png',fig_fname), '-r300')
        saveas(xfig, sprintf('%s.png',fig_fname), 'png')
        
        close(xfig);
        
    end
end

end