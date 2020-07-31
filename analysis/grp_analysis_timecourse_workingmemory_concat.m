function[] = grp_analysis_timecourse_workingmemory_concat(args, grp_mvpaout, dirs)

%---------------------------------------------------------------------------
%*************** group: timecourse of the working memory (category|subcategory)
%---------------------------------------------------------------------------

xph = args.xphase;
fprintf('\n(+) timecourse of working memory: %s\n', args.phase_name{xph});

%% ============= UNPACK PARAMETERS
xindex          = args.index{xph};% param index from study

xparam          = xindex.param;
xcond_name      = xparam.conds_names;
n_condition     = length(xcond_name);
dur_stim        = xparam.dur_stim;
dur_manip       = xparam.dur_manipulation;
n_tc_trs        = xparam.n_tc_trs;
xalpha          = 0.05;
line_w          = 2; 

condition_names = {'maintain','replace category','replace subcategory','target suppress','global clear'};

%*************** output basename
base_name       = args.analysis_basename;

%% ============= PLOT PARAMS
if strcmp(args.level, 'subcategory')
    for xcond = 1:n_condition
        
        if xcond==2, n_targ = 5; elseif xcond==3, n_targ = 4; else n_targ = 3; end
        
        xcolor{xcond}{1}      = [238, 20, 91]/255;% targ
        xcolor{xcond}{2}      = [160, 46, 83]/255;% nonsubtarg
        xcolor{xcond}{n_targ} = [144, 144, 144]/255;% non_target
        
        if (xcond==2) || (xcond==3)
                xcolor{xcond}{3} = [0, 188, 182]/255;% newtarg      
        end
        if xcond==2
            xcolor{xcond}{4} = [34, 127, 126]/255;% new_nonsubtarg
        end
    end
elseif strcmp(args.level,'category')
    for xcond = 1:n_condition
        
        if xcond==2, n_targ = 3; else n_targ = 2; end
        
        xcolor{xcond}{1}      = [238, 20, 91]/255;% targ
        xcolor{xcond}{n_targ} = [144, 144, 144]/255;% non_target
        
        if xcond==2
            xcolor{xcond}{2} = [0, 188, 182]/255;% newtarg  
        end
    end
end

on_stim         = args.shift_TRs + 1;
on_operation    = args.shift_TRs + 1 + dur_stim;
on_fixation     = on_operation + dur_manip;

x_tick          = 1:n_tc_trs;
x_ticklable     = 0:(n_tc_trs-1);
x_lim           = [x_tick(1) x_tick(end)];% + [-0.5 0.5];
y_lim           = [0 1.5];

%*************** legend
if strcmp(args.level,'subcategory')
    for i = [1 4 5];
        legend_names{i} = {'target','nontargets','baseline'};
    end
    
    legend_names{2} = {'target','nontargets','new target','new nontargets','baseline'};
    legend_names{3} = {'target','nontargets','new target','baseline'};
    
elseif strcmp(args.level,'category')
    for i = [1 3:5];
        legend_names{i} = {'target','baseline'};
    end
    
    legend_names{2} = {'target','new target','baseline'};
    
end

xcond_color{1} = [255, 0, 0]/255;
xcond_color{2} = [255, 127, 0]/255;
xcond_color{3} = [0, 255, 0]/255;
xcond_color{4} = [0, 0, 255]/255;
xcond_color{5} = [148, 0, 211]/255;

xbase_color    = [144, 144, 144]/255;
xonset_color   = 'r';

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%*************** CONCATENATED TIMECOURSE: RANDOM & FIXED
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%*************** working memory contents
% condition: 1_maintain,2_replace_category,3_replace_item,4_target_suppress,5_global_clear
% xtr: 12 presentation + 5:9 fixation = 21 max
% decode.timecourse.condition{xcond}.target{xtarg}.evidence{xtr}

%*************** working memory contents
%--------------- subcategory level
% 1. maintain            : {1} subtarget, {2} nonsubtarget(2), {3} mean(nontarget(6))
% 2. replace category    : {1} subtarget, {2} nonsubtarget(2), {3} new_subtarget, {4} new_nonsubtarget(2), {5} mean(nontarget(3))
% 3. replace subcategory : {1} subtarget, {2} nonsubtarget(1), {3} new_subtarget(1), {4} mean(nontarget(6))
% 4. target suppress     : {1} subtarget, {2} nonsubtarget(2), {3} mean(nontargets(6))
% 5. global clear        : {1} subtarget, {2} nonsubtarget(2), {3} mean(nontargets(6))

%% ============= TIMECOURSE
% decode.timecourse.condition{xcond}.target{xtarg}.evidence{xtr}

for xcond = 1:n_condition
    if strcmp(args.level,'subcategory')
        if xcond==2, n_targ = 5; elseif xcond==3, n_targ = 4; else n_targ = 3; end
    elseif strcmp(args.level,'category')
        if xcond==2, n_targ = 3; else n_targ = 2; end
    end

    for xtr = 1:n_tc_trs
        for xtarg = 1:n_targ
            xevidence = []; xevidence_random = [];
            
            %*************** fixed effect
            for xsub = args.g_sub
                xevidence = horzcat(xevidence, ...
                    grp_mvpaout{xsub}.decode.timecourse.condition{xcond}.target{xtarg}.evidence{xtr});
            end
            
            timecourse_fixed.evidence{xcond}{xtarg}{xtr} = xevidence;
            timecourse_fixed.mean{xcond}(xtarg, xtr)     = mean(xevidence);
            timecourse_fixed.std{xcond}(xtarg, xtr)      = std(xevidence);
            
            %*************** if tr is empty
            if isnan(timecourse_fixed.mean{xcond}(xtarg, xtr))
                timecourse_fixed.mean{xcond}(xtarg, xtr) = 0;
                timecourse_fixed.std{xcond}(xtarg, xtr)  = 0;
            end
            
            %*************** random effect
            for xsub = args.g_sub
                xevidence_random = horzcat(xevidence_random, ...
                    mean(grp_mvpaout{xsub}.decode.timecourse.condition{xcond}.target{xtarg}.evidence{xtr}));
            end
            
            timecourse_random.evicence{xcond}{xtarg}{xtr} = xevidence_random;
            timecourse_random.mean{xcond}(xtarg, xtr)     = mean(xevidence_random);
            timecourse_random.std{xcond}(xtarg, xtr)      = std(xevidence_random);
        end
    end
end

%% ============= STATS
for xcond = 1:n_condition
    if strcmp(args.level,'subcategory')
        if xcond==2, n_targ = 5; elseif xcond==3, n_targ = 4; else n_targ = 3; end
    elseif strcmp(args.level,'category')
        if xcond==2, n_targ = 3; else n_targ = 2; end
    end

    for xtr = 1:n_tc_trs

        timecourse_random.t_test{xcond}{xtr}.pvalue = zeros(n_targ, n_targ);
        
        for xcol = 1:n_targ
            for xrow = 1:n_targ
                xtarg_col = timecourse_random.evicence{xcond}{xcol}{xtr};%array_1
                xtarg_row = timecourse_random.evicence{xcond}{xrow}{xtr};%array_2
                
                %*************** ttest
                [xhypothesis, xpvalue, ~, xstats] = ttest(xtarg_col, xtarg_row, 'Alpha', 0.05);
                
                timecourse_random.t_test{xcond}{xtr}.pvalue(xrow, xcol)     = xpvalue;
                timecourse_random.t_test{xcond}{xtr}.hypothesis(xrow, xcol) = xhypothesis;
                timecourse_random.t_test{xcond}{xtr}.stats{xrow, xcol}      = xstats;
            end
        end
    end
end

%% ============= TIMECOURSE: TARGET-NONTARGER
%*************** working memory contents
% condition: 1_maintain,2_replace_category,3_replace_item,4_target_suppress,5_global_clear
% xtr: 12 presentation + 5:9 fixation = 21 max
% decode.timecourse.condition{xcond}.target{xtarg}.evidence{xtr}

%*************** working memory contents
%--------------- subcategory level
% 1. maintain            : {1} subtarget, {2} nonsubtarget(2), {3} mean(nontarget(6))
% 2. replace category    : {1} subtarget, {2} nonsubtarget(2), {3} new_subtarget, {4} new_nonsubtarget(2), {5} mean(nontarget(3))
% 3. replace subcategory : {1} subtarget, {2} nonsubtarget(1), {3} new_subtarget(1), {4} mean(nontarget(6))
% 4. target suppress     : {1} subtarget, {2} nonsubtarget(2), {3} mean(nontargets(6))
% 5. global clear        : {1} subtarget, {2} nonsubtarget(2), {3} mean(nontargets(6))
%--------------- category level
% 1. maintain            : {1} target, {2} mean(nontargets)
% 2. replace category    : {1} target, {2} new category, {3} nontargets
% 3. replace subcategory : {1} target, {2} mean(nontargets)
% 4. target suppress     : {1} target, {2} mean(nontargets)
% 5. global clear        : {1} target, {2} mean(nontargets)

for xcond = 1:n_condition
    if strcmp(args.level,'subcategory')
        if xcond==2, xbase = 5; elseif xcond==3, xbase = 4; else xbase = 3; end
    elseif strcmp(args.level,'category')
        if xcond==2, xbase = 3; else xbase = 2; end
    end
    
    xtarg = 1;

    for xtr = 1:n_tc_trs
        xevidence = timecourse_random.evicence{xcond}{xtarg}{xtr} -...
            timecourse_random.evicence{xcond}{xbase}{xtr};
        
        timecourse_target.evicence{xcond}{xtr}  = xevidence;
        timecourse_target.mean{xcond}(xtr)      = mean(xevidence);
        timecourse_target.std{xcond}(xtr)       = std(xevidence);
    end
end

%% ============= STATS
for xtr = 1:n_tc_trs
    
    timecourse_target.t_test{xtr}.pvalue = zeros(n_condition-1, n_condition-1);
    
    for xcol = 1:n_condition-1
        for xrow = 2:n_condition
            xtarg_col = timecourse_target.evicence{xcol}{xtr};%array_1
            xtarg_row = timecourse_target.evicence{xrow}{xtr};%array_2
            
            %*************** ttest
            [xhypothesis, xpvalue, ~, xstats] = ttest(xtarg_row, xtarg_col, 'Alpha', 0.05);
            
            timecourse_target.t_test{xtr}.pvalue(xrow, xcol)     = xpvalue;
            timecourse_target.t_test{xtr}.hypothesis(xrow, xcol) = xhypothesis;
            timecourse_target.t_test{xtr}.stats{xrow, xcol}      = xstats;
            
        end
    end
end

%% ============= SAVE
grp_timecourse.fixed  = timecourse_fixed; %#ok<*STRNU>
grp_timecourse.random = timecourse_random; 
grp_timecourse.target = timecourse_target; 

fname = fullfile(dirs.mvpa.group.out{xph}, sprintf('grp_wmem_timecourse_%s.mat', base_name));
save(fname, 'grp_timecourse');

%% ============= PLOTTING
LogRegs_fig = figure;
set(LogRegs_fig, 'Position', [0 0 1500 250])

for xcond = 1:n_condition
    if strcmp(args.level,'subcategory')
        if xcond==2, n_targ = 5; elseif xcond==3, n_targ = 4; else n_targ = 3; end
    elseif strcmp(args.level,'category')
        if xcond==2, n_targ = 3; else n_targ = 2; end
    end
    
    subplot(1, n_condition, xcond)
    
    %*************** mean lines
    for xtarg = 1:n_targ
        xmean = timecourse_fixed.mean{xcond}(xtarg, x_tick);
        xerr  = timecourse_fixed.std{xcond}(xtarg, x_tick);
        errorbar(x_tick, xmean, xerr, '-o', 'Color', xcolor{xcond}{xtarg});
        hold on;
    end
    
    %*************** legend
    legend([legend_names{xcond}])
    
    plot([on_stim on_stim], y_lim,'--k', 'Color', xbase_color)
    plot([on_operation on_operation], y_lim,'--k', 'Color', xbase_color)
    plot([on_fixation on_fixation], y_lim,'--k', 'Color', xbase_color)
    
    set(gca,'xlim', x_lim, 'ylim', y_lim);
    set(gca,'XTick', x_tick)
    set(gca,'XTickLabel', x_ticklable)
    
    title(sprintf('%s classifier: %s (N=%s)', args.level, xcond_name{xcond}, num2str(length(args.g_sub))));
    xlabel('Volume (tr)');
    ylabel('classifier evidence');
end

fig_fname = fullfile(dirs.mvpa.group.out{xph}, sprintf('grp_plot_wmem_timecourse_%s', base_name));

savefig(LogRegs_fig, sprintf('%s.fig', fig_fname));
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 15 2.5])
print('-djpeg', sprintf('%s.jpg',fig_fname), '-r100')
saveas(LogRegs_fig, sprintf('%s.jpg',fig_fname), 'jpg')

close(LogRegs_fig);

%% ============= CONTINUOUS PLOT: FIXED EFFECT
for xcond = 1:n_condition
    clear xmean xstd fity
    
    LogRegs_fig = figure;
    set(LogRegs_fig, 'Position', [0 0 1000 500])
    
    if strcmp(args.level,'subcategory')
        if xcond==2, n_targ = 5; elseif xcond==3, n_targ = 4; else n_targ = 3; end
    elseif strcmp(args.level,'category')
        if xcond==2, n_targ = 3; else n_targ = 2; end
    end
    
    %*************** mean line plots
    fitx = linspace(1, n_tc_trs, n_tc_trs*10);

    for xtarg = 1:n_targ
        xmean{xtarg} = timecourse_fixed.mean{xcond}(xtarg, x_tick);
        xstd{xtarg}  = timecourse_fixed.std{xcond}(xtarg, x_tick); %#ok<*AGROW>
        
        fity{xtarg}  = interp1(x_tick, xmean{xtarg}, fitx,'spline');
        
        plot(fitx, fity{xtarg}, '-','Color', xcolor{xcond}{xtarg}, 'LineWidth', 2); hold on;
    end
    
    %*************** legend
    legend([legend_names{xcond}])
    
    %*************** std error-bar filling
    for xtarg = 1:n_targ
        clear xerr
        xerr(1,:) = xmean{xtarg} - xstd{xtarg};
        xerr(2,:) = xmean{xtarg} + xstd{xtarg};
        
        for i = 1:2, fit_err(:,i) = interp1(x_tick, xerr(i,:), fitx,'spline'); end
        
        in_between = [fit_err(:,1)', fliplr(fit_err(:,2)')];
        h = fill([fitx fliplr(fitx)], in_between, xcolor{xcond}{xtarg});
        set(h,'facealpha', .2)
        
        plot(fitx, fit_err(:,1), '-', 'Color', xcolor{xcond}{xtarg})
        plot(fitx, fit_err(:,2), '-', 'Color', xcolor{xcond}{xtarg})
    end
    
    %*************** real onset lines
    plot([dur_stim dur_stim]+1, y_lim, '-', 'Color', xonset_color)
    plot([(dur_stim + dur_manip) (dur_stim + dur_manip)]+1, y_lim, '-', 'Color', xonset_color)
    
    %*************** shifted onset lines
    plot([on_stim on_stim], y_lim,'--', 'Color', xbase_color)
    plot([on_operation on_operation], y_lim,'--', 'Color', xbase_color)
    plot([on_fixation on_fixation], y_lim,'--', 'Color', xbase_color)
    
    set(gca,'xlim', x_lim, 'ylim', y_lim);
    set(gca,'XTick', x_tick)
    set(gca,'YTick', y_lim(1):0.1:y_lim(end))
    set(gca,'XTickLabel', x_ticklable)
    
    title(sprintf('%s classifier: %s (N=%s)', args.level, condition_names{xcond}, num2str(length(args.g_sub))));
    xlabel('Volume (tr)');
    ylabel('classifier evidence');
    
    h1 = text(1.5, 0.2, '(1) stim onset', 'Color', xonset_color);
    h2 = text(dur_stim + 1.5, 0.2, '(2) operation onset', 'Color', xonset_color);
    h3 = text(dur_stim + dur_manip + 1.5, 0.2, '(3) fixation onset', 'Color', xonset_color);
    set(h1,'Rotation', 90); set(h2,'Rotation', 90); set(h3,'Rotation', 90);
    
    text(on_stim - 1, 0.1, sprintf('(1) shift %str', num2str(args.shift_TRs)));
    text(on_operation - 1, 0.1, sprintf('(2) shift %str', num2str(args.shift_TRs)));
    text(on_fixation - 1, 0.1, sprintf('(3) shift %str', num2str(args.shift_TRs)));
    
    %% *************** save fig
    fig_fname = fullfile(dirs.mvpa.group.out{xph}, ...
        sprintf('grp_plot_continous_wn_timecourse_fixed_%s_%s', ...
        xcond_name{xcond}, base_name));
    
    savefig(LogRegs_fig, sprintf('%s.fig', fig_fname));
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 10 5])
    print('-djpeg', sprintf('%s.jpg',fig_fname), '-r100')
    saveas(LogRegs_fig, sprintf('%s.jpg',fig_fname), 'jpg')
    
    close(LogRegs_fig);
    
end

%% ============= CONTINUOUS PLOT: RANDOM EFFECT
y_lim   = [0 1.5];
figrect = [0 0 1700 1000];

for xcond = 1:n_condition
    clear xmean xstd fity
    
    LogRegs_fig = figure;
    set(LogRegs_fig, 'Position', figrect)
    
    if strcmp(args.level,'subcategory')
        if xcond==2, n_targ = 5; elseif xcond==3, n_targ = 4; else n_targ = 3; end
    elseif strcmp(args.level,'category')
        if xcond==2, n_targ = 3; else n_targ = 2; end
    end
    
    %*************** mean line plots
    fitx = linspace(1, n_tc_trs, n_tc_trs*10);
    
    for xtarg = 1:n_targ
        xmean{xtarg} = timecourse_random.mean{xcond}(xtarg, x_tick);
        xstd{xtarg}  = timecourse_random.std{xcond}(xtarg, x_tick); %#ok<*AGROW>
        
        fity{xtarg}  = interp1(x_tick, xmean{xtarg}, fitx,'spline');
        
        plot(fitx, fity{xtarg}, '-','Color', xcolor{xcond}{xtarg}, 'LineWidth', 2); hold on;
    end
    
    %*************** legend
    legend([legend_names{xcond}], 'Location', 'SouthEast')
    
    %*************** stats stars
    n = 0;
    for xcol = 1:n_targ
        for xrow = 1:n_targ
            if xcol < xrow
                xheight = 1.4 - (n * 0.03);
                text(-4, xheight, sprintf('%s vs. %s', legend_names{xcond}{xcol}, ...
                    legend_names{xcond}{xrow}));
                plot(x_lim, [xheight xheight], ':', 'Color', [0.75 0.75 0.75])
                
                for xtr = 1:n_tc_trs
                    xpvalue = timecourse_random.t_test{xcond}{xtr}.pvalue(xrow, xcol);
                    
                    if xpvalue <= xalpha
                        text(xtr-0.1, xheight, '*', 'FontSize', 20);
                    end
                    
                    plot([xtr xtr], y_lim, ':', 'Color', [0.75 0.75 0.75])
                    
                end
                
                n = n + 1;
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
        h = fill([fitx fliplr(fitx)], in_between, xcolor{xcond}{xtarg});
        set(h,'facealpha', .4)
        
        plot(fitx, fit_err(:,1), '-', 'Color', xcolor{xcond}{xtarg})
        plot(fitx, fit_err(:,2), '-', 'Color', xcolor{xcond}{xtarg})
    end
    
    %*************** real onset lines
    plot([dur_stim dur_stim]+1, y_lim, '-', 'Color', xonset_color)
    plot([(dur_stim + dur_manip) (dur_stim + dur_manip)]+1, y_lim, '-', 'Color', xonset_color)
    
    %*************** shifted onset lines
    plot([on_stim on_stim], y_lim,'--', 'Color', xbase_color)
    plot([on_operation on_operation], y_lim,'--', 'Color', xbase_color)
    plot([on_fixation on_fixation], y_lim,'--', 'Color', xbase_color)
    
    set(gca,'xlim', x_lim, 'ylim', y_lim);
    set(gca,'XTick', x_tick)
    set(gca,'YTick', y_lim(1):0.1:y_lim(end))
    set(gca,'XTickLabel', x_ticklable)
    
    title(sprintf('%s classifier: %s (N=%s)', args.level, condition_names{xcond}, num2str(length(args.g_sub))));
    xlabel('Volume (tr)');
    ylabel('classifier evidence');
    
    h1 = text(1.5, 0.2, '(1) stim onset', 'Color', xonset_color);
    h2 = text(dur_stim + 1.5, 0.2, '(2) operation onset', 'Color', xonset_color);
    h3 = text(dur_stim + dur_manip + 1.5, 0.2, '(3) fixation onset', 'Color', xonset_color);
    set(h1,'Rotation', 90); set(h2,'Rotation', 90); set(h3,'Rotation', 90);
    
    text(on_stim - 1, 0.1, sprintf('(1) shift %str', num2str(args.shift_TRs)));
    text(on_operation - 1, 0.1, sprintf('(2) shift %str', num2str(args.shift_TRs)));
    text(on_fixation - 1, 0.1, sprintf('(3) shift %str', num2str(args.shift_TRs)));
    
    %% *************** save fig
    fig_fname = fullfile(dirs.mvpa.group.out{xph}, ...
        sprintf('grp_plot_continous_wn_timecourse_random_%s_%s', ...
        xcond_name{xcond}, base_name));
    
    savefig(LogRegs_fig, sprintf('%s.fig', fig_fname));
    set(gcf,'PaperUnits','inches','PaperPosition', figrect/100)
    print('-djpeg', sprintf('%s.jpg',fig_fname), '-r100')
    saveas(LogRegs_fig, sprintf('%s.jpg',fig_fname), 'jpg')
    
    close(LogRegs_fig);
end

%% ============= CONTINUOUS PLOT: TARGET RANDOM EFFECT 
%*************** MAINTAIN VS. SUPPRESS VS. CLEAR 
%*************** evidence = target - baseline 
std_fill  = 0;
sel_conds{1} = 1:n_condition;
sel_conds{2} = [1 4 5];
sel_conds{3} = 2:3;

for xsel = 1:length(sel_conds)
    xsel_conds = sel_conds{xsel};
    
    y_lim   = [-0.3 0.8];
    figrect = [0 0 2000 1000];
    
    LogRegs_fig = figure;
    set(LogRegs_fig, 'Position', figrect)
    
    clear xmean xstd fity
    
    %*************** mean line plots
    fitx = linspace(1, n_tc_trs, n_tc_trs*10);
    
    for xcond = xsel_conds
        xmean{xcond} = timecourse_target.mean{xcond}(x_tick);
        xstd{xcond}  = timecourse_target.std{xcond}(x_tick); %#ok<*AGROW>
        
        fity{xcond}  = interp1(x_tick, xmean{xcond}, fitx,'spline');
        
        plot(fitx, fity{xcond}, '-','Color', xcond_color{xcond}, 'LineWidth', line_w); hold on;
    end
    
    %*************** legend
    legend(condition_names{xsel_conds}, 'Location', 'SouthEast')
    
    %*************** stats stars
    
    n = 0;
    for i = 1:(length(xsel_conds)-1)
        xcol = xsel_conds(i);
        for j = i+1:length(xsel_conds)
            xrow = xsel_conds(j);
            xheight = y_lim(2) - 0.05 - (n * 0.03);
            text(-4, xheight, sprintf('%s vs. %s', condition_names{xcol}, ...
                condition_names{xrow}));
            plot(x_lim, [xheight xheight], ':', 'Color', [0.75 0.75 0.75])
            
            for xtr = 1:n_tc_trs
                xpvalue = timecourse_target.t_test{xtr}.pvalue(xrow, xcol);
                
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
        for xcond = xsel_conds
            clear xerr
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
    
    %*************** real onset lines
    plot([dur_stim dur_stim]+1, y_lim, '-','Color', xonset_color)
    plot([(dur_stim + dur_manip) (dur_stim + dur_manip)]+1, y_lim, '-', 'Color', xonset_color)
    
    %*************** shifted onset lines
    plot([on_stim on_stim], y_lim,'--', 'Color', xbase_color)
    plot([on_operation on_operation], y_lim,'--', 'Color', xbase_color)
    plot([on_fixation on_fixation], y_lim,'--', 'Color', xbase_color)
    
    set(gca,'xlim', x_lim, 'ylim', y_lim);
    set(gca,'XTick', x_tick)
    set(gca,'YTick', y_lim(1):0.1:y_lim(end))
    set(gca,'XTickLabel', x_ticklable)
    
    title(sprintf('%s classifier: target-baseline (N=%s)', args.level,num2str(length(args.g_sub))));
    xlabel('Volume (tr)');
    ylabel('classifier evidence (target-baseline)');
    
    h1 = text(1.5, y_lim(1), '(1) stim onset', 'Color', xonset_color);
    h2 = text(dur_stim + 1.5, y_lim(1), '(2) operation onset', 'Color', xonset_color);
    h3 = text(dur_stim + dur_manip + 1.5, y_lim(1), '(3) fixation onset', 'Color', xonset_color);
    set(h1,'Rotation', 90); set(h2,'Rotation', 90); set(h3,'Rotation', 90);
    
    text(on_stim - 1, y_lim(1)+0.1, sprintf('(1) shift %str', num2str(args.shift_TRs)));
    text(on_operation - 1, y_lim(1)+0.1, sprintf('(2) shift %str', num2str(args.shift_TRs)));
    text(on_fixation - 1, y_lim(1)+0.1, sprintf('(3) shift %str', num2str(args.shift_TRs)));
    
    %% *************** save fig
    if std_fill
        fig_fname = fullfile(dirs.mvpa.group.out{xph}, ...
            sprintf('grp_plot_continous_wn_timecourse_sel_target_%d_%s', xsel, base_name));
    else
        fig_fname = fullfile(dirs.mvpa.group.out{xph}, ...
            sprintf('grp_plot_continous_wn_timecourse_iine_sel_target_%d_%s', xsel, base_name));
    end
    
    savefig(LogRegs_fig, sprintf('%s.fig', fig_fname));
    set(gcf,'PaperUnits','inches','PaperPosition', figrect/100)
    print('-djpeg', sprintf('%s.jpg',fig_fname), '-r100')
    saveas(LogRegs_fig, sprintf('%s.jpg',fig_fname), 'jpg')
    
    close(LogRegs_fig);
end

end