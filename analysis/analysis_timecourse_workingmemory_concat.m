function[] = analysis_timecourse_workingmemory_concat(args, mvpaout, dirs)

%---------------------------------------------------------------------------
%*************** timecourse of the working memory (category|subcategory)
%---------------------------------------------------------------------------

xph = args.xphase;
fprintf('\n(+) timecourse of working memory: %s\n', args.phase_name{xph});

%% ============= UNPACK PARAMETERS
decode          = mvpaout.decode;
xindex          = args.index{xph};% param index from study

xparam          = xindex.param;
xcond_name      = xparam.conds_names;
n_condition     = length(xcond_name);
dur_stim        = xparam.dur_stim;
dur_manip       = xparam.dur_manipulation;
n_tc_trs        = xparam.n_tc_trs;

condition_names = {'maintain','replace category','replace subcategory','target suppress','global clear'};

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%*************** CONCATENATED
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

%*************** output basename
basename        = args.analysis_basename;

%% ============= PLOT PARAMS
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

on_stim         = args.shift_TRs + 1;
on_operation    = args.shift_TRs + 1 + dur_stim;
on_fixation     = on_operation + dur_manip;

x_tick          = 1:n_tc_trs;
x_ticklable     = 0:(n_tc_trs-1);
x_lim           = [x_tick(1) x_tick(end)];% + [-0.5 0.5];
y_lim           = [0 1.5];

%% ============= TIMECOURSE
% decode.timecourse.condition{xcond}.target{xtarg}.evidence{xtr}

for xcond = 1:n_condition
    if xcond==2, n_targ = 5; elseif xcond==3, n_targ = 4; else n_targ = 3; end
    
    for xtr = 1:n_tc_trs
        for xtarg = 1:n_targ
            timecourse.mean{xcond}(xtarg, xtr) = ...
                mean(decode.timecourse.condition{xcond}.target{xtarg}.evidence{xtr});
            timecourse.std{xcond}(xtarg, xtr) = ...
                std(decode.timecourse.condition{xcond}.target{xtarg}.evidence{xtr});
            
            %*************** if tr is empty
            if isnan(timecourse.mean{xcond}(xtarg, xtr))
                timecourse.mean{xcond}(xtarg, xtr) = 0;
                timecourse.std{xcond}(xtarg, xtr)  = 0;
            end
        end
    end
end

%% ============= SAVE
fname = fullfile(dirs.mvpa.parse{xph}, sprintf('wmem_timecourse_%s.mat', basename));
save(fname, 'timecourse');

%% ============= PLOTTING
LogRegs_fig = figure;
set(LogRegs_fig, 'Position', [0 0 1500 250])

for xcond = 1:n_condition
    if xcond==2, n_targ = 5; elseif xcond==3, n_targ = 4; else n_targ = 3; end
    
    subplot(1, n_condition, xcond)
    
    %*************** mean lines
    for xtarg = 1:n_targ
        xmean = timecourse.mean{xcond}(xtarg, x_tick);
        xerr  = timecourse.std{xcond}(xtarg, x_tick);
        errorbar(x_tick, xmean, xerr, '-o', 'Color', xcolor{xcond}{xtarg});
        hold on;
    end
    
    %*************** legend
    if xcond==2
        legend('target','nontargets','new target','new nontargets','baseline');
    elseif xcond==3
        legend('target','nontargets','new target','baseline');
    else
        legend('target','nontargets','baseline');
    end
    
    plot([on_stim on_stim], y_lim,'--k', 'Color', xcolor{xcond}{n_targ})
    plot([on_operation on_operation], y_lim,'--k', 'Color', xcolor{xcond}{n_targ})
    plot([on_fixation on_fixation], y_lim,'--k', 'Color', xcolor{xcond}{n_targ})
    
    set(gca,'xlim', x_lim, 'ylim', y_lim);
    set(gca,'XTick', x_tick)
    set(gca,'XTickLabel', x_ticklable)
    
    title(sprintf('%s classifier: %s (%s)', args.level, xcond_name{xcond}, flip(strtok(flip(args.subject_id),'_'))));
    xlabel('Volume (tr)');
    ylabel('classifier evidence');
end

fig_fname = fullfile(dirs.mvpa.parse{xph}, sprintf('plot_wmem_timecourse_%s', basename));

savefig(LogRegs_fig, sprintf('%s.fig', fig_fname));
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 15 2.5])
print('-djpeg', sprintf('%s.jpg',fig_fname), '-r100')
saveas(LogRegs_fig, sprintf('%s.jpg',fig_fname), 'jpg')

close(LogRegs_fig);

%% ============= CONTINUOUS PLOT
for xcond = 1:n_condition
    LogRegs_fig = figure;
    set(LogRegs_fig, 'Position', [0 0 1000 500])
    
    clear xmean xstd fity
    
    if xcond==2, n_targ = 5; elseif xcond==3, n_targ = 4; else n_targ = 3; end
    
    %*************** mean line plots
    fitx = linspace(1, n_tc_trs, n_tc_trs*10);
    
    for xtarg = 1:n_targ
        xmean{xtarg} = timecourse.mean{xcond}(xtarg, x_tick);
        xstd{xtarg}  = timecourse.std{xcond}(xtarg, x_tick); %#ok<*AGROW>
        
        fity{xtarg}  = interp1(x_tick, xmean{xtarg}, fitx,'spline');
        
        plot(fitx, fity{xtarg}, '-','Color', xcolor{xcond}{xtarg}, 'LineWidth', 2); hold on;
    end
    
    %*************** legend
    if xcond==2
        legend('target','nontargets','new target','new nontargets','baseline');
    elseif xcond==3
        legend('target','nontargets','new target','baseline');
    else
        legend('target','nontargets','baseline');
    end
    
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
    plot([dur_stim dur_stim]+1, y_lim, '-', 'Color', xcolor{xcond}{1})
    plot([(dur_stim + dur_manip) (dur_stim + dur_manip)]+1, y_lim, '-', 'Color', xcolor{xcond}{1})
    
    %*************** shifted onset lines
    plot([on_stim on_stim], y_lim,'--', 'Color', xcolor{xcond}{n_targ})
    plot([on_operation on_operation], y_lim,'--', 'Color', xcolor{xcond}{n_targ})
    plot([on_fixation on_fixation], y_lim,'--', 'Color', xcolor{xcond}{n_targ})
    
    set(gca,'xlim', x_lim, 'ylim', y_lim);
    set(gca,'XTick', x_tick)
    set(gca,'YTick', y_lim(1):0.1:y_lim(end))
    set(gca,'XTickLabel', x_ticklable)
    
    title(sprintf('%s classifier: %s (%s)', args.level, condition_names{xcond}, flip(strtok(flip(args.subject_id),'_'))));
    xlabel('Volume (tr)');
    ylabel('classifier evidence');
    
    h1 = text(1.5, 0.2, '(1) stim onset', 'Color', xcolor{xcond}{1});
    h2 = text(dur_stim + 1.5, 0.2, '(2) operation onset', 'Color', xcolor{xcond}{1});
    h3 = text(dur_stim + dur_manip + 1.5, 0.2, '(3) fixation onset', 'Color', xcolor{xcond}{1});
    set(h1,'Rotation', 90); set(h2,'Rotation', 90); set(h3,'Rotation', 90);
    
    text(on_stim - 1, 0.1, sprintf('(1) shift %str', num2str(args.shift_TRs)));
    text(on_operation - 1, 0.1, sprintf('(2) shift %str', num2str(args.shift_TRs)));
    text(on_fixation - 1, 0.1, sprintf('(3) shift %str', num2str(args.shift_TRs)));
    
    %% ============= SAVE FIGURE
    fig_fname = fullfile(dirs.mvpa.parse{xph}, ...
        sprintf('plot_continous_wn_timecourse_%s_%s', ...
        xcond_name{xcond}, basename));
    
    savefig(LogRegs_fig, sprintf('%s.fig', fig_fname));
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 10 5])
    print('-djpeg', sprintf('%s.jpg',fig_fname), '-r100')
    saveas(LogRegs_fig, sprintf('%s.jpg',fig_fname), 'jpg')
    
    close(LogRegs_fig);
    
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%*************** SEPARATED CATEGORY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%*************** working memory contents
% condition: 1_maintain,2_replace_category,3_replace_item,4_target_suppress,5_global_clear
% xtr: 12 presentation + 5:9 fixation = 21 max
% decode.timecourse_each.condition{xcond}.target{xtarg}.evidence{xtr}

%--------------- subcategory level
% 1. maintain            : {1} subtarget, {2} nonsubtarget(2)
% 2. replace category    : {1} subtarget, {2} nonsubtarget(2), {3} new_subtarget(1), {4} new_nonsubtarget(2)
% 3. replace subcategory : {1} subtarget, {2} nonsubtarget(1), {3} new_subtarget(1)
% 4. target suppress     : {1} subtarget, {2} nonsubtarget(2) 
% 5. global clear        : {1} subtarget, {2} nonsubtarget(2)

%*************** output basename
basename        = sprintf('each_%s', args.analysis_basename);

%% ============= PLOT PARAMS
clear xcolor

for xcond = 1:n_condition
    
    xcolor{xcond}{1}      = [238, 20, 91]/255;% targ
    xcolor{xcond}{2}      = [160, 46, 83]/255;% nonsubtarg
    
    if (xcond==2) || (xcond==3)
        xcolor{xcond}{3} = [0, 188, 182]/255;% newtarg
    end
    if xcond==2
        xcolor{xcond}{4} = [34, 127, 126]/255;% new_nonsubtarg
    end
end

on_stim         = args.shift_TRs + 1;
on_operation    = args.shift_TRs + 1 + dur_stim;
on_fixation     = on_operation + dur_manip;

x_tick          = 1:n_tc_trs;
x_ticklable     = 0:(n_tc_trs-1);
x_lim           = [x_tick(1) x_tick(end)];% + [-0.5 0.5];
y_lim           = [0 1.5];

%% ============= TIMECOURSE
% decode.timecourse.condition{xcond}.target{xtarg}.evidence{xtr}

for xcond = 1:n_condition
    if xcond==2, n_targ = 4; elseif xcond==3, n_targ = 3; else n_targ = 2; end
    
    for xtr = 1:n_tc_trs
        for xtarg = 1:n_targ
            timecourse.mean{xcond}(xtarg, xtr) = ...
                mean(decode.timecourse_each.condition{xcond}.target{xtarg}.evidence{xtr});
            timecourse.std{xcond}(xtarg, xtr) = ...
                std(decode.timecourse_each.condition{xcond}.target{xtarg}.evidence{xtr});
            
            %*************** if tr is empty
            if isnan(timecourse.mean{xcond}(xtarg, xtr))
                timecourse.mean{xcond}(xtarg, xtr) = 0;
                timecourse.std{xcond}(xtarg, xtr)  = 0;
            end
        end
    end
end

%% ============= SAVE
fname = fullfile(dirs.mvpa.parse{xph}, sprintf('wmem_timecourse_%s.mat', basename));
save(fname, 'timecourse');

%% ============= PLOTTING
LogRegs_fig = figure;
set(LogRegs_fig, 'Position', [0 0 1500 250])

for xcond = 1:n_condition
    if xcond==2, n_targ = 4; elseif xcond==3, n_targ = 3; else n_targ = 2; end
    
    subplot(1, n_condition, xcond)
    
    %*************** mean lines
    for xtarg = 1:n_targ
        xmean = timecourse.mean{xcond}(xtarg, x_tick);
        xerr  = timecourse.std{xcond}(xtarg, x_tick);
        errorbar(x_tick, xmean, xerr, '-o', 'Color', xcolor{xcond}{xtarg});
        hold on;
    end
    
    %*************** legend
    if xcond==2
        legend('target','nontargets','new target','new nontargets','baseline');
    elseif xcond==3
        legend('target','nontargets','new target','baseline');
    else
        legend('target','nontargets','baseline');
    end
    
    plot([on_stim on_stim], y_lim,'--k', 'Color', xcolor{xcond}{n_targ})
    plot([on_operation on_operation], y_lim,'--k', 'Color', xcolor{xcond}{n_targ})
    plot([on_fixation on_fixation], y_lim,'--k', 'Color', xcolor{xcond}{n_targ})
    
    set(gca,'xlim', x_lim, 'ylim', y_lim);
    set(gca,'XTick', x_tick)
    set(gca,'XTickLabel', x_ticklable)
    
    title(sprintf('%s classifier: %s (%s)', args.level, xcond_name{xcond}, flip(strtok(flip(args.subject_id),'_'))));
    xlabel('Volume (tr)');
    ylabel('classifier evidence');
end

fig_fname = fullfile(dirs.mvpa.parse{xph}, sprintf('plot_wmem_timecourse_%s', basename));

savefig(LogRegs_fig, sprintf('%s.fig', fig_fname));
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 15 2.5])
print('-djpeg', sprintf('%s.jpg',fig_fname), '-r100')
saveas(LogRegs_fig, sprintf('%s.jpg',fig_fname), 'jpg')

close(LogRegs_fig);

%% ============= CONTINUOUS PLOT
for xcond = 1:n_condition
    LogRegs_fig = figure;
    set(LogRegs_fig, 'Position', [0 0 1000 500])
    
    clear xmean xstd fity
    
    if xcond==2, n_targ = 4; elseif xcond==3, n_targ = 3; else n_targ = 2; end
    
    %*************** mean line plots
    fitx = linspace(1, n_tc_trs, n_tc_trs*10);
    
    for xtarg = 1:n_targ
        xmean{xtarg} = timecourse.mean{xcond}(xtarg, x_tick);
        xstd{xtarg}  = timecourse.std{xcond}(xtarg, x_tick); %#ok<*AGROW>
        
        fity{xtarg}  = interp1(x_tick, xmean{xtarg}, fitx,'spline');
        
        plot(fitx, fity{xtarg}, '-','Color', xcolor{xcond}{xtarg}, 'LineWidth', 2); hold on;
    end
    
    %*************** legend
    if xcond==2
        legend('target','nontargets','new target','new nontargets','baseline');
    elseif xcond==3
        legend('target','nontargets','new target','baseline');
    else
        legend('target','nontargets','baseline');
    end
    
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
    plot([dur_stim dur_stim]+1, y_lim, '-', 'Color', xcolor{xcond}{1})
    plot([(dur_stim + dur_manip) (dur_stim + dur_manip)]+1, y_lim, '-', 'Color', xcolor{xcond}{1})
    
    %*************** shifted onset lines
    plot([on_stim on_stim], y_lim,'--', 'Color', xcolor{xcond}{n_targ})
    plot([on_operation on_operation], y_lim,'--', 'Color', xcolor{xcond}{n_targ})
    plot([on_fixation on_fixation], y_lim,'--', 'Color', xcolor{xcond}{n_targ})
    
    set(gca,'xlim', x_lim, 'ylim', y_lim);
    set(gca,'XTick', x_tick)
    set(gca,'YTick', y_lim(1):0.1:y_lim(end))
    set(gca,'XTickLabel', x_ticklable)
    
    title(sprintf('%s classifier: %s (%s)', args.level, condition_names{xcond}, flip(strtok(flip(args.subject_id),'_'))));
    xlabel('Volume (tr)');
    ylabel('classifier evidence');
    
    h1 = text(1.5, 0.2, '(1) stim onset', 'Color', xcolor{xcond}{1});
    h2 = text(dur_stim + 1.5, 0.2, '(2) operation onset', 'Color', xcolor{xcond}{1});
    h3 = text(dur_stim + dur_manip + 1.5, 0.2, '(3) fixation onset', 'Color', xcolor{xcond}{1});
    set(h1,'Rotation', 90); set(h2,'Rotation', 90); set(h3,'Rotation', 90);
    
    text(on_stim - 1, 0.1, sprintf('(1) shift %str', num2str(args.shift_TRs)));
    text(on_operation - 1, 0.1, sprintf('(2) shift %str', num2str(args.shift_TRs)));
    text(on_fixation - 1, 0.1, sprintf('(3) shift %str', num2str(args.shift_TRs)));
    
    %% ============= SAVE FIGURE
    fig_fname = fullfile(dirs.mvpa.parse{xph}, ...
        sprintf('plot_continous_wn_timecourse_%s_%s', ...
        xcond_name{xcond}, basename));
    
    savefig(LogRegs_fig, sprintf('%s.fig', fig_fname));
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 10 5])
    print('-djpeg', sprintf('%s.jpg',fig_fname), '-r100')
    saveas(LogRegs_fig, sprintf('%s.jpg',fig_fname), 'jpg')
    
    close(LogRegs_fig);
    
end

end