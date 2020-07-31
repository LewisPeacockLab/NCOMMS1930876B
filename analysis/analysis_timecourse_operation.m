function[] = analysis_timecourse_operation(args, mvpaout, dirs)

%---------------------------------------------------------------------------
%*************** timecourse of operation during study
%---------------------------------------------------------------------------

xph = args.xphase;
fprintf('\n(+) timecourse of operation: %s\n', args.phase_name{xph});

%% ============= UNPACK PARAMETERS
decode           = mvpaout.decode;
xindex           = args.index{xph};% param index from study

xparam           = xindex.param;
xcond_name       = xparam.conds_names;
n_condition      = length(xcond_name);

dur_stim         = xparam.dur_stim;% 6 TR 
dur_manipulation = xparam.dur_manipulation;% 6 TR
it_trs           = args.tc_tr;%timecourse trs

condition_names  = {'maintain','replace category',...
    'replace subcategory','target suppress','global clear'};

%*************** output basename
basename         = sprintf('%s_%s', args.analysis_basename, num2str(args.xpenalty));

%% ============= PLOT PARAMS
xcolor{1}  = [238, 20, 91]/255;% targ
xcolor{2}  = [144, 144, 144]/255;% non_target

x_tick     = 1:it_trs;
x_lim      = [x_tick(1) x_tick(end)-1];% + [-0.5 0.5];
y_lim      = [0 1.5];
n_ticks    = x_tick(end);
line_w     = 2;
line_w_err = 0.1;

%*************** separated operations
xcond_color{1} = [255, 0, 0]/255;
xcond_color{2} = [255, 127, 0]/255;
xcond_color{3} = [0, 255, 0]/255;
xcond_color{4} = [0, 0, 255]/255;
xcond_color{5} = [148, 0, 211]/255;

%% ============= TIMECOURSE STRUCTURE
%*************** working memory operation
%*************** timecourse/timecourse_cat(concatenating runs)
%*************** timecourse{1}: target
%*************** timecourse{2}: nontarget baseline
for xcond = 1:n_condition
    for xtr = 1:it_trs
        for xtarg = 1:2
            timecourse{xcond}.mean{xtarg}(xtr) = ...
                mean(decode.timecourse.operation{xcond}.target{xtarg}.evidence{xtr}); %#ok<*AGROW>
            timecourse{xcond}.std{xtarg}(xtr)  = ...
                std(decode.timecourse.operation{xcond}.target{xtarg}.evidence{xtr});
            
            if isnan(timecourse{xcond}.mean{1}(xtr))
                timecourse{xcond}.mean{xtarg}(xtr) = 0;
                timecourse{xcond}.std{xtarg}(xtr)  = 0;
            end
        end
        
    end
end

%*************** all working memory operation
%*************** decode.timecourse.operation{xcond}.decoded_operation{xoper}.evidence{xtr}

for xcond = 1:n_condition
    for xtr = 1:it_trs
        for xtarg = 1:n_condition
            all_timecourse{xcond}.mean{xtarg}(xtr) = ...
                mean(decode.timecourse.operation{xcond}.decoded_operation{xtarg}.evidence{xtr}); %#ok<*AGROW>
            all_timecourse{xcond}.std{xtarg}(xtr)  = ...
                std(decode.timecourse.operation{xcond}.decoded_operation{xtarg}.evidence{xtr});
            
            if isnan(all_timecourse{xcond}.mean{1}(xtr))
                all_timecourse{xcond}.mean{xtarg}(xtr) = 0;
                all_timecourse{xcond}.std{xtarg}(xtr)  = 0;
            end
        end
        
    end
end

%% ============= SAVE
fname = fullfile(dirs.mvpa.parse{xph}, sprintf('operation_timecourse_%s.mat', basename));
save(fname, 'timecourse','all_timecourse');

%% ============= PLOTTING
LogRegs_fig = figure;
set(LogRegs_fig, 'Position', [0 0 1500 250])

for xcond = 1:n_condition
    subplot(1,5, xcond)
    errorbar(x_tick, timecourse{xcond}.mean{1}, timecourse{xcond}.std{1})
    
    set(gca,'xlim', x_lim, 'ylim', y_lim);
    set(gca,'XTick', x_tick)
    
    title(sprintf('Timecourse %s', xcond_name{xcond}));
    xlabel('Volume (tr)');
    ylabel('classifier evidence');
end

fig_fname = fullfile(dirs.mvpa.parse{xph}, sprintf('plot_timecourse_%s', basename));

savefig(LogRegs_fig, sprintf('%s.fig', fig_fname));
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 15 2.5])
print('-djpeg', sprintf('%s.jpg',fig_fname), '-r100')
saveas(LogRegs_fig, sprintf('%s.jpg',fig_fname), 'jpg')

close(LogRegs_fig);

%% ============= CONTINUOUS PLOT
n_targ = 2;%target/nontarget

for xcond = 1:n_condition
    LogRegs_fig = figure;
    if strcmp(args.rest,'nonrest'), fig_rect = [0 0 500 500];
    else                            fig_rect = [0 0 850 500]; end
    
    set(LogRegs_fig, 'Position', fig_rect)
    
    fitx = linspace(1, n_ticks, n_ticks*10);
    
    for xtarg = 1:n_targ
        %*************** continuous mean
        fity{xtarg} = interp1(x_tick, timecourse{xcond}.mean{xtarg}, fitx, 'spline');
        
        %*************** continuous errorbar
        xerr{xtarg}(1,:) = timecourse{xcond}.mean{xtarg} - timecourse{xcond}.std{xtarg};
        xerr{xtarg}(2,:) = timecourse{xcond}.mean{xtarg} + timecourse{xcond}.std{xtarg};
        
        for i = 1:2 
            fit_err{xtarg}(:,i) = interp1(x_tick, xerr{xtarg}(i,:), fitx,'spline'); 
        end
    end
    
    %*************** mean plotting
    for xtarg = 1:n_targ
        plot(fitx, fity{xtarg}, '-', 'Color', xcolor{xtarg}, 'LineWidth', line_w); hold on;
    end
    
    %*************** legend
    legend('target','baseline');
    
    %*************** errorbar plotting
    for xtarg = 1:n_targ
        in_between = [fit_err{xtarg}(:,1)', fliplr(fit_err{xtarg}(:,2)')];
        h = fill([fitx fliplr(fitx)], in_between, xcolor{xtarg});
        set(h,'facealpha', .2)
        
        for i = 1:2
            plot(fitx, fit_err{xtarg}(:,i), '-', 'Color', xcolor{xtarg},'LineWidth',line_w_err); 
        end
    end
    
    set(gca,'xlim', x_lim, 'ylim', y_lim);
    set(gca,'XTick', x_tick)
    set(gca,'YTick', y_lim(1):0.1:y_lim(end))
    
    title(sprintf('Timecourse: %s', condition_names{xcond}));
    xlabel('Volume (tr)');
    ylabel('classifier evidence');
    
    %*************** onset lines
    plot([dur_stim dur_stim], y_lim,'--', 'Color', xcolor{2})
    plot([dur_stim+dur_manipulation, dur_stim+dur_manipulation], y_lim,'--', 'Color', xcolor{2})
    
    text(1, 0.05, 'stim onset');
    text(dur_stim, 0.05, 'operation onset');
    text(dur_stim + dur_manipulation, 0.05, 'fixation onset');
    
    %% ============= SAVE FIGURE
    fig_fname = fullfile(dirs.mvpa.parse{xph}, ...
        sprintf('plot_continous_operation_timecourse_%s_%s', ...
        xcond_name{xcond}, basename));
    
    savefig(LogRegs_fig, sprintf('%s.fig', fig_fname));
    set(gcf,'PaperUnits','inches','PaperPosition',fig_rect/100)
    print('-djpeg', sprintf('%s.jpg',fig_fname), '-r100')
    saveas(LogRegs_fig, sprintf('%s.jpg',fig_fname), 'jpg')
    
    close(LogRegs_fig);
    
end

%% ============= ALL OPERATION: CONTINUOUS PLOT
n_targ = n_condition;

for xcond = 1:n_condition
    LogRegs_fig = figure;
    if strcmp(args.rest,'nonrest'), fig_rect = [0 0 500 500];
    else                            fig_rect = [0 0 850 500]; end
    
    set(LogRegs_fig, 'Position', fig_rect)
    
    fitx = linspace(1, n_ticks, n_ticks*10);
    
    for xtarg = 1:n_targ
        %*************** continuous mean
        fity{xtarg} = interp1(x_tick, all_timecourse{xcond}.mean{xtarg}, fitx,'spline');
        
        %*************** continuous errorbar
        xerr{xtarg}(1,:) = all_timecourse{xcond}.mean{xtarg} - all_timecourse{xcond}.std{xtarg};
        xerr{xtarg}(2,:) = all_timecourse{xcond}.mean{xtarg} + all_timecourse{xcond}.std{xtarg};
        
        for i = 1:2, fit_err{xtarg}(:,i) = interp1(x_tick, xerr{xtarg}(i,:), fitx,'spline'); end
    end
    
    %*************** mean plotting
    for xtarg = 1:n_targ
        if xtarg==xcond
            plot(fitx, fity{xtarg}, '-', 'Color', xcond_color{xtarg}, 'LineWidth', line_w + 1); hold on;
        else
            plot(fitx, fity{xtarg}, '-', 'Color', xcond_color{xtarg}, 'LineWidth', line_w); hold on;
        end
    end
    
    %*************** legend
    legend(condition_names);
    
    %*************** errorbar plotting
    for xtarg = 1:n_targ
        in_between = [fit_err{xtarg}(:,1)', fliplr(fit_err{xtarg}(:,2)')];
        h = fill([fitx fliplr(fitx)], in_between, xcond_color{xtarg});
        set(h,'facealpha', .2)
        
        for i = 1:2
            plot(fitx, fit_err{xtarg}(:,i), '-', 'Color', xcond_color{xtarg},'LineWidth', line_w_err); 
        end
    end
    
    set(gca,'xlim', x_lim, 'ylim', y_lim);
    set(gca,'XTick', x_tick)
    set(gca,'YTick', y_lim(1):0.1:y_lim(end))
    
    title(sprintf('Timecourse: %s', condition_names{xcond}));
    xlabel('Volume (tr)');
    ylabel('classifier evidence');
    
    %*************** onset lines
    plot([dur_stim dur_stim], y_lim,'--', 'Color', xcolor{2})
    plot([dur_stim+dur_manipulation, dur_stim+dur_manipulation], y_lim,'--', 'Color', xcolor{2})
    
    text(1, 0.05, 'stim onset');
    text(dur_stim, 0.05, 'operation onset');
    text(dur_stim + dur_manipulation, 0.05, 'fixation onset');
    
    %% ============= SAVE FIGURE
    fig_fname = fullfile(dirs.mvpa.parse{xph}, ...
        sprintf('plot_continous_operation_all_timecourse_%s_%s', ...
        xcond_name{xcond}, basename));
    
    savefig(LogRegs_fig, sprintf('%s.fig', fig_fname));
    set(gcf,'PaperUnits','inches','PaperPosition',fig_rect/100)
    print('-djpeg', sprintf('%s.jpg',fig_fname), '-r100')
    saveas(LogRegs_fig, sprintf('%s.jpg',fig_fname), 'jpg')
    
    close(LogRegs_fig);
    
end

end