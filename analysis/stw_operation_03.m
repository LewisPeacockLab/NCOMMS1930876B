function[] = clearmem_stw_operation_03(args, dirs)
% step 3: RSA timecourse: test on study: 1st + 2nd level

%% ============= UNPACK ARGS.
xph               = args.xphase;
subject_list      = args.subject_list;
xsub_groups       = args.filtered_subs;
n_subs            = length(xsub_groups);
xout_dir          = dirs.mvpa.group.out{xph};
%*************** output basename
basename          = args.analysis_basename;

%% ============= UNPACK PARAMETERS
%*************** regressors: shift or beta
xheader          = args.index{xph}.header;%from study
xmatrix          = args.index{xph}.matrix;
xparam           = args.index{xph}.param;

%*************** unpack parameters
condition_names  = {'maintain','replace category',...
    'replace subcategory','target suppress','global clear'};
n_condition      = length(condition_names);
n_runs           = xparam.n_runs;
n_trials         = xparam.n_trials;
dur_stim         = xparam.dur_stim;
dur_manipulation = xparam.dur_manipulation;% 6 TR
fix_win          = 23:31;
xchance          = (1/n_condition);

%*************** sliding time window param
it_wins          = 1:(args.tc_tr_disp+1);
stw_name         = {'grp_stw_all','grp_stw_filtered'};

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%============= 2ND LEVEL PARSED MVPAOUT STRUCTURE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fname = fullfile(xout_dir, sprintf('grp_mvpaout_%s.mat', basename));

%% ============= 1ST LEVEL
for xsub = xsub_groups
    %*************** setup subject & directories
    args.subject_id = subject_list(xsub).name;
    dirs            = setup_directory(dirs, args);
    
    fprintf('(+) 2nd level parsed mvpaout: s%s_%s\n', num2str(xsub), args.subject_id);
    
    %% ============= LOAD MVPAOUT
    %*************** stw timecourse
    mvpaout_name = sprintf('%s/mvpaout_%s.mat', dirs.mvpa.parse{xph}, basename);
    xparse       = load(mvpaout_name);%'tc_stw'
    
    grp_mvpaout{xsub} = xparse; %#ok<*AGROW,*NASGU>
    
end

%%*************** 2ND LEVEL SAVE
save(fname, 'grp_mvpaout', '-v7.3');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%============= 2ND LEVEL TIMECOURSE OF EVIDENCE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ============= CLASSIFIER ACCURACY IN FIXATION WINDOW
for ff = 1:2%1_all, 2_filtered (n~=n+1 operation)
    for xcond = 1:n_condition % desired
        for xtr = fix_win
            for it = 1:n_subs
                xsub = xsub_groups(it);
                xevi = [];
                for xtarg = 1:n_condition % guessed
                    if ff==1, t_evi = grp_mvpaout{xsub}.tc_stw.operation{xcond}.decoded_operation{xtarg}.tr{xtr}';
                    else,     t_evi = grp_mvpaout{xsub}.tc_stw.filter.operation{xcond}.decoded_operation{xtarg}.tr{xtr}'; end
                    
                    xevi = horzcat(xevi, t_evi);
                end
                
                [~, xwhich] = max(xevi,[],2);
                
                grp_stw{ff}.fixation_operation.accuracy(it, xcond) = ...
                    mean(ismember(xwhich, xcond));
            end
        end
    end
end

%% *************** one-sample ttest
xfile    = fullfile(xout_dir, sprintf('STW_stats_%s_n%s.txt', basename, num2str(n_subs)));
xout_txt = fopen(xfile, 'w+');

fprintf(xout_txt,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
fprintf(xout_txt,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
fprintf(xout_txt,'* OPERATION STW CLASSIFICATION: fixation window\n');
fprintf(xout_txt,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
fprintf(xout_txt,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n\n');
fprintf(xout_txt,'written from: clearmem_stw_operation_03.m\n\n');

fprintf(xout_txt,'======================================================\n');
fprintf(xout_txt,'* one-sample ttest\n');
fprintf(xout_txt,'======================================================\n\n');
fprintf(xout_txt,'chance level = %1.2f\n', xchance);

for ff = 1:2
    fprintf(xout_txt,'\n=====================================\n');
    fprintf(xout_txt,'* %s\n', stw_name{ff});
    fprintf(xout_txt,'=====================================\n\n');
    
    for xcond = 1:n_condition
        clear xacc xp xstats
        xacc = grp_stw{ff}.fixation_operation.accuracy(:, xcond);
        
        [~, xp, ~, xstats] = ttest(xacc, xchance);
        
        fprintf(xout_txt,'*************************************\n');
        fprintf(xout_txt,'*** %s: M = %1.2f, SE = %1.3f\nT(%s) = %1.2f, P=%1.3f\n\n', ...
            condition_names{xcond}, mean(xacc), std(xacc)/sqrt(n_subs), ...
            num2str(xstats.df), xstats.tstat, xp);
    end
    
    clear xacc xp xstats
    xacc = mean(grp_stw{ff}.fixation_operation.accuracy, 2);
    
    [~, xp, ~, xstats] = ttest(xacc, xchance);
    
    fprintf(xout_txt,'*************************************\n');
    fprintf(xout_txt,'*** average: M = %1.2f, SE = %1.3f\nT(%s) = %1.2f, P=%1.3f\n\n', ...
        mean(xacc), std(xacc)/sqrt(n_subs), ...
        num2str(xstats.df), xstats.tstat, xp);
    
end

%% ============= SEPARATED OPERATIONS
for ff = 1:2%1_all, 2_filtered (n~=n+1 operation)
    for xcond = 1:n_condition
        for xtarg = 1:n_condition
            for xtr = it_wins
                grp_stw{ff}.operation{xcond}.decoded_operation{xtarg}.tr{xtr} = [];
            end
        end
    end
end

for xcond = 1:n_condition
    for xtarg = 1:n_condition
        for xtr = it_wins
            for xsub = xsub_groups
                grp_stw{1}.operation{xcond}.decoded_operation{xtarg}.tr{xtr} = ...
                    horzcat(grp_stw{1}.operation{xcond}.decoded_operation{xtarg}.tr{xtr}, ...
                    mean(grp_mvpaout{xsub}.tc_stw.operation{xcond}.decoded_operation{xtarg}.tr{xtr}));
                
                %*************** filtered
                grp_stw{2}.operation{xcond}.decoded_operation{xtarg}.tr{xtr} = ...
                    horzcat(grp_stw{2}.operation{xcond}.decoded_operation{xtarg}.tr{xtr}, ...
                    mean(grp_mvpaout{xsub}.tc_stw.filter.operation{xcond}.decoded_operation{xtarg}.tr{xtr}));
            end
        end
    end
end

for xcond = 1:n_condition
    for xtarg = 1:n_condition
        for xtr = it_wins
            for ff = 1:2
                grp_stw{ff}.operation{xcond}.mean(xtarg, xtr) = ...
                    mean(grp_stw{ff}.operation{xcond}.decoded_operation{xtarg}.tr{xtr});
                grp_stw{ff}.operation{xcond}.se(xtarg, xtr) = ...
                    std(grp_stw{ff}.operation{xcond}.decoded_operation{xtarg}.tr{xtr})/sqrt(length(xsub_groups));
            end
        end
    end
end

%% *************** fixation bin
for ff = 1:2%2_filtered (removed trials when N=N+1 operations)
    for xcond = 1:n_condition
        for xtarg = 1:n_condition
            t_evi = [];
            for xtr = fix_win
                t_evi = horzcat(t_evi, grp_stw{ff}.operation{xcond}.decoded_operation{xtarg}.tr{xtr}');
            end
            
            grp_stw{ff}.target_fixation{xcond}(:, xtarg) = mean(t_evi, 2);
        end
    end
end

%% ============= TARGET-ONLY
for ff = 1:2
    for xcond = 1:n_condition
        for xtr = it_wins
            grp_stw{ff}.target_operation{xcond}.tr{xtr} = ...
                grp_stw{ff}.operation{xcond}.decoded_operation{xcond}.tr{xtr};
            
            grp_stw{ff}.targets_operation.mean(xcond, xtr) = ...
                mean(grp_stw{ff}.target_operation{xcond}.tr{xtr});
            grp_stw{ff}.targets_operation.se(xcond, xtr) = ...
                std(grp_stw{ff}.target_operation{xcond}.tr{xtr})/sqrt(length(xsub_groups));
        end
    end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%============= PLOTTING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%============= PLOT PARAMS
x_lim      = [it_wins(1) it_wins(end)];% + [-0.5 0.5];
y_lim      = [0 1];
x_tick     = it_wins;
y_ticks    = y_lim(1):0.1:y_lim(end);

n_ticks    = args.tc_tr_disp;
line_w     = 2;
rect_w     = 650;
rect_h     = 500;
font_size  = 10;

peak_tr    = 17:22;
%*************** separated operations
xcond_color  = args.cond_color;

%% ============= SEPARATED OPERATIONS
fig_rect = [0 0 (rect_w * 3) (rect_h * 2)];

for ff = 1:2%1_all trials, 2_filtered (remove if N=N+1 operation)
    
    xfig = figure;
    set(xfig, 'Position', fig_rect)
    
    for xcond = 1:(n_condition + 1)
        clear xmean xstd fity
        
        subplot(2, 3, xcond)
        
        %*************** mean line plots
        fitx = linspace(1, n_ticks, n_ticks*10);
        
        for xtarg = 1:n_condition
            if xcond ~= (n_condition + 1)
                xmean{xtarg} = grp_stw{ff}.operation{xcond}.mean(xtarg, x_tick);
                xse{xtarg}   = grp_stw{ff}.operation{xcond}.se(xtarg, x_tick); %#ok<*AGROW>
            else
                xmean{xtarg} = grp_stw{ff}.targets_operation.mean(xtarg, x_tick);
                xse{xtarg}   = grp_stw{ff}.targets_operation.se(xtarg, x_tick);
            end
            
            fity{xtarg} = interp1(x_tick, xmean{xtarg}, fitx,'spline');
            
            if xtarg==xcond
                plot(fitx, fity{xtarg}, '-', 'Color', xcond_color{xtarg}, 'LineWidth', line_w + 2); hold on;
            else
                plot(fitx, fity{xtarg}, '-', 'Color', xcond_color{xtarg}, 'LineWidth', line_w); hold on;
            end
        end
        
        %*************** legend
        if xcond == (n_condition+1)
            xlegend     = condition_names;
            lg          = legend(xlegend);
            lg.Location = 'SouthEast';
            lg.FontSize = font_size;
            
            legend(xlegend,'AutoUpdate','off')
        end
        grid on
        
        %*************** std error-bar filling
        for xtarg = 1:n_condition
            clear xerr h
            xerr(1,:) = xmean{xtarg} - xse{xtarg};
            xerr(2,:) = xmean{xtarg} + xse{xtarg};
            
            for i = 1:2, fit_err(:,i) = interp1(x_tick, xerr(i,:), fitx,'spline'); end
            
            in_between = [fit_err(:,1)', fliplr(fit_err(:,2)')];
            h = fill([fitx fliplr(fitx)], in_between, xcond_color{xtarg});
            h.EdgeColor = xcond_color{xtarg};
            h.EdgeAlpha = .2;
            h.FaceAlpha = .2;
        end
        
        %*************** peak filling
        p = fill([peak_tr fliplr(peak_tr)], ...
            [ones(1,length(peak_tr))*y_lim(1) ones(1,length(peak_tr))*y_lim(2)], 'g');
        p.EdgeColor = 'g';
        p.EdgeAlpha = .1;
        p.FaceAlpha = .1;
        
        %*************** set plot
        set(gca,'xlim', x_lim, 'ylim', y_lim);
        set(gca,'XTick', x_tick, 'YTick', y_ticks, 'XTickLabel', x_tick, 'FontSize', font_size)
        
        if xcond ~= (n_condition+1)
            title(sprintf('STW: %s (N=%s)', condition_names{xcond}, num2str(length(xsub_groups))));
        else
            title(sprintf('TARGETS (N=%s)', num2str(length(xsub_groups))));
        end
        xlabel('Volume (tr)');
        ylabel('classifier evidence');
        
        %*************** onset lines
        plot([dur_stim dur_stim] + 1, y_lim,'-', 'Color', 'r')
        plot([dur_stim+dur_manipulation, dur_stim+dur_manipulation] + 1, y_lim,'-', 'Color', 'r')
        
        text(1, y_lim(1)+0.05, 'stim', 'Color','r', 'FontSize', font_size);
        text(dur_stim + 1, y_lim(1)+0.05, 'operation','Color','r', 'FontSize', font_size);
        text(dur_stim + dur_manipulation + 1, y_lim(1)+0.05, 'fixation','Color','r', 'FontSize', font_size);
    end
    
    %% ============= SAVE FIGURE
    fig_fname = fullfile(dirs.mvpa.group.out{xph}, ...
        sprintf('%s_%s', stw_name{ff}, basename));
    
    savefig(xfig, sprintf('%s.fig', fig_fname));
    set(gcf,'PaperUnits','inches','PaperPosition',fig_rect/100)
    print('-dpng', sprintf('%s.png',fig_fname), '-r300')
    saveas(xfig, sprintf('%s.png',fig_fname), 'png')
    
    close(xfig);
    
end
end