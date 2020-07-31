function[] = grp_analysis_timecourse_workingmemory_category(args, grp_mvpaout, dirs)

%---------------------------------------------------------------------------
%*************** group: timecourse of the working memory (category|subcategory)
%---------------------------------------------------------------------------
% with 2017

xph = args.xphase;
fprintf('\n(+) category-based timecourse of working memory: %s\n', args.phase_name{xph});

xsub_groups = args.filtered_subs;%args.g_sub;

%% ============= UNPACK PARAMETERS
xindex          = args.index{xph};% param index from study

xparam          = xindex.param;
xcond_name      = xparam.conds_names;
n_condition     = length(xcond_name);
n_category      = xparam.n_category;
dur_stim        = xparam.dur_stim;
dur_manip       = xparam.dur_manipulation;
n_tc_trs        = xparam.n_tc_trs;
xalpha          = 0.05;

condition_names = {'maintain','replace category','replace subcategory','target suppress','global clear'};
category_names  = {'face','fruit','scene'};
class_index     = 1:n_category;

%*************** BASELINE CORRECTION
baseline_trs    = args.baseline_trs;

%*************** output basename
base_name       = sprintf('%s_n%s', args.analysis_basename, num2str(length(xsub_groups)));

%% ============= PLOT PARAMS

for xcond = 1:n_condition
        
    xcolor{xcond}{1}     = [238, 20, 91]/255;% targ
    
    if xcond~=2
        xcolor{xcond}{2} = [144, 144, 144]/255;% non_target_1
        xcolor{xcond}{3} = [72, 72, 72]/255;% non_target_2
    else
        xcolor{xcond}{2} = [0, 188, 182]/255;% newtarg
        xcolor{xcond}{3} = [144, 144, 144]/255;% non_target_1
    end
end

on_stim      = args.shift_TRs + 1;
on_operation = args.shift_TRs + 1 + dur_stim;
on_fixation  = on_operation + dur_manip;

x_tick       = 1:n_tc_trs;
x_ticklable  = 0:(n_tc_trs-1);
x_lim        = [x_tick(1) x_tick(end)];% + [-0.5 0.5];

%*************** legend
for i = [1 3:5]
    legend_names{i} = {'targ','ntarg1','ntarg 2'};
end

legend_names{2} = {'targ','new targ','ntarg'};

xcond_color  = args.cond_color;
xcate_color  = args.cate_color;

xbase_color  = args.base_color;
xonset_color = args.onset_color;

%% ============= PARSE DATA
fname = fullfile(dirs.mvpa.group.out{xph}, sprintf('grp_cate_based_wmem_timecourse_%s.mat', base_name));

if args.parse_timecourse_cate
    
    %% ============= CORRELATION CONCAT TIMECOURSE: PER CATEGORY timecourse
    %*************** category/trial based
    % decode.concat_timecourse.condition{xcond}.targ{xcate}(xheader,tr)
    
    for xcond = 1:n_condition
        for xcate = 1:n_category
            clear xcorr xpvalue
            for xsub = 1:length(xsub_groups)
                it_sub = xsub_groups(xsub);
                clear xevidence
                
                % xevidence{1}:targ, xevidence{2} & xevidence{3}:nontarg
                % xcond=2: xevidence{new}{1}:targ, xevidence{new}{2}:newtarg & xevidence{new}{3}:nontarg
                if xcond~=2
                    for i = 1:n_category
                        xevidence{i} = ...
                            grp_mvpaout{it_sub}.decode.concat_timecourse.condition{xcond}.targ{xcate}(i+1,:);
                    end
                else
                    for xnew_cate = class_index(~ismember(class_index, xcate))
                        for i = 1:n_category
                            xevidence{xnew_cate}{i} = ...
                                grp_mvpaout{it_sub}.decode.concat_timecourse.condition{xcond}.targ{xcate}.newtarg{xnew_cate}(i+1,:);
                        end
                    end
                end
                
                %*************** targ vs. (nontarg/newtarg)
                if xcond~=2
                    for j = 1:2 % pair-wise
                        [xR, xP] = corrcoef(xevidence{1}, xevidence{j + 1});
                        
                        xcorr{j}(xsub)   = xR(1,2);
                        xpvalue{j}(xsub) = xP(1,2);
                    end
                else
                    for xnew_cate = class_index(~ismember(class_index, xcate))
                        for j = 1:2
                            [xR, xP] = corrcoef(xevidence{xnew_cate}{1}, xevidence{xnew_cate}{j + 1});
                            
                            xcorr{xnew_cate}{j}(xsub)   = xR(1,2);
                            xpvalue{xnew_cate}{j}(xsub) = xP(1,2);
                        end
                    end
                end
            end
            
            if xcond~=2
                for j = 1:2
                    timecourse_corr.condition{xcond}.targ{xcate}.ntarg{j}.corr      = xcorr{j};
                    timecourse_corr.condition{xcond}.targ{xcate}.ntarg{j}.pvalue    = xpvalue{j};
                    timecourse_corr.condition{xcond}.targ{xcate}.ntarg{j}.mean_corr = mean(xcorr{j});
                    timecourse_corr.condition{xcond}.targ{xcate}.ntarg{j}.se_corr   = std(xcorr{j})/sqrt(length(xsub_groups));
                end
            else
                for j = 1:2
                    for xnew_cate = class_index(~ismember(class_index, xcate))
                        timecourse_corr.condition{xcond}.targ{xcate}.newtarg{xnew_cate}{j}.corr      = xcorr{xnew_cate}{j};
                        timecourse_corr.condition{xcond}.targ{xcate}.newtarg{xnew_cate}{j}.pvalue    = xpvalue{xnew_cate}{j};
                        timecourse_corr.condition{xcond}.targ{xcate}.newtarg{xnew_cate}{j}.mean_corr = mean(xcorr{xnew_cate}{j});
                        timecourse_corr.condition{xcond}.targ{xcate}.newtarg{xnew_cate}{j}.se_corr   = std(xcorr{xnew_cate}{j})/sqrt(length(xsub_groups));
                    end
                end
            end
        end
    end
    
    %% ============= BASELINE CORRELATION CONCAT TIMECOURSE: PER CATEGORY timecourse
    %*************** category/trial based
    % decode.concat_timecourse.condition{xcond}.targ{xcate}(xheader,tr)
    
    for xcond = 1:n_condition
        for xcate = 1:n_category
            clear xcorr xpvalue trial_array rand_trial t_evidence xevidence
            for xsub = 1:length(xsub_groups)
                it_sub = xsub_groups(xsub);
                
                %*************** randomization
                if xcond~=2
                    trial_array = ...
                        unique(grp_mvpaout{it_sub}.decode.concat_timecourse.condition{xcond}.targ{xcate}(1,:));
                    % randomize trials for all 3 evidences
                    for rr = 1:n_category, rand_trial{rr} = shuffle(trial_array); end
                    
                    for rr = 1:n_category
                        xevidence{rr} = [];
                        
                        for xx = rand_trial{rr}
                            xunit = grp_mvpaout{it_sub}.decode.concat_timecourse.condition{xcond}.targ{xcate}(1,:)==xx;
                            xevidence{rr} = horzcat(xevidence{rr}, ...
                                grp_mvpaout{it_sub}.decode.concat_timecourse.condition{xcond}.targ{xcate}(rr+1,xunit));
                        end
                    end
                else
                    for xnew_cate = class_index(~ismember(class_index, xcate))
                        trial_array{xnew_cate} = ...
                            unique(grp_mvpaout{it_sub}.decode.concat_timecourse.condition{xcond}.targ{xcate}.newtarg{xnew_cate}(1,:));
                        for rr = 1:n_category, rand_trial{xnew_cate}{rr} = shuffle(trial_array{xnew_cate}); end
                        
                        for rr = 1:n_category
                            xevidence{xnew_cate}{rr} = [];
                            
                            for xx = rand_trial{xnew_cate}{rr}
                                xunit = grp_mvpaout{it_sub}.decode.concat_timecourse.condition{xcond}.targ{xcate}.newtarg{xnew_cate}(1,:)==xx;
                                xevidence{xnew_cate}{rr} = horzcat(xevidence{xnew_cate}{rr}, ...
                                    grp_mvpaout{it_sub}.decode.concat_timecourse.condition{xcond}.targ{xcate}.newtarg{xnew_cate}(rr+1,xunit));
                            end
                        end
                    end
                end
                
                %*************** correlation
                % xevidence{1}:targ, xevidence{2} & xevidence{3}:nontarg
                % xcond=2: xevidence{new}{1}:targ, xevidence{new}{2}:newtarg & xevidence{new}{3}:nontarg
                %*************** targ vs. (nontarg/newtarg)
                if xcond~=2
                    for j = 1:2 % pair-wise
                        [xR, xP] = corrcoef(xevidence{1}, xevidence{j + 1});
                        
                        xcorr{j}(xsub)   = xR(1,2);
                        xpvalue{j}(xsub) = xP(1,2);
                    end
                else
                    for xnew_cate = class_index(~ismember(class_index, xcate))
                        for j = 1:2
                            [xR, xP] = corrcoef(xevidence{xnew_cate}{1}, xevidence{xnew_cate}{j + 1});
                            
                            xcorr{xnew_cate}{j}(xsub)   = xR(1,2);
                            xpvalue{xnew_cate}{j}(xsub) = xP(1,2);
                        end
                    end
                end
            end
            
            if xcond~=2
                for j = 1:2
                    timecourse_corr.random.condition{xcond}.targ{xcate}.ntarg{j}.corr      = xcorr{j};
                    timecourse_corr.random.condition{xcond}.targ{xcate}.ntarg{j}.pvalue    = xpvalue{j};
                    timecourse_corr.random.condition{xcond}.targ{xcate}.ntarg{j}.mean_corr = mean(xcorr{j});
                    timecourse_corr.random.condition{xcond}.targ{xcate}.ntarg{j}.se_corr   = std(xcorr{j})/sqrt(length(xsub_groups));
                end
            else
                for j = 1:2
                    for xnew_cate = class_index(~ismember(class_index, xcate))
                        timecourse_corr.random.condition{xcond}.targ{xcate}.newtarg{xnew_cate}{j}.corr      = xcorr{xnew_cate}{j};
                        timecourse_corr.random.condition{xcond}.targ{xcate}.newtarg{xnew_cate}{j}.pvalue    = xpvalue{xnew_cate}{j};
                        timecourse_corr.random.condition{xcond}.targ{xcate}.newtarg{xnew_cate}{j}.mean_corr = mean(xcorr{xnew_cate}{j});
                        timecourse_corr.random.condition{xcond}.targ{xcate}.newtarg{xnew_cate}{j}.se_corr   = std(xcorr{xnew_cate}{j})/sqrt(length(xsub_groups));
                    end
                end
            end
        end
    end
    
    %% ============= CORRELATION: PLOTTING
    y_lim = [-0.7 0];
    
    for xcond = 1:n_condition
        
        if xcond~=2, fig_rect = [0 0 600 400];
        else,        fig_rect = [0 0 1200 400]; end
        
        LogRegs_fig = figure;
        set(LogRegs_fig, 'Position', fig_rect)
        
        clear xcorr xse
        
        if xcond~=2
            for xcate = 1:n_category
                for j = 1:2
                    xcorr{xcond}(xcate,j) = timecourse_corr.condition{xcond}.targ{xcate}.ntarg{j}.mean_corr;
                    xse{xcond}(xcate,j)   = timecourse_corr.condition{xcond}.targ{xcate}.ntarg{j}.se_corr;
                    
                    xcorr{xcond}(xcate,j+2) = timecourse_corr.random.condition{xcond}.targ{xcate}.ntarg{j}.mean_corr;
                    xse{xcond}(xcate,j+2)   = timecourse_corr.random.condition{xcond}.targ{xcate}.ntarg{j}.se_corr;
                end
            end
        else
            for xcate = 1:n_category
                for j = 1:2
                    t_new = class_index(~ismember(class_index, xcate));
                    
                    for nn = 1:2
                        xnew_cate = t_new(nn);
                        
                        xunit = nn + 2 * (xcate-1);
                        xcorr{xcond}(xunit,j) = ...
                            timecourse_corr.condition{xcond}.targ{xcate}.newtarg{xnew_cate}{j}.mean_corr;
                        xse{xcond}(xunit,j) = ...
                            timecourse_corr.condition{xcond}.targ{xcate}.newtarg{xnew_cate}{j}.se_corr;
                        xcorr{xcond}(xunit,j+2) = ...
                            timecourse_corr.random.condition{xcond}.targ{xcate}.newtarg{xnew_cate}{j}.mean_corr;
                        xse{xcond}(xunit,j+2) = ...
                            timecourse_corr.random.condition{xcond}.targ{xcate}.newtarg{xnew_cate}{j}.se_corr;
                    end
                end
            end
        end
        
        ax = axes;
        h  = bar(xcorr{xcond},'BarWidth',1);
        if xcond~=2
            for i = 1:2, set(h(i), 'FaceColor', xcolor{xcond}{1}); end
            for i = 3:4, set(h(i), 'FaceColor', xcolor{xcond}{2}); end
        else
            for i = 1:2, set(h(i), 'FaceColor', xcolor{xcond}{i}); end
            for i = 3:4, set(h(i), 'FaceColor', xcolor{xcond}{3}); end
        end
        hold on;
        
        % Finding the number of groups and the number of bars in each group
        ngroups    = size(xcorr{xcond}, 1);
        nbars      = size(xcorr{xcond}, 2);
        % Calculating the width for each bar group
        groupwidth = min(0.8, nbars/(nbars + 1.5));
        % Set the position of each error bar in the centre of the main bar
        % Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
        for i = 1:nbars
            % Calculate center of each bar
            x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
            errorbar(x, xcorr{xcond}(:,i), xse{xcond}(:,i), 'k', 'linestyle', 'none');
        end
        
        ax.YGrid = 'on';
        ax.GridLineStyle = '-';
        ax.YLim = y_lim;
        
        if xcond~=2
            xticks(ax, 1:3);
            xticklabels(ax, category_names);
        else
            xticks(ax, 1:6);
            xticklabels(ax, sort([category_names category_names]));
        end
        set(ax, 'Ydir', 'reverse')
        
        title(sprintf('%s (N=%s)', xcond_name{xcond}, num2str(length(xsub_groups))));
        xlabel ('Target category');
        ylabel ('correlation coefficient');
        
        if xcond~=2
            xlegend = {'ntarg1','ntarg2','b(ntarg1)','b(ntarg2)'};
            legend('ntarg1','ntarg2','b(ntarg1)','b(ntarg2)');
        else
            xlegend = {'new','ntarg','b(new)','b(ntarg)'};
            legend('new','ntarg','b(new)','b(ntarg)');
        end

        lg.Location    = 'BestOutside';
%         lg             = legend(xlegend);
%         legend(xlegend,'AutoUpdate','off')
        
        hold off;
        
        %% ============= SAVE FIGURE
        fig_fname = fullfile(dirs.mvpa.group.out{xph}, ...
            sprintf('grp_correl_category_timecourse_%s_%s', xcond_name{xcond}, base_name));
        
        savefig(LogRegs_fig, sprintf('%s.fig', fig_fname));
        set(gcf,'PaperUnits','inches','PaperPosition',fig_rect/100)
        print('-djpeg', sprintf('%s.jpg',fig_fname), '-r100')
        saveas(LogRegs_fig, sprintf('%s.jpg',fig_fname), 'jpg')
        
        close(LogRegs_fig);
        
    end

    %% ============= TIMECOURSE: RANDOM EFFECT PER CATEGORY
    %*************** category/trial based
    % decode.cate_timecourse.condition{xcond}.category{xcate}.target{xx}.evidence{xtr}
    % xx: 1_target, 2_nontarg1/newtarg, 3_nontarg2/nontarg
    % decode.cate_timecourse.condition{xcond}.category{xcate}.newcategory{xnew_targ}.target{xx}.evidence{xtr}
    
    for xcond = 1:n_condition
        for xcate = 1:n_category
            for xtr = 1:n_tc_trs
                for xx = 1:3
                    clear xevidence
                    
                    if xcond~=2, xevidence = [];
                    else, for i=1:3, xevidence{i} = []; end
                    end
                    
                    %*************** random effect
                    for xsub = xsub_groups
                        if xcond~=2
                            xevidence = horzcat(xevidence, ...
                                mean(grp_mvpaout{xsub}.decode.cate_timecourse.condition{xcond}.category{xcate}.target{xx}.evidence{xtr}));
                        else
                            for xnew_cate = class_index(~ismember(class_index, xcate))
                                xevidence{xnew_cate} = horzcat(xevidence{xnew_cate}, ...
                                    mean(grp_mvpaout{xsub}.decode.cate_timecourse.condition{xcond}.category{xcate}.newcategory{xnew_cate}.target{xx}.evidence{xtr}));
                            end
                        end
                    end
                    
                    if xcond~=2
                        timecourse_cate.evidence{xcond}.cate{xcate}.targ{xx}{xtr} = xevidence;
                        timecourse_cate.mean{xcond}.cate{xcate}.targ(xx, xtr)     = mean(xevidence);
                        timecourse_cate.se{xcond}.cate{xcate}.targ(xx, xtr)       = std(xevidence)/sqrt(length(xsub_groups));
                    else
                        for xnew_cate = class_index(~ismember(class_index, xcate))
                            timecourse_cate.evidence{xcond}.cate{xcate}.newtarg{xnew_cate}{xx}{xtr} = xevidence{xnew_cate};
                            timecourse_cate.mean{xcond}.cate{xcate}.newtarg{xnew_cate}(xx, xtr)     = mean(xevidence{xnew_cate});
                            timecourse_cate.se{xcond}.cate{xcate}.newtarg{xnew_cate}(xx, xtr)       = std(xevidence{xnew_cate})/sqrt(length(xsub_groups));
                        end
                    end
                end
            end
        end
    end
    
    %% ============= STATS: RANDOM
    for xcond = 1:n_condition
        for xcate = 1:n_category
            for xtr = 1:n_tc_trs
                n_targ = 3;
                
                if xcond~=2
                    timecourse_cate.t_test{xcond}.cate{xcate}.tr{xtr}.pvalue = nan(n_targ, n_targ);
                else
                    for xnew_cate = class_index(~ismember(class_index, xcate))
                        timecourse_cate.t_test{xcond}.cate{xcate}.newtarg{xnew_cate}.tr{xtr}.pvalue = nan(n_targ, n_targ);
                    end
                end
                
                for xcol = 1:(n_targ-1)
                    for xrow = (xcol+1):n_targ
                        if xcond~=2
                            xtarg_col = timecourse_cate.evidence{xcond}.cate{xcate}.targ{xcol}{xtr};%array_1
                            xtarg_row = timecourse_cate.evidence{xcond}.cate{xcate}.targ{xrow}{xtr};%array_2
                            
                            %*************** ttest
                            [xhypothesis, xpvalue, ~, xstats] = ttest(xtarg_col, xtarg_row, 'Alpha', 0.05);
                            
                            timecourse_cate.t_test{xcond}.cate{xcate}.tr{xtr}.pvalue(xrow, xcol)     = xpvalue;
                            timecourse_cate.t_test{xcond}.cate{xcate}.tr{xtr}.hypothesis(xrow, xcol) = xhypothesis;
                            timecourse_cate.t_test{xcond}.cate{xcate}.tr{xtr}.stats{xrow, xcol}      = xstats;
                        else
                            for xnew_cate = class_index(~ismember(class_index, xcate))
                                xtarg_col = timecourse_cate.evidence{xcond}.cate{xcate}.newtarg{xnew_cate}{xcol}{xtr};%array_1
                                xtarg_row = timecourse_cate.evidence{xcond}.cate{xcate}.newtarg{xnew_cate}{xrow}{xtr};%array_2
                                
                                %*************** ttest
                                [xhypothesis, xpvalue, ~, xstats] = ttest(xtarg_col, xtarg_row, 'Alpha', 0.05);
                                
                                timecourse_cate.t_test{xcond}.cate{xcate}.newtarg{xnew_cate}.tr{xtr}.pvalue(xrow, xcol)     = xpvalue;
                                
                                timecourse_cate.t_test{xcond}.cate{xcate}.newtarg{xnew_cate}.tr{xtr}.hypothesis(xrow, xcol) = xhypothesis;
                                timecourse_cate.t_test{xcond}.cate{xcate}.newtarg{xnew_cate}.tr{xtr}.stats{xrow, xcol}      = xstats;
                            end
                        end
                    end
                end
            end
        end
    end
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%============= BASELINE CORRECTION: RANDOM
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % timecourse_cate.evidence{xcond}.cate{xcate}.targ{xx}{xtr}
    % timecourse_cate.baseline.evidence{xcond}.cate{xcate}.targ{xx}{xtr}
    % timecourse_cate.evidence{xcond}.cate{xcate}.newtarg{xnew_cate}{xx}{xtr}
    % timecourse_cate.baseline.evidence{xcond}.cate{xcate}.newtarg{xnew_cate}{xx}{xtr}
    n_targ = 3;
    
    for xcond = 1:n_condition
        for xcate = 1:n_category
            clear xbase_evidence
            if xcond~=2
                %*************** baseline_trs: 6TRs
                xbase_evidence = zeros(n_targ, baseline_trs);
                for xx = 1:n_targ
                    for xtr = 1:baseline_trs
                        xbase_evidence(xx, xtr) = mean(timecourse_cate.evidence{xcond}.cate{xcate}.targ{xx}{xtr});
                    end
                    
                    mean_baseline(xx) = mean(xbase_evidence(xx,:));
                end
                
                %*************** corrected evidence
                for xtr = 1:n_tc_trs
                    for xx = 1:n_targ
                        xcorrect_evidence = ...
                            timecourse_cate.evidence{xcond}.cate{xcate}.targ{xx}{xtr} - mean_baseline(xx);
                        
                        timecourse_cate.baseline.evidence{xcond}.cate{xcate}.targ{xx}{xtr} = xcorrect_evidence;
                        timecourse_cate.baseline.mean{xcond}.cate{xcate}.targ(xx, xtr)     = mean(xcorrect_evidence);
                        timecourse_cate.baseline.se{xcond}.cate{xcate}.targ(xx, xtr)       = std(xcorrect_evidence)/sqrt(length(xsub_groups));
                    end
                end
            else
                for xnew_cate = class_index(~ismember(class_index, xcate))
                    %*************** baseline_trs: 6TRs
                    xbase_evidence = zeros(n_targ, baseline_trs);
                    for xx = 1:n_targ
                        for xtr = 1:baseline_trs
                            xbase_evidence(xx, xtr) = mean(timecourse_cate.evidence{xcond}.cate{xcate}.newtarg{xnew_cate}{xx}{xtr});
                        end
                        
                        mean_baseline(xx) = mean(xbase_evidence(xx,:));
                    end
                    
                    %*************** corrected evidence
                    for xtr = 1:n_tc_trs
                        for xx = 1:n_targ
                            xcorrect_evidence = ...
                                timecourse_cate.evidence{xcond}.cate{xcate}.newtarg{xnew_cate}{xx}{xtr} - mean_baseline(xx);
                            
                            timecourse_cate.baseline.evidence{xcond}.cate{xcate}.newtarg{xnew_cate}{xx}{xtr} = xcorrect_evidence;
                            timecourse_cate.baseline.mean{xcond}.cate{xcate}.newtarg{xnew_cate}(xx, xtr)     = mean(xcorrect_evidence);
                            timecourse_cate.baseline.se{xcond}.cate{xcate}.newtarg{xnew_cate}(xx, xtr)       = std(xcorrect_evidence)/sqrt(length(xsub_groups));
                        end
                    end
                end
            end
        end
    end
    
    %% ============= STATS: RANDOM BASELINE CORRECTED
    for xcond = 1:n_condition
        for xcate = 1:n_category
            for xtr = 1:n_tc_trs
                n_targ = 3;
                
                if xcond~=2
                    timecourse_cate.baseline.t_test{xcond}.cate{xcate}.tr{xtr}.pvalue = nan(n_targ, n_targ);
                else
                    for xnew_cate = class_index(~ismember(class_index, xcate))
                        timecourse_cate.baseline.t_test{xcond}.cate{xcate}.newtarg{xnew_cate}.tr{xtr}.pvalue = nan(n_targ, n_targ);
                    end
                end
                
                for xcol = 1:(n_targ-1)
                    for xrow = (xcol+1):n_targ
                        if xcond~=2
                            xtarg_col = timecourse_cate.baseline.evidence{xcond}.cate{xcate}.targ{xcol}{xtr};%array_1
                            xtarg_row = timecourse_cate.baseline.evidence{xcond}.cate{xcate}.targ{xrow}{xtr};%array_2
                            
                            %*************** ttest
                            [xhypothesis, xpvalue, ~, xstats] = ttest(xtarg_col, xtarg_row, 'Alpha', 0.05);
                            
                            timecourse_cate.baseline.t_test{xcond}.cate{xcate}.tr{xtr}.pvalue(xrow, xcol)     = xpvalue;
                            timecourse_cate.baseline.t_test{xcond}.cate{xcate}.tr{xtr}.hypothesis(xrow, xcol) = xhypothesis;
                            timecourse_cate.baseline.t_test{xcond}.cate{xcate}.tr{xtr}.stats{xrow, xcol}      = xstats;
                        else
                            for xnew_cate = class_index(~ismember(class_index, xcate))
                                xtarg_col = timecourse_cate.baseline.evidence{xcond}.cate{xcate}.newtarg{xnew_cate}{xcol}{xtr};%array_1
                                xtarg_row = timecourse_cate.baseline.evidence{xcond}.cate{xcate}.newtarg{xnew_cate}{xrow}{xtr};%array_2
                                
                                %*************** ttest
                                [xhypothesis, xpvalue, ~, xstats] = ttest(xtarg_col, xtarg_row, 'Alpha', 0.05);
                                
                                timecourse_cate.baseline.t_test{xcond}.cate{xcate}.newtarg{xnew_cate}.tr{xtr}.pvalue(xrow, xcol)     = xpvalue;
                                timecourse_cate.baseline.t_test{xcond}.cate{xcate}.newtarg{xnew_cate}.tr{xtr}.hypothesis(xrow, xcol) = xhypothesis;
                                timecourse_cate.baseline.t_test{xcond}.cate{xcate}.newtarg{xnew_cate}.tr{xtr}.stats{xrow, xcol}      = xstats;
                            end
                        end
                    end
                end
            end
        end
    end
    
    %% ============= TIMECOURSE: TARGET ONLY IN CATEGORY FOR EACH CONDITION
    for xcond = 2:n_condition
        for xcate = 1:n_category
            
            %*************** target evidence
            if xcond~=2
                for xtr = 1:n_tc_trs
                    clear xevidence
                    xbs = timecourse_cate.baseline.evidence{1}.cate{xcate}.targ{1}{xtr};
                    xevidence = timecourse_cate.baseline.evidence{xcond}.cate{xcate}.targ{1}{xtr} - xbs;
                    
                    timecourse_cate.target.condition.evidence{xcond}.cate{xcate}.targ{xtr} = xevidence;
                    timecourse_cate.target.condition.mean{xcond}.cate{xcate}.targ(xtr)     = mean(xevidence);
                    timecourse_cate.target.condition.se{xcond}.cate{xcate}.targ(xtr)       = std(xevidence)/sqrt(length(xsub_groups));
                end
            else
                for xtr = 1:n_tc_trs
                    clear xevidence t_evi
                    for i = 1:2
                        t_new      = class_index(~ismember(class_index, xcate));
                        xnew_cate  = t_new(i);
                        t_evi(i,:) = timecourse_cate.baseline.evidence{xcond}.cate{xcate}.newtarg{xnew_cate}{1}{xtr};
                    end
                    
                    xevidence = mean(t_evi);
                    
                    timecourse_cate.target.condition.evidence{xcond}.cate{xcate}.targ{xtr} = xevidence;
                    timecourse_cate.target.condition.mean{xcond}.cate{xcate}.targ(xtr)     = mean(xevidence);
                    timecourse_cate.target.condition.se{xcond}.cate{xcate}.targ(xtr)       = std(xevidence)/sqrt(length(xsub_groups));
                end
            end
        end
    end
    
    %% ============= STATS: TARGET ONLY IN CATEGORY FOR EACH CONDITION
    n_targ = n_category;
    
    for xcond = 2:n_condition
        for xcate = 1:n_category
            for xtr = 1:n_tc_trs
                
                timecourse_cate.target.condition.t_test{xcond}.targ{xtr}.pvalue = nan(n_targ, n_targ);
                
                for xcol = 1:(n_targ-1)
                    for xrow = (xcol+1):n_targ
                        xtarg_col = timecourse_cate.target.condition.evidence{xcond}.cate{xcol}.targ{xtr};
                        xtarg_row = timecourse_cate.target.condition.evidence{xcond}.cate{xrow}.targ{xtr};%array_2
                        
                        %*************** ttest
                        [xhypothesis, xpvalue, ~, xstats] = ttest(xtarg_col, xtarg_row, 'Alpha', 0.05);
                        
                        timecourse_cate.target.condition.t_test{xcond}.targ{xtr}.pvalue(xrow, xcol)     = xpvalue;
                        timecourse_cate.target.condition.t_test{xcond}.targ{xtr}.hypothesis(xrow, xcol) = xhypothesis;
                        timecourse_cate.target.condition.t_test{xcond}.targ{xtr}.stats{xrow, xcol}      = xstats;
                    end
                end
            end
        end
    end
    
    %% ============= TIMECOURSE: TARGET ONLY IN CONDITION FOR EACH CATEGORY
    for xcond = 2:n_condition
        for xcate = 1:n_category
            for xtr = 1:n_tc_trs
                clear xevidence
                
                xbs = timecourse_cate.target.condition.evidence{1}.cate{xcate}.targ{xtr};
                xevidence = timecourse_cate.target.condition.evidence{xcond}.cate{xcate}.targ{xtr} - xbs;
                
                timecourse_cate.target.cate.evidence{xcate}.condition{xcond}.targ{xtr} = xevidence;
                timecourse_cate.target.cate.mean{xcate}.condition{xcond}.targ(xtr)     = mean(xevidence);
                timecourse_cate.target.cate.se{xcate}.condition{xcond}.targ(xtr)       = std(xevidence)/sqrt(length(xsub_groups));
            end
        end
    end
    
    %% ============= STATS: TARGET ONLY IN CONDITION FOR EACH CATEGORY
    n_targ = n_condition;
    
    for xcate = 1:n_category
        for xcond = 2:n_condition
            for xtr = 1:n_tc_trs
                
                timecourse_cate.target.cate.t_test{xcate}.targ{xtr}.pvalue = nan(n_targ, n_targ);
                
                for xcol = 1:(n_targ-1)
                    for xrow = (xcol+1):n_targ
                        xtarg_col = timecourse_cate.target.cate.evidence{xcate}.condition{xcol}.targ{xtr};
                        xtarg_row = timecourse_cate.target.cate.evidence{xcate}.condition{xrow}.targ{xtr};%array_2
                        
                        %*************** ttest
                        [xhypothesis, xpvalue, ~, xstats] = ttest(xtarg_col, xtarg_row, 'Alpha', 0.05);
                        
                        timecourse_cate.target.cate.t_test{xcate}.targ{xtr}.pvalue(xrow, xcol)     = xpvalue;
                        timecourse_cate.target.cate.t_test{xcate}.targ{xtr}.hypothesis(xrow, xcol) = xhypothesis;
                        timecourse_cate.target.cate.t_test{xcate}.targ{xtr}.stats{xrow, xcol}      = xstats;
                    end
                end
            end
        end
    end
    
    %% ============= SAVE
    grp_timecourse.correlation = timecourse_corr;
    grp_timecourse.category    = timecourse_cate;  %#ok<*STRNU>
    
    save(fname, 'grp_timecourse');
    
else
    
    load(fname);%grp_timecourse
    
    timecourse_corr = grp_timecourse.correlation; %#ok<*NASGU>
    timecourse_cate = grp_timecourse.category; 
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%============= PLOTTING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ============= CONTINUOUS PLOT: BASELINE CORRECTED RANDOM EFFECT
y_lim    = [-0.5 0.6];
fig_rect = [0 0 3000 600];  
n_targ   = 3;

for xcond = [1 3:n_condition]%1:n_condition
    LogRegs_fig = figure;
    set(LogRegs_fig, 'Position', fig_rect)
        
    for xcate = 1:n_category
        clear xmean xstd fity xlegend
        
        subplot(1, 3, xcate)
        
        %*************** mean line plots
        fitx = linspace(1, n_tc_trs, n_tc_trs*10);
        
        for xtarg = 1:n_targ
            xmean{xtarg} = timecourse_cate.baseline.mean{xcond}.cate{xcate}.targ(xtarg, x_tick); %#ok<*AGROW>
            xstd{xtarg}  = timecourse_cate.baseline.se{xcond}.cate{xcate}.targ(xtarg, x_tick);
            
            fity{xtarg}  = interp1(x_tick, xmean{xtarg}, fitx,'spline');
            
            plot(fitx, fity{xtarg}, '-','Color', xcolor{xcond}{xtarg}, 'LineWidth', 2); hold on;
        end
        
        %*************** legend
        xnontarg = class_index(~ismember(class_index, xcate));
%         xlegend  = {category_names{xcate}, category_names{xnontarg(1)}, category_names{xnontarg(2)}};
        
        legend(category_names{xcate}, category_names{xnontarg(1)}, category_names{xnontarg(2)}, 'Location', 'SouthEast')
        legend(category_names{xcate}, category_names{xnontarg(1)}, category_names{xnontarg(2)}, 'AutoUpdate', 'off') 
        
        %*************** stats stars
        n = 0;
        for xcol = 1:(n_targ-1)
            for xrow = (xcol+1):n_targ
                xheight = (y_lim(2)-0.05) - (n * 0.02);
                text(-4, xheight, sprintf('%s vs. %s', legend_names{xcond}{xcol}, ...
                    legend_names{xcond}{xrow}));
                plot(x_lim, [xheight xheight], ':', 'Color', [0.75 0.75 0.75])
                
                for xtr = 1:n_tc_trs
                    xpvalue = timecourse_cate.baseline.t_test{xcond}.cate{xcate}.tr{xtr}.pvalue(xrow, xcol);
                    
                    if xpvalue <= xalpha
                        text(xtr-0.1, xheight, '*', 'FontSize', 20);
                    end
                    
                    plot([xtr xtr], y_lim, ':', 'Color', [0.75 0.75 0.75])
                    
                end
                
                n = n + 1;
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
        
        %*************** baseline 0
        plot(x_lim, [0 0], '--', 'Color', xbase_color)
        
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
        
        title(sprintf('%s: %s target (N=%s)', ...
            condition_names{xcond}, category_names{xcate}, num2str(length(xsub_groups))));
        xlabel('Volume (tr)');
        ylabel('classifier evidence');
        
        h1 = text(1.5, y_lim(1) + 0.05, '(1) stim onset', 'Color', xonset_color);
        h2 = text(dur_stim + 1.5, y_lim(1) + 0.05, '(2) operation onset', 'Color', xonset_color);
        h3 = text(dur_stim + dur_manip + 1.5, y_lim(1) + 0.05, '(3) fixation onset', 'Color', xonset_color);
        set(h1,'Rotation', 90); set(h2,'Rotation', 90); set(h3,'Rotation', 90);
        
        text(on_stim - 1, y_lim(1) + 0.05, sprintf('(1) shift %str', num2str(args.shift_TRs)));
        text(on_operation - 1, y_lim(1) + 0.05, sprintf('(2) shift %str', num2str(args.shift_TRs)));
        text(on_fixation - 1, y_lim(1) + 0.05, sprintf('(3) shift %str', num2str(args.shift_TRs)));
    end
    
    %% *************** save fig
    fig_fname = fullfile(dirs.mvpa.group.out{xph}, ...
        sprintf('grp_plot_cate_based_continous_wn_timecourse_random_%s_%s', ...
        xcond_name{xcond}, base_name));
    
    savefig(LogRegs_fig, sprintf('%s.fig', fig_fname));
    set(gcf,'PaperUnits','inches','PaperPosition', fig_rect/100)
    print('-djpeg', sprintf('%s.jpg',fig_fname), '-r100')
    saveas(LogRegs_fig, sprintf('%s.jpg',fig_fname), 'jpg')
    
    close(LogRegs_fig);
    
end

%% ============= CONTINUOUS PLOT: COND2: BASELINE CORRECTED RANDOM EFFECT
xcond    = 2;
fig_rect = [0 0 3000 1200];  

LogRegs_fig = figure;
set(LogRegs_fig, 'Position', fig_rect)

for xcate = 1:n_category
    clear xmean xstd fity xlegend
    
    %*************** mean line plots
    fitx = linspace(1, n_tc_trs, n_tc_trs*10);
    
    for it_new = 1:2
        t_new     = class_index(~ismember(class_index, xcate));
        xnew_cate = t_new(it_new);
        xnontarg  = class_index(~ismember(class_index, [xcate xnew_cate]));
        
        subplot(2, 3, xcate + 3 * (it_new-1))
        hold on;
        
        for xtarg = 1:n_targ
            xmean{xtarg} = timecourse_cate.baseline.mean{xcond}.cate{xcate}.newtarg{xnew_cate}(xtarg, x_tick);
            xstd{xtarg}  = timecourse_cate.baseline.se{xcond}.cate{xcate}.newtarg{xnew_cate}(xtarg, x_tick);
            
            fity{xtarg}  = interp1(x_tick, xmean{xtarg}, fitx,'spline');
            
            plot(fitx, fity{xtarg}, '-','Color', xcolor{xcond}{xtarg}, 'LineWidth', 2);
        end
        
        %*************** legend        
        legend(category_names{xcate}, category_names{xnew_cate}, category_names{xnontarg}, 'Location', 'SouthEast')
        legend(category_names{xcate}, category_names{xnew_cate}, category_names{xnontarg}, 'AutoUpdate','off')
        
        %*************** stats stars
        n = 0;
        for xcol = 1:(n_targ-1)
            for xrow = (xcol+1):n_targ
                xheight = (y_lim(2)-0.05) - (n * 0.02);
                text(-4, xheight, sprintf('%s vs. %s', legend_names{xcond}{xcol}, ...
                    legend_names{xcond}{xrow}));
                plot(x_lim, [xheight xheight], ':', 'Color', [0.75 0.75 0.75])
                
                for xtr = 1:n_tc_trs
                    xpvalue = ...
                        timecourse_cate.baseline.t_test{xcond}.cate{xcate}.newtarg{xnew_cate}.tr{xtr}.pvalue(xrow, xcol);
                    
                    if xpvalue <= xalpha
                        text(xtr-0.1, xheight, '*', 'FontSize', 20);
                    end
                    
                    plot([xtr xtr], y_lim, ':', 'Color', [0.75 0.75 0.75])
                    
                end
                
                n = n + 1;
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
        
        %*************** baseline 0
        plot(x_lim, [0 0], '--', 'Color', xbase_color)
    
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
        
        title(sprintf('%s: %s target (N=%s)', ...
            condition_names{xcond}, category_names{xcate}, num2str(length(xsub_groups))));
        xlabel('Volume (tr)');
        ylabel('classifier evidence');
        
        h1 = text(1.5, y_lim(1) + 0.05, '(1) stim onset', 'Color', xonset_color);
        h2 = text(dur_stim + 1.5, y_lim(1) + 0.05, '(2) operation onset', 'Color', xonset_color);
        h3 = text(dur_stim + dur_manip + 1.5, y_lim(1) + 0.05, '(3) fixation onset', 'Color', xonset_color);
        set(h1,'Rotation', 90); set(h2,'Rotation', 90); set(h3,'Rotation', 90);
        
        text(on_stim - 1, y_lim(1) + 0.05, sprintf('(1) shift %str', num2str(args.shift_TRs)));
        text(on_operation - 1, y_lim(1) + 0.05, sprintf('(2) shift %str', num2str(args.shift_TRs)));
        text(on_fixation - 1, y_lim(1) + 0.05, sprintf('(3) shift %str', num2str(args.shift_TRs)));
    end
end

%%*************** save fig
fig_fname = fullfile(dirs.mvpa.group.out{xph}, ...
    sprintf('grp_plot_cate_based_continous_wn_timecourse_random_%s_%s', ...
    xcond_name{xcond}, base_name));

savefig(LogRegs_fig, sprintf('%s.fig', fig_fname));
set(gcf,'PaperUnits','inches','PaperPosition', fig_rect/100)
print('-djpeg', sprintf('%s.jpg',fig_fname), '-r100')
saveas(LogRegs_fig, sprintf('%s.jpg',fig_fname), 'jpg')

close(LogRegs_fig);

%% ============= PLOT: TARGET ONLY IN CATEGORY FOR EACH CONDITION
y_lim    = [-0.5 0.6];
fig_rect = [0 0 1200 800]; 
n_targ   = n_category;

for xcond = 1:n_condition
    clear xmean xstd fity
    
    LogRegs_fig = figure;
    set(LogRegs_fig, 'Position', fig_rect)

    for xcate = 1:n_category                
        %*************** mean line plots
        fitx = linspace(1, n_tc_trs, n_tc_trs*10);
        
        xmean{xcate} = timecourse_cate.target.condition.mean{xcond}.cate{xcate}.targ(x_tick);
        xstd{xcate}  = timecourse_cate.target.condition.se{xcond}.cate{xcate}.targ(x_tick);
        
        fity{xcate}  = interp1(x_tick, xmean{xcate}, fitx,'spline');
        
        plot(fitx, fity{xcate}, '-','Color', xcate_color{xcate}, 'LineWidth', 2); 
        hold on;
    end
    
    %*************** legend
    legend('face','fruit','scene', 'Location', 'SouthEast')
    legend('face','fruit','scene', 'AutoUpdate', 'off')
    
    %*************** stats stars
    n = 0;
    for xcol = 1:(n_targ-1)
        for xrow = (xcol+1):n_targ
            xheight = (y_lim(2)-0.05) - (n * 0.02);
            text(-4, xheight, sprintf('%s vs. %s', category_names{xcol}, category_names{xrow}));
            plot(x_lim, [xheight xheight], ':', 'Color', [0.75 0.75 0.75])
            
            for xtr = 1:n_tc_trs
                xpvalue = timecourse_cate.target.condition.t_test{xcond}.targ{xtr}.pvalue(xrow, xcol);
                
                if xpvalue <= xalpha
                    text(xtr-0.1, xheight, '*', 'FontSize', 20);
                end
                
                plot([xtr xtr], y_lim, ':', 'Color', [0.75 0.75 0.75])
            end
            
            n = n + 1;
        end
    end
    
    %*************** std error-bar filling
    for xtarg = 1:n_targ
        clear xerr
        xerr(1,:) = xmean{xtarg} - xstd{xtarg};
        xerr(2,:) = xmean{xtarg} + xstd{xtarg};
        
        for i = 1:2, fit_err(:,i) = interp1(x_tick, xerr(i,:), fitx,'spline'); end
        
        in_between = [fit_err(:,1)', fliplr(fit_err(:,2)')];
        h = fill([fitx fliplr(fitx)], in_between, xcate_color{xtarg});
        set(h,'facealpha', .4)
        
        plot(fitx, fit_err(:,1), '-', 'Color', xcate_color{xtarg})
        plot(fitx, fit_err(:,2), '-', 'Color', xcate_color{xtarg})
    end
    
    %*************** baseline 0
    plot(x_lim, [0 0], '--', 'Color', xbase_color)
    
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
    
    title(sprintf('%s: target evidence (N=%s)', ...
        condition_names{xcond}, num2str(length(xsub_groups))));
    xlabel('Volume (tr)');
    ylabel('classifier evidence');
    
    h1 = text(1.5, y_lim(1) + 0.05, '(1) stim onset', 'Color', xonset_color);
    h2 = text(dur_stim + 1.5, y_lim(1) + 0.05, '(2) operation onset', 'Color', xonset_color);
    h3 = text(dur_stim + dur_manip + 1.5, y_lim(1) + 0.05, '(3) fixation onset', 'Color', xonset_color);
    set(h1,'Rotation', 90); set(h2,'Rotation', 90); set(h3,'Rotation', 90);
    
    text(on_stim - 1, y_lim(1) + 0.05, sprintf('(1) shift %str', num2str(args.shift_TRs)));
    text(on_operation - 1, y_lim(1) + 0.05, sprintf('(2) shift %str', num2str(args.shift_TRs)));
    text(on_fixation - 1, y_lim(1) + 0.05, sprintf('(3) shift %str', num2str(args.shift_TRs)));
    
    %% *************** save fig
    fig_fname = fullfile(dirs.mvpa.group.out{xph}, ...
        sprintf('grp_plot_cate_based_target_wn_timecourse_%s_%s', ...
        xcond_name{xcond}, base_name));
    
    savefig(LogRegs_fig, sprintf('%s.fig', fig_fname));
    set(gcf,'PaperUnits','inches','PaperPosition', fig_rect/100)
    print('-djpeg', sprintf('%s.jpg',fig_fname), '-r100')
    saveas(LogRegs_fig, sprintf('%s.jpg',fig_fname), 'jpg')
    
    close(LogRegs_fig);
    
end

%% ============= PLOT: TARGET ONLY IN CONDITION FOR EACH CATEGORY
y_lim    = [-0.6 0.3];
fig_rect = [0 0 2000 400];
n_targ   = n_condition;

LogRegs_fig = figure;
set(LogRegs_fig, 'Position', fig_rect)

for xcate = 1:3
    
    subplot(1,3,xcate);
    
    for xsel_conds = [2 4 5]
        
        clear xmean xstd fity
        
        %*************** mean line plots
        fitx = linspace(1, n_tc_trs, n_tc_trs*10);
        
        for xcond = xsel_conds
            xmean{xcond} = timecourse_cate.target.cate.mean{xcate}.condition{xcond}.targ(x_tick);
            xstd{xcond}  = timecourse_cate.target.cate.se{xcate}.condition{xcond}.targ(x_tick);
            
            fity{xcond}  = interp1(x_tick, xmean{xcond}, fitx,'spline');
            
            plot(fitx, fity{xcond}, '-','Color', xcond_color{xcond}, 'LineWidth', 2);
            hold on;
        end
        
        %*************** legend
        %             clear xlegend
        %             for i=1:length(xsel_conds); xlegend{i} = condition_names{xsel_conds(i)}; end
%         lg             = legend('replace','suppress','clear');
%         lg.Location    = 'BestOutside';
%         legend('replace','suppress','clear', 'AutoUpdate', 'off')
        
        %*************** stats stars
        n = 0;
        for i = 1:(length(xsel_conds)-1)
            xcol = xsel_conds(i);
            for j = i+1:length(xsel_conds)
                xrow = xsel_conds(j);
                xheight = (y_lim(2)-0.03) - (n * 0.02);
                text(-4, xheight, sprintf('%s vs. %s', xcond_name{xcol}, xcond_name{xrow}));
                plot(x_lim, [xheight xheight], ':', 'Color', [0.75 0.75 0.75])
                
                for xtr = 1:n_tc_trs
                    xpvalue = timecourse_cate.target.cate.t_test{xcate}.targ{xtr}.pvalue(xrow, xcol);
                    
                    if xpvalue <= xalpha
                        text(xtr-0.1, xheight, '*', 'FontSize', 20);
                    end
                    
                    plot([xtr xtr], y_lim, ':', 'Color', [0.75 0.75 0.75])
                end
                
                n = n + 1;
            end
        end
        
        %*************** std error-bar filling
        for xtarg = xsel_conds
            clear xerr
            xerr(1,:) = xmean{xtarg} - xstd{xtarg};
            xerr(2,:) = xmean{xtarg} + xstd{xtarg};
            
            for i = 1:2, fit_err(:,i) = interp1(x_tick, xerr(i,:), fitx,'spline'); end
            
            in_between = [fit_err(:,1)', fliplr(fit_err(:,2)')];
            h = fill([fitx fliplr(fitx)], in_between, xcond_color{xtarg});
            set(h,'facealpha', .4)
            
            plot(fitx, fit_err(:,1), '-', 'Color', xcond_color{xtarg})
            plot(fitx, fit_err(:,2), '-', 'Color', xcond_color{xtarg})
        end
        
        %*************** baseline 0
        plot(x_lim, [0 0], '--', 'Color', xbase_color)
        
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
        
        title(sprintf('%s: target evidence (N=%s)', ...
            category_names{xcate}, num2str(length(xsub_groups))));
        xlabel('Volume (tr)');
        ylabel('classifier evidence');
        
        h1 = text(1.5, y_lim(1) + 0.05, '(1) stim onset', 'Color', xonset_color);
        h2 = text(dur_stim + 1.5, y_lim(1) + 0.05, '(2) operation onset', 'Color', xonset_color);
        h3 = text(dur_stim + dur_manip + 1.5, y_lim(1) + 0.05, '(3) fixation onset', 'Color', xonset_color);
        set(h1,'Rotation', 90); set(h2,'Rotation', 90); set(h3,'Rotation', 90);
        
        text(on_stim - 1, y_lim(1) + 0.05, sprintf('(1) shift %str', num2str(args.shift_TRs)));
        text(on_operation - 1, y_lim(1) + 0.05, sprintf('(2) shift %str', num2str(args.shift_TRs)));
        text(on_fixation - 1, y_lim(1) + 0.05, sprintf('(3) shift %str', num2str(args.shift_TRs)));
        
    end
    
    %% *************** save fig
    fig_fname = fullfile(dirs.mvpa.group.out{xph}, ...
        sprintf('grp_plot_cate_based_target_wn_timecourse_%s', ...
        xsel, base_name));
    
    savefig(LogRegs_fig, sprintf('%s.fig', fig_fname));
    set(gcf,'PaperUnits','inches','PaperPosition', fig_rect/100)
    print('-djpeg', sprintf('%s.jpg',fig_fname), '-r100')
    saveas(LogRegs_fig, sprintf('%s.jpg',fig_fname), 'jpg')
    
%     close(LogRegs_fig);
    
end

end