function[] = grp_analysis_timecourse_workingmemory(args, grp_mvpaout, dirs)

%---------------------------------------------------------------------------
%*************** group: timecourse of the working memory (category|subcategory)
%---------------------------------------------------------------------------

xph = args.xphase;
fprintf('\n(+) timecourse of working memory: %s\n', args.phase_name{xph});

xsub_groups = args.filtered_subs;%args.g_sub;

%*************** checking version
assert(strcmp(version,'9.7.0.1190202 (R2019b)'),'[*] need higher version than R2017a');

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
xblk_alpha      = xalpha/n_stat_blks;

%*************** BASELINE CORRECTION
baseline_trs    = args.baseline_trs;

%*************** output basename
base_name       = sprintf('%s_n%s', args.analysis_basename, num2str(length(xsub_groups)));

%% ============= PLOT PARAMS
if strcmp(args.level, 'subcategory')
    for xcond = 1:n_condition
        
        if xcond==2, n_targ = 5; elseif xcond==3, n_targ = 4; else, n_targ = 3; end
        
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
        
        if xcond==2, n_targ = 3; else, n_targ = 2; end
        
        xcolor{xcond}{1}      = [238, 20, 91]/255;% targ
        xcolor{xcond}{n_targ} = [144, 144, 144]/255;% non_target
        
        if xcond==2
            xcolor{xcond}{2} = [0, 188, 182]/255;% newtarg  
        end
    end
end

on_stim      = args.shift_TRs + 1;
on_operation = args.shift_TRs + 1 + dur_stim;
on_fixation  = on_operation + dur_manip;

%*************** legend
if strcmp(args.level,'subcategory')
    for i = [1 4 5]
        legend_names{i} = {'target','nontargets','baseline'};
    end
    
    legend_names{2} = {'target','nontargets','new target','new nontargets','baseline'};
    legend_names{3} = {'target','nontargets','new target','baseline'};
    
elseif strcmp(args.level,'category')
    for i = [1 3:5]
        legend_names{i} = {'target','baseline'};
    end
    
    legend_names{2} = {'target','new target','baseline'};
    
end

xcond_color  = args.cond_color;
% xcate_color  = args.cate_color;

xbase_color  = args.base_color;
xonset_color = args.onset_color;

%% ============= PARSE TIMECOURSE 
fname = fullfile(dirs.mvpa.group.out{xph}, sprintf('grp_wmem_timecourse_blk%d_%s.mat', n_tr_blks, base_name));

if args.parse_timecourse
    %% ============= TIMECOURSE: RANDOM
    %*************** working memory contents
    % condition: 1_maintain,2_replace_category,3_replace_item,4_target_suppress,5_global_clear
    % xtr: 12 presentation + 5:9 fixation = 21 max
    
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
    clear timecourse
    
    for xcond = 1:n_condition
        if strcmp(args.level,'subcategory')
            if xcond==2, n_targ = 5; elseif xcond==3, n_targ = 4; else, n_targ = 3; end
        elseif strcmp(args.level,'category')
            if xcond==2, n_targ = 3; else, n_targ = 2; end
        end
        
        for xtr = 1:n_tc_trs
            for xtarg = 1:n_targ
                
                %*************** random effect
                xevidence = []; xevidence_sync = [];
                
                for xsub = xsub_groups
                    xevidence = horzcat(xevidence, ...
                        mean(grp_mvpaout{xsub}.decode.timecourse.condition{xcond}.target{xtarg}.evidence{xtr}));
                    
                    if xtr <= dur_sync
                        xevidence_sync = horzcat(xevidence_sync, ...
                            mean(grp_mvpaout{xsub}.decode.timecourse.sync.condition{xcond}.target{xtarg}.evidence{xtr}));
                    end
                end
                
                timecourse.evidence{xcond}.targ{xtarg}.tr{xtr}      = xevidence;
                timecourse.mean{xcond}(xtarg, xtr)                  = mean(xevidence);
                timecourse.se{xcond}(xtarg, xtr)                    = std(xevidence)/sqrt(length(xsub_groups));
                
                if xtr <= dur_sync
                    timecourse.sync.evidence{xcond}.targ{xtarg}.tr{xtr} = xevidence_sync;
                    timecourse.sync.mean{xcond}(xtarg, xtr)             = mean(xevidence_sync);
                    timecourse.sync.se{xcond}(xtarg, xtr)               = std(xevidence_sync)/sqrt(length(xsub_groups));
                end
            end
        end
    end
    
    %% ============= STATS: RANDOM
    for xcond = 1:n_condition
        if strcmp(args.level,'subcategory')
            if xcond==2, n_targ = 5; elseif xcond==3, n_targ = 4; else, n_targ = 3; end
        elseif strcmp(args.level,'category')
            if xcond==2, n_targ = 3; else, n_targ = 2; end
        end
        
        for xblk = 1:n_stat_blks
            
            timecourse.t_test{xcond}.blk{xblk}.pvalue = nan(n_targ, n_targ);
            
            it_trs = (1:n_tr_blks) + (n_tr_blks * (xblk-1));
            
            for xcol = 1:(n_targ-1)
                for xrow = (xcol+1):n_targ
                    
                    xtarg_col = []; xtarg_row = [];
                    
                    for xtr = it_trs
                        xtarg_col = vertcat(xtarg_col, timecourse.evidence{xcond}.targ{xcol}.tr{xtr});
                        xtarg_row = vertcat(xtarg_row, timecourse.evidence{xcond}.targ{xrow}.tr{xtr});
                    end
                    
                    %*************** ttest
                    [xhypothesis, xpvalue, ~, xstats] = ttest(mean(xtarg_col), mean(xtarg_row), 'Alpha', xblk_alpha);
                    
                    timecourse.t_test{xcond}.blk{xblk}.pvalue(xrow, xcol)     = xpvalue;
                    timecourse.t_test{xcond}.blk{xblk}.hypothesis(xrow, xcol) = xhypothesis;
                    timecourse.t_test{xcond}.blk{xblk}.stats{xrow, xcol}      = xstats;
                end
            end
            
            %*************** sync
            if xblk <= (dur_sync/n_tr_blks)
                
                timecourse.sync.t_test{xcond}.blk{xblk}.pvalue = nan(n_targ, n_targ);
                
                it_trs = (1:n_tr_blks) + (n_tr_blks * (xblk-1));
                
                for xcol = 1:(n_targ-1)
                    for xrow = (xcol+1):n_targ
                        
                        xtarg_col = []; xtarg_row = [];
                        
                        for xtr = it_trs
                            xtarg_col = vertcat(xtarg_col, timecourse.sync.evidence{xcond}.targ{xcol}.tr{xtr});
                            xtarg_row = vertcat(xtarg_row, timecourse.sync.evidence{xcond}.targ{xrow}.tr{xtr});
                        end
                        
                        %*************** ttest
                        [xhypothesis, xpvalue, ~, xstats] = ttest(mean(xtarg_col), mean(xtarg_row), 'Alpha', xblk_alpha);
                        
                        timecourse.t_test{xcond}.blk{xblk}.pvalue(xrow, xcol)     = xpvalue;
                        timecourse.t_test{xcond}.blk{xblk}.hypothesis(xrow, xcol) = xhypothesis;
                        timecourse.t_test{xcond}.blk{xblk}.stats{xrow, xcol}      = xstats;
                    end
                end
            end
        end
    end
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%============= BASELINE CORRECTION: RANDOM
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    clear timecourse_bcorr
    for xcond = 1:n_condition
        if strcmp(args.level,'subcategory')
            if xcond==2, n_targ = 5; elseif xcond==3, n_targ = 4; else, n_targ = 3; end
        elseif strcmp(args.level,'category')
            if xcond==2, n_targ = 3; else, n_targ = 2; end
        end
        
        %*************** baseline_trs: 6TRs
        clear mean_baseline
        xbase_evidence = zeros(n_targ, baseline_trs);
        for xtarg = 1:n_targ
            for xtr = 1:baseline_trs
                xbase_evidence(xtarg, xtr) = ...
                    mean(timecourse.evidence{xcond}.targ{xtarg}.tr{xtr});
            end
            
            mean_baseline(xtarg) = mean(xbase_evidence(xtarg,:));
        end
        
        baseline_corr{xcond} = mean_baseline;
        
        %*************** corrected evidence
        for xtr = 1:n_tc_trs
            for xtarg = 1:n_targ
                clear xcorrect_evidence
                xcorrect_evidence = timecourse.evidence{xcond}.targ{xtarg}.tr{xtr} - ...
                    mean_baseline(xtarg);
                    
                timecourse_bcorr.evidence{xcond}.targ{xtarg}.tr{xtr} = xcorrect_evidence;
                timecourse_bcorr.mean{xcond}(xtarg, xtr)             = mean(xcorrect_evidence);
                timecourse_bcorr.se{xcond}(xtarg, xtr)               = std(xcorrect_evidence)/sqrt(length(xsub_groups));
                
                %*************** sync
                clear xcorrect_evidence
                if xtr <= dur_sync
                    xcorrect_evidence = timecourse.sync.evidence{xcond}.targ{xtarg}.tr{xtr} - ...
                        mean_baseline(xtarg);
                    
                    timecourse_bcorr.sync.evidence{xcond}.targ{xtarg}.tr{xtr} = xcorrect_evidence;
                    timecourse_bcorr.sync.mean{xcond}(xtarg, xtr)             = mean(xcorrect_evidence);
                    timecourse_bcorr.sync.se{xcond}(xtarg, xtr)               = std(xcorrect_evidence)/sqrt(length(xsub_groups));
                end
            end
        end
    end
   
    %% ============= STATS: RANDOM_BASELINE CORRECTED
    for xcond = 1:n_condition
        if strcmp(args.level,'subcategory')
            if xcond==2, n_targ = 5; elseif xcond==3, n_targ = 4; else, n_targ = 3; end
        elseif strcmp(args.level,'category')
            if xcond==2, n_targ = 3; else, n_targ = 2; end
        end
        
        for xblk = 1:n_stat_blks
            
            timecourse_bcorr.t_test{xcond}.blk{xblk}.pvalue = nan(n_targ, n_targ);
            
            it_trs = (1:n_tr_blks) + (n_tr_blks * (xblk-1));
            
            for xcol = 1:(n_targ-1)
                for xrow = (xcol+1):n_targ
                    
                    xtarg_col = []; xtarg_row = [];
                    
                    for xtr = it_trs
                        xtarg_col = vertcat(xtarg_col, timecourse_bcorr.evidence{xcond}.targ{xcol}.tr{xtr});
                        xtarg_row = vertcat(xtarg_row, timecourse_bcorr.evidence{xcond}.targ{xrow}.tr{xtr});
                    end
                    
                    %*************** ttest
                    [xhypothesis, xpvalue, ~, xstats] = ttest(mean(xtarg_col), mean(xtarg_row), 'Alpha', xblk_alpha);
                    
                    timecourse_bcorr.t_test{xcond}.blk{xblk}.pvalue(xrow, xcol)     = xpvalue;
                    timecourse_bcorr.t_test{xcond}.blk{xblk}.hypothesis(xrow, xcol) = xhypothesis;
                    timecourse_bcorr.t_test{xcond}.blk{xblk}.stats{xrow, xcol}      = xstats;
                end
            end
            
            %*************** sync
            if xblk <= (dur_sync/n_tr_blks)
                
                timecourse_bcorr.sync.t_test{xcond}.blk{xblk}.pvalue = nan(n_targ, n_targ);
                
                it_trs = (1:n_tr_blks) + (n_tr_blks * (xblk-1));
                
                for xcol = 1:(n_targ-1)
                    for xrow = (xcol+1):n_targ
                        
                        xtarg_col = []; xtarg_row = [];
                        
                        for xtr = it_trs
                            xtarg_col = vertcat(xtarg_col, timecourse_bcorr.sync.evidence{xcond}.targ{xcol}.tr{xtr});
                            xtarg_row = vertcat(xtarg_row, timecourse_bcorr.sync.evidence{xcond}.targ{xrow}.tr{xtr});
                        end
                        
                        %*************** ttest
                        [xhypothesis, xpvalue, ~, xstats] = ttest(mean(xtarg_col), mean(xtarg_row), 'Alpha', xblk_alpha);
                        
                        timecourse_bcorr.sync.t_test{xcond}.blk{xblk}.pvalue(xrow, xcol)     = xpvalue;
                        timecourse_bcorr.sync.t_test{xcond}.blk{xblk}.hypothesis(xrow, xcol) = xhypothesis;
                        timecourse_bcorr.sync.t_test{xcond}.blk{xblk}.stats{xrow, xcol}      = xstats;
                    end
                end
            end
        end
    end
   
    %% ============= BASELINE CORRECTION: TARGET-ITEM EVIDENCE ONLY
    for xcond = 1:n_condition
        for xtr = 1:n_tc_trs
            xevidence = timecourse_bcorr.evidence{xcond}.targ{1}.tr{xtr};
            
            timecourse_bcorr.targets.evidence{xcond}.tr{xtr} = xevidence;
            timecourse_bcorr.targets.mean{xcond}(xtr)        = mean(xevidence);
            timecourse_bcorr.targets.se{xcond}(xtr)          = std(xevidence)/sqrt(length(xsub_groups));
            
            %*************** sync
            if xtr <= dur_sync
                xevidence = timecourse_bcorr.sync.evidence{xcond}.targ{1}.tr{xtr};
                
                timecourse_bcorr.targets.sync.evidence{xcond}.tr{xtr} = xevidence;
                timecourse_bcorr.targets.sync.mean{xcond}(xtr)        = mean(xevidence);
                timecourse_bcorr.targets.sync.se{xcond}(xtr)          = std(xevidence)/sqrt(length(xsub_groups));
            end
        end
    end
    
    %% ============= STATS: BASELINE CORRECTION: TARGET-ITEM EVIDENCE ONLY
    for xblk = 1:n_stat_blks
        
        timecourse_bcorr.targets.t_test.blk{xblk}.pvalue = nan(n_condition-1, n_condition-1);
        
        it_trs = (1:n_tr_blks) + (n_tr_blks * (xblk-1));
        
        for xcol = 1:n_condition-1
            for xrow = xcol+1:n_condition
                
                xtarg_col = []; xtarg_row = [];
                
                for xtr = it_trs
                    xtarg_col = vertcat(xtarg_col, timecourse_bcorr.targets.evidence{xcol}.tr{xtr});
                    xtarg_row = vertcat(xtarg_row, timecourse_bcorr.targets.evidence{xrow}.tr{xtr});
                end
                
                %*************** ttest
                [xhypothesis, xpvalue, ~, xstats] = ttest(mean(xtarg_col), mean(xtarg_row), 'Alpha', xblk_alpha);
                
                timecourse_bcorr.targets.t_test.blk{xblk}.pvalue(xrow, xcol)     = xpvalue;
                timecourse_bcorr.targets.t_test.blk{xblk}.hypothesis(xrow, xcol) = xhypothesis;
                timecourse_bcorr.targets.t_test.blk{xblk}.stats{xrow, xcol}      = xstats;
            end
        end
    end
    
    %*************** sync
    for xblk = 1:(dur_sync/n_tr_blks)
        
        timecourse_bcorr.targets.sync.t_test.blk{xblk}.pvalue = nan(n_condition-1, n_condition-1);
        
        it_trs = (1:n_tr_blks) + (n_tr_blks * (xblk-1));
        
        for xcol = 1:n_condition-1
            for xrow = xcol+1:n_condition
                
                xtarg_col = []; xtarg_row = [];
                
                for xtr = it_trs
                    xtarg_col = vertcat(xtarg_col, timecourse_bcorr.targets.sync.evidence{xcol}.tr{xtr});
                    xtarg_row = vertcat(xtarg_row, timecourse_bcorr.targets.sync.evidence{xrow}.tr{xtr});
                end
                
                %*************** ttest
                [xhypothesis, xpvalue, ~, xstats] = ttest(mean(xtarg_col), mean(xtarg_row), 'Alpha', xblk_alpha);
                
                timecourse_bcorr.targets.sync.t_test.blk{xblk}.pvalue(xrow, xcol)     = xpvalue;
                timecourse_bcorr.targets.sync.t_test.blk{xblk}.hypothesis(xrow, xcol) = xhypothesis;
                timecourse_bcorr.targets.sync.t_test.blk{xblk}.stats{xrow, xcol}      = xstats;
            end
        end
    end
    
    %% ============= STATS: ANOVA: 5/4 operations
    sel_cond{1} = 1:n_condition;%all 5
    sel_cond{2} = [1 2 4 5];
        
    for xsel = 1:length(sel_cond)
        for xblk = 1:n_stat_blks
            clear xevidence xdata xheader xtable
            
            %*************** collecting evidence
            xevidence  = [];
            
            it_trs = (1:n_tr_blks) + (n_tr_blks * (xblk-1));
            
            for xcond = 1:length(sel_cond{xsel})
                
                it_cond = sel_cond{xsel}(xcond);
                
                t_evidence = [];
                for xtr = it_trs
                    t_evidence = vertcat(t_evidence, timecourse_bcorr.targets.evidence{it_cond}.tr{xtr});
                end
                
                xevidence(:, xcond) = mean(t_evidence)';
                
                xheader{xcond+1} = condition_names{it_cond};
            end
            
            %*************** data table for ANOVA
            xdata(:,1) = xsub_groups';
            xdata(:, 2:size(xevidence, 2)+1) = xevidence;
            xheader{1} = 'SN';
            xtable     = array2table(xdata, 'VariableNames', xheader);
            
            xanova     = {'target_evidence'};
            xmeasures  = table((1:length(sel_cond{xsel}))','VariableNames', xanova);
            xrepmeas   = fitrm(xtable,'Maintain-Clear~1','WithinDesign',xmeasures);
            xanova_out = ranova(xrepmeas);
            
            var_names  = [xanova xanova_out.Properties.VariableNames];
            row_names  = xanova_out.Properties.RowNames;
            
            tmp_val    = [row_names num2cell(xanova_out.Variables)];
            xanva_table = cell2table(tmp_val, 'VariableNames', var_names);
            
            xtx_out = sprintf('ANOVA: F(%s, %s)=%4.4f, p=%1.4f', ...
                num2str(xanova_out.DF(1)), num2str(xanova_out.DF(2)),...
                xanova_out.F(1), xanova_out.pValue(1));
            
            timecourse_bcorr.targets.anova{xsel}.blk{xblk}.pvalue = xanova_out.pValue(1);
            timecourse_bcorr.targets.anova{xsel}.blk{xblk}.table  = xanva_table;
            timecourse_bcorr.targets.anova{xsel}.blk{xblk}.stats  = xtx_out;
        end
        
        %*************** sync
        for xblk = 1:(dur_sync/n_tr_blks)
            clear xevidence xdata xheader xtable
            
            %*************** collecting evidence
            xevidence  = [];
            
            it_trs = (1:n_tr_blks) + (n_tr_blks * (xblk-1));
            
            for xcond = 1:length(sel_cond{xsel})
                
                it_cond = sel_cond{xsel}(xcond);
                
                t_evidence = [];
                for xtr = it_trs
                    t_evidence = vertcat(t_evidence, timecourse_bcorr.targets.sync.evidence{it_cond}.tr{xtr});
                end
                
                xevidence(:, xcond) = mean(t_evidence)';
                
                xheader{xcond+1} = condition_names{it_cond};
            end
            
            %*************** data table for ANOVA
            xdata(:,1) = xsub_groups';
            xdata(:, 2:size(xevidence, 2)+1) = xevidence;
            xheader{1} = 'SN';
            xtable     = array2table(xdata, 'VariableNames', xheader);
            
            xanova     = {'target_evidence'};
            xmeasures  = table((1:length(sel_cond{xsel}))','VariableNames', xanova);
            xrepmeas   = fitrm(xtable,'Maintain-Clear~1','WithinDesign',xmeasures);
            xanova_out = ranova(xrepmeas);
            
            var_names  = [xanova xanova_out.Properties.VariableNames];
            row_names  = xanova_out.Properties.RowNames;
            
            tmp_val    = [row_names num2cell(xanova_out.Variables)];
            xanva_table = cell2table(tmp_val, 'VariableNames', var_names);
            
            xtx_out = sprintf('ANOVA: F(%s, %s)=%4.4f, P=%1.4f', ...
                num2str(xanova_out.DF(1)), num2str(xanova_out.DF(2)),...
                xanova_out.F(1), xanova_out.pValue(1));
            
            timecourse_bcorr.targets.sync.anova{xsel}.blk{xblk}.pvalue = xanova_out.pValue(1);
            timecourse_bcorr.targets.sync.anova{xsel}.blk{xblk}.table  = xanva_table;
            timecourse_bcorr.targets.sync.anova{xsel}.blk{xblk}.stats  = xtx_out;
        end
    end
    
    %% ============= REPLACE: BIN BASED
    % baseline corrected
    %--------------- category level
    % 2. replace category    : {1} target, {2} new category, {3} nontargets
    %--------------- subcategory level
    % 2. replace category    : {1} subtarget, {2} nonsubtarget(2), {3} new_subtarget, {4} new_nonsubtarget(2), {5} mean(nontarget(3))
    % 3. replace subcategory : {1} subtarget, {2} nonsubtarget(1), {3} new_subtarget(1), {4} mean(nontarget(6))
    
    if strcmp(args.level,'category')
        it_conds = 2; n_targ = 3;
    elseif strcmp(args.level,'subcategory')
        it_conds = [2 3]; 
    end
    
    for xcond = it_conds
        
        if strcmp(args.level,'subcategory')
            if xcond==2, n_targ = 5; elseif xcond==3, n_targ = 4; end
        end
        
        for xbin = 1:2
            
            it_trs = args.bin_windows(xbin,1):args.bin_windows(xbin,2);
            
            for xtarg = 1:n_targ
                
                %*************** random effect
                xevidence = [];
                
                for xsub = xsub_groups
                    
                    t_evidence = [];
                    
                    for xtr = it_trs
                        t_evidence = horzcat(t_evidence, ...
                            mean(grp_mvpaout{xsub}.decode.timecourse.condition{xcond}.target{xtarg}.evidence{xtr}) - ...
                            baseline_corr{xcond}(xtarg));
                    end
                    
                    xevidence = horzcat(xevidence, mean(t_evidence));
                end
                
                timebin.evidence{xcond}.targ{xtarg}.bin{xbin} = xevidence;
                timebin.mean{xcond}(xtarg, xbin)              = mean(xevidence);
                timebin.se{xcond}(xtarg, xbin)                = std(xevidence)/sqrt(length(xsub_groups));
                
            end
        end
    end
    
    %% *************** STATS
    for xcond = it_conds
        if strcmp(args.level,'subcategory')
            if xcond==2, n_targ = 5; elseif xcond==3, n_targ = 4; end
        end
        
        for xbin = 1:2
            
            timebin.t_test{xcond}.bin{xbin}.pvalue = nan(n_targ, n_targ);
            
            for xcol = 1:(n_targ-1)
                for xrow = (xcol+1):n_targ
                    
                    xtarg_col = timebin.evidence{xcond}.targ{xcol}.bin{xbin};
                    xtarg_row = timebin.evidence{xcond}.targ{xrow}.bin{xbin};
                    
                    %*************** ttest
                    [~, xpvalue, ~, xstats] = ttest(xtarg_col, xtarg_row, 'Alpha', 0.05);
                    
                    timebin.t_test{xcond}.bin{xbin}.pvalue(xrow, xcol) = xpvalue;
                    timebin.t_test{xcond}.bin{xbin}.stats{xrow, xcol}  = xstats;
                end
            end
        end
    end
    
    %% ============= SAVE
    grp_timecourse.tbin     = timebin;%baseline corrected
    grp_timecourse.tc       = timecourse;
    grp_timecourse.tc_bcorr = timecourse_bcorr; %#ok<*STRNU>
    
    save(fname, 'grp_timecourse');
    
else
    load(fname);%'grp_timecourse'
    
    timecourse       = grp_timecourse.tc; %#ok<*NASGU,*NODEF>
    timecourse_bcorr = grp_timecourse.tc_bcorr;
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BAR PLOT FOR REPLACE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(args.level,'category')
    it_conds = 2; n_targ = 3;
elseif strcmp(args.level,'subcategory')
    it_conds = [2 3];
end

y_lim    = [-0.3 0.6];
fig_rect = [0 0 1000 500*length(it_conds)];

xfig = figure;
set(xfig, 'Position', fig_rect)

for i = 1:length(it_conds)
        
    xcond = it_conds(i);
    
    if strcmp(args.level,'subcategory')
        if xcond==2, n_targ = 5; elseif xcond==3, n_targ = 4; end
    end
    
    for xbin = 1:2
        clear b
        
        subplot(length(it_conds), 2, xbin + 2*(i-1))
        
        xmean = timebin.mean{xcond}(:, xbin);
        xse   = timebin.se{xcond}(:, xbin);
       
        for xtarg = 1:n_targ
            b{xtarg} = bar(xtarg, xmean(xtarg)); hold on
            set(b{xtarg}, 'facecolor', xcolor{xcond}{xtarg});
        end
        
        for xtarg = 1:n_targ
            errorbar(xtarg, xmean(xtarg), xse(xtarg),'k.')
        end
        
        %*************** ttest
        n = 1;
        for xcol = 1:(n_targ-1)
            for xrow = (xcol+1):n_targ
                xpvalue = timebin.t_test{xcond}.bin{xbin}.pvalue(xrow, xcol);
                xstats  = timebin.t_test{xcond}.bin{xbin}.stats{xrow, xcol};
                
                if xpvalue <= xalpha
                    yy_sig = y_lim(2) - (0.05 * n);
                    plot([xcol xrow], [yy_sig yy_sig], '-k')
                    n = n + 1;
                    
                    if (xpvalue <= xalpha) && (xpvalue > (xalpha/5))
                        xsig = '*';
                    elseif (xpvalue <= (xalpha/5)) && (xpvalue > (xalpha/50))
                        xsig = '**';
                    elseif xpvalue <= (xalpha/50)%0.001
                        xsig = '***';
                    end
                    
                    text(xcol+((xrow-xcol)/2), yy_sig, xsig, 'FontSize', 20);
                    text(xrow, yy_sig + 0.01, ...
                        sprintf('T(%s)=%4.4f, P=%4.4f', ...
                        num2str(xstats.df),xstats.tstat,  xpvalue), 'FontSize', 10);
                end
            end
        end
        
        %*************** legend
        xlegend        = legend_names{xcond};
        lg             = legend(xlegend);
        lg.Location    = 'southeastoutside';%'SouthEast';
        legend(xlegend,'AutoUpdate','off')
        
        set(gca,'xlim',[1 n_targ]+[-0.5 0.5],'ylim', y_lim);
        set(gca,'XTick', 1:n_targ, 'YTick', y_lim(1):0.1:y_lim(end));
        set(gca,'XTickLabel', legend_names{xcond})
        
        title(sprintf('%s: %s: %s-%s', ...
            condition_names{xcond}, args.level, ...
            num2str(args.bin_windows(xbin,1)), num2str(args.bin_windows(xbin,2))));
        ylabel('classifier evidence');
    end
end

%*************** save fig
fig_fname = fullfile(dirs.mvpa.group.out{xph}, ...
    sprintf('grp_bar_wn_bin_%s_blk%d_%s', ...
    condition_names{xcond}, n_tr_blks,base_name));

savefig(xfig, sprintf('%s.fig', fig_fname));
set(gcf,'PaperUnits','inches','PaperPosition', fig_rect/100)
print('-dpng', '-r300', sprintf('%s.png',fig_fname))
saveas(xfig, sprintf('%s.png',fig_fname), 'png')

close(xfig);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TIMECOURSE PLOTS: BASELINE CORRECTED 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x_tick       = 1:n_tc_trs;
x_ticklable  = 1:n_tc_trs;
x_lim        = [x_tick(1)-1 x_tick(end)+1];% + [-0.5 0.5];
if strcmp(args.level, 'category'), y_lim = [-0.3 0.6];
else,                              y_lim = [-0.2 0.4]; end
fig_rect     = [0 0 1200 800]; 

%% ============= CONTINUOUS PLOT: BASELINE CORRECTED RANDOM EFFECT
% target vs. nontarget / target vs. newtarget vs. nontarget

for xcond = 1:n_condition
    clear xmean xstd fity tt pp
    
    xfig = figure;
    set(xfig, 'Position', fig_rect)
    
    if strcmp(args.level,'subcategory')
        if xcond==2, n_targ = 5; elseif xcond==3, n_targ = 4; else, n_targ = 3; end
    elseif strcmp(args.level,'category')
        if xcond==2, n_targ = 3; else, n_targ = 2; end
    end
    
    %*************** mean line plots
    fitx = linspace(1, n_tc_trs, n_tc_trs*10);
    
    for xtarg = 1:n_targ
        xmean{xtarg} = timecourse_bcorr.mean{xcond}(xtarg, x_tick);
        xstd{xtarg}  = timecourse_bcorr.se{xcond}(xtarg, x_tick); %#ok<*AGROW>
        
        fity{xtarg}  = interp1(x_tick, xmean{xtarg}, fitx,'spline');
        
        plot(fitx, fity{xtarg}, '-','Color', xcolor{xcond}{xtarg}, 'LineWidth', 2); hold on;
    end
    
    %*************** legend
    xlegend        = legend_names{xcond};
    lg             = legend(xlegend);
    lg.Location    = 'SouthEast';
    legend(xlegend,'AutoUpdate','off')
    
    %*************** stats stars
    n = 0; 
    for xcol = 1:(n_targ-1)
        for xrow = (xcol+1):n_targ
            for ss = 1:2, tt{ss} = []; pp{ss} = []; end
            
            xheight = (y_lim(2)-0.05) - (n * 0.02);
            text(-3, xheight, sprintf('%s vs. %s', legend_names{xcond}{xcol}, ...
                legend_names{xcond}{xrow}));
            plot(x_lim, [xheight xheight], ':', 'Color', [0.75 0.75 0.75])
            
            for xblk = 1:n_stat_blks
                
                xtr = (xblk-1) * n_tr_blks;
                plot([xtr xtr]+0.5, y_lim, '--', 'Color', [0.75 0.75 0.75])
                
                xpvalue = timecourse_bcorr.t_test{xcond}.blk{xblk}.pvalue(xrow, xcol);
                xstats  = timecourse_bcorr.t_test{xcond}.blk{xblk}.stats{xrow, xcol};
                
                if (xpvalue <= (0.1/n_stat_blks)) && (xpvalue > (xblk_alpha))
                    xsig = '+';
                elseif (xpvalue <= xblk_alpha) && (xpvalue > (xblk_alpha/5))
                    xsig = '*';
                elseif (xpvalue <= (xblk_alpha/5)) && (xpvalue > (xblk_alpha/50))
                    xsig = '**';
                elseif xpvalue <= (xblk_alpha/50)%0.001
                    xsig = '***';
                end
                
                if xpvalue <= (0.1/n_stat_blks)
                    text(xtr + (n_tr_blks/2)+0.3, xheight, xsig, 'FontSize', 20);
                    xx = xtr + 0.7 + (0.4 * n); 
                    
                    h = text(xx, y_lim(end) - 0.3, sprintf('T=%4.4f, P=%4.4f', ...
                        xstats.tstat,  xpvalue), 'FontSize', 8);
                    set(h,'Rotation',90);
                    
                    if xstats.tstat>0, ss = 1; else, ss = 2; end
                    tt{ss} = horzcat(tt{ss}, xstats.tstat);
                    pp{ss} = horzcat(pp{ss}, xpvalue);
                end
            end
            n = n + 1;
            
            if ~isempty(tt{1})
                text(2, 0.2 - (n * 0.05), ...
                    sprintf('Ts(%s)>%4.4f, Ps<%4.4f', num2str(xstats.df), min(tt{1}), max(pp{1})), 'FontSize', 10);
            end
            if  ~isempty(tt{2})
                text(10, 0.2 - (n * 0.05), ...
                    sprintf('Ts(%s)<%4.4f, Ps<%4.4f', num2str(xstats.df), max(tt{2}), max(pp{2})), 'FontSize', 10);
            end
        end
    end
    
    %*************** std error-bar filling
    for xtarg = 1:n_targ
        clear xerr fit_err
        xerr(1,:) = xmean{xtarg} - xstd{xtarg};
        xerr(2,:) = xmean{xtarg} + xstd{xtarg};
        
        for i = 1:2, fit_err(:,i) = interp1(x_tick, xerr(i,:), fitx,'spline'); end
        
        in_between = [fit_err(:,1)', fliplr(fit_err(:,2)')];
        h = fill([fitx fliplr(fitx)], in_between, xcolor{xcond}{xtarg});
        h.FaceAlpha = .2;
        h.EdgeAlpha = .2;
        h.EdgeColor = xcolor{xcond}{xtarg};
    end
    
    %*************** real onset lines
    plot([dur_stim dur_stim]+1, y_lim, '-', 'Color', xonset_color)
    plot([(dur_stim + dur_manip) (dur_stim + dur_manip)]+1, y_lim, '-', 'Color', xonset_color)
    
    set(gca,'xlim', x_lim, 'ylim', y_lim);
    set(gca,'XTick', x_tick, 'YTick', y_lim(1):0.1:y_lim(end))
    set(gca,'XTickLabel', x_ticklable)
    
    title(sprintf('%s classifier: %s (N=%s, p<%1.4f)', ...
        args.level, condition_names{xcond}, num2str(length(xsub_groups)), xblk_alpha));
    xlabel('Volume (tr)');
    ylabel('classifier evidence');
    
    h1 = text(1.5, y_lim(1)+0.05, '(1) stim onset', 'Color', xonset_color);
    h2 = text(dur_stim + 1.5, y_lim(1)+0.05, '(2) operation onset', 'Color', xonset_color);
    h3 = text(dur_stim + dur_manip + 1.5, y_lim(1)+0.05, '(3) fixation onset', 'Color', xonset_color);
    set(h1,'Rotation', 90); set(h2,'Rotation', 90); set(h3,'Rotation', 90);
    
    %% *************** save fig
    fig_fname = fullfile(dirs.mvpa.group.out{xph}, ...
        sprintf('grp_plot_wn_tc_%s_blk%d_%s', ...
        condition_names{xcond}, n_tr_blks,base_name));
    
    savefig(xfig, sprintf('%s.fig', fig_fname));
    set(gcf,'PaperUnits','inches','PaperPosition', fig_rect/100)
    print('-dpng', '-r300', sprintf('%s.png',fig_fname))
    saveas(xfig, sprintf('%s.png',fig_fname), 'png')
    
    close(xfig);
    
end

%% ============= BLOCKED CONTINUOUS PLOT: RAW TARGET RANDOM EFFECT 
%*************** ANOVA per block
%*************** evidence = target 
sel_conds{1} = 1:n_condition;
sel_conds{2} = [1 2 4 5];
sel_conds{3} = [1 4];
sel_conds{4} = [2 3];
sel_conds{5} = [1 2];
sel_conds{6} = [4 5];
sel_conds{7} = [2 4];
sel_conds{8} = [2 5];
    
for xfill = 1%0:1
    for xsel = 1:length(sel_conds)
        xsel_conds = sel_conds{xsel};
        
        xfig = figure;
        set(xfig, 'Position', fig_rect)
        
        clear xmean xstd fity tt pp
        
        %*************** mean line plots
        fitx = linspace(1, n_tc_trs, n_tc_trs*10);
        
        for xcond = xsel_conds
            xmean{xcond} = timecourse_bcorr.targets.mean{xcond}(x_tick);
            xstd{xcond}  = timecourse_bcorr.targets.se{xcond}(x_tick); %#ok<*AGROW>
            
            fity{xcond}  = interp1(x_tick, xmean{xcond}, fitx,'spline');
            
            plot(fitx, fity{xcond}, '-','Color', xcond_color{xcond}, 'LineWidth', line_w); hold on;
        end
        
        %*************** legend
        clear xlegend
        for i=1:length(xsel_conds); xlegend{i} = condition_names{xsel_conds(i)}; end
        lg             = legend(xlegend);
        lg.Location    = 'SouthEast';
        legend(xlegend,'AutoUpdate','off')
        
        %*************** stats stars
        n = 0;
        for i = 1:(length(xsel_conds)-1)
            xcol = xsel_conds(i);
            for j = i+1:length(xsel_conds)
                for ss = 1:2, tt{ss} = []; pp{ss} = []; end
                
                xrow    = xsel_conds(j);
                xheight = (y_lim(2)-0.05) - (n * 0.02);
                text(-2, xheight, sprintf('%d vs. %d', xrow, xcol));
                plot(x_lim, [xheight xheight], ':', 'Color', [0.75 0.75 0.75])
                
                for xblk = 1:n_stat_blks
                    
                    xtr = (xblk-1) * n_tr_blks;
                    plot([xtr xtr]+0.5, y_lim, '--', 'Color', [0.75 0.75 0.75])
                    
                    xpvalue = timecourse_bcorr.targets.t_test.blk{xblk}.pvalue(xrow, xcol);
                    xstats  = timecourse_bcorr.targets.t_test.blk{xblk}.stats{xrow, xcol};
                    
                    if (xpvalue <= (0.1/n_stat_blks)) && (xpvalue > (xblk_alpha))
                        xsig = '+';
                    elseif (xpvalue <= xblk_alpha) && (xpvalue > (xblk_alpha/5))
                        xsig = '*';
                    elseif (xpvalue <= (xblk_alpha/5)) && (xpvalue > (xblk_alpha/50))
                        xsig = '**';
                    elseif xpvalue <= (xblk_alpha/50)%0.001
                        xsig = '***';
                    end
                    
                    if xpvalue <= (0.1/n_stat_blks)
                        text(xtr + (n_tr_blks/2)+0.3, xheight, xsig, 'FontSize', 20);
                        xx = xtr + 0.7 + (0.4 * n);
                        
                        if length(xsel_conds) < 3
                            h = text(xx, y_lim(end) - 0.3, sprintf('T=%4.4f, P=%4.4f', ...
                                xstats.tstat,  xpvalue), 'FontSize', 8);
                            set(h,'Rotation',90);
                            
                            if xstats.tstat>0, ss = 1; else, ss = 2; end
                            tt{ss} = horzcat(tt{ss}, xstats.tstat);
                            pp{ss} = horzcat(pp{ss}, xpvalue);
                        end
                    end
                end
                n = n + 1;
                
                if length(xsel_conds) < 3
                    if ~isempty(tt{1})
                        text(2, 0.2 - (n * 0.05), ...
                            sprintf('Ts(%s)>%4.4f, Ps<%4.4f', num2str(xstats.df), min(tt{1}), max(pp{1})), 'FontSize', 10);
                    end
                    if  ~isempty(tt{2})
                        text(10, 0.2 - (n * 0.05), ...
                            sprintf('Ts(%s)<%4.4f, Ps<%4.4f', num2str(xstats.df), max(tt{2}), max(pp{2})), 'FontSize', 10);
                    end
                end
            end
        end
        
        %*************** anova
        if (xsel == 1) || (xsel == 2)
            xheight = y_lim(1) + 0.15;
            text(-2, xheight, sprintf('ANOVA: %d operations', length(length(xsel_conds))));
            
            for xblk = 1:n_stat_blks
                
                xpvalue = timecourse_bcorr.targets.anova{xsel}.blk{xblk}.pvalue;
                xstats  = timecourse_bcorr.targets.anova{xsel}.blk{xblk}.stats;
                xtr     = (xblk-1) * n_tr_blks;
                
                if (xpvalue <= (0.1/n_stat_blks)) && (xpvalue > (xblk_alpha))
                    xsig = '+';
                elseif (xpvalue <= xblk_alpha) && (xpvalue > (xblk_alpha/5))
                    xsig = '*';
                elseif (xpvalue <= (xblk_alpha/5)) && (xpvalue > (xblk_alpha/50))
                    xsig = '**';
                elseif xpvalue <= (xblk_alpha/50)%0.001
                    xsig = '***';
                end
                
                if xpvalue <= (0.1/n_stat_blks)
                    text(xtr+(n_tr_blks/2)+0.3, xheight, xsig, 'FontSize', 20);
                    h = text(xtr+(n_tr_blks/2)+0.3, -0.1, xstats, 'FontSize', 7);
                    set(h,'Rotation',90);
                end
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
                h.FaceAlpha = .2;
                h.EdgeAlpha = .2;
                h.EdgeColor = xcond_color{xcond};
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
        
        title(sprintf('%s classifier: targets (N=%s, p<%1.4f, TR/blk=%s)', ...
            args.level,num2str(length(xsub_groups)), xblk_alpha, num2str(n_tr_blks)));
        xlabel('Volume (tr)');
        ylabel('classifier evidence (target)');
        
        %*************** real onset
        h1 = text(1.5, y_lim(1)+0.01, '(1) stim onset', 'Color', xonset_color);
        h2 = text(dur_stim + 1.5, y_lim(1)+0.01, '(2) operation onset', 'Color', xonset_color);
        h3 = text(dur_stim + dur_manip + 1.5, y_lim(1)+0.01, '(3) fixation onset', 'Color', xonset_color);
        set(h1,'Rotation', 90); set(h2,'Rotation', 90); set(h3,'Rotation', 90);
        
        %% *************** save fig
        if xfill
            fig_fname = fullfile(dirs.mvpa.group.out{xph}, ...
                sprintf('grp_plot_wn_targets_sel%d_blk%d_%s', xsel, n_tr_blks, base_name));
        else
            fig_fname = fullfile(dirs.mvpa.group.out{xph}, ...
                sprintf('grp_plot_wn_targets_line_sel%d_blk%d_%s', xsel, n_tr_blks, base_name));
        end
        
        savefig(xfig, sprintf('%s.fig', fig_fname));
        set(gcf,'PaperUnits','inches','PaperPosition', fig_rect/100)
        print('-dpng', '-r300', sprintf('%s.png',fig_fname))
        saveas(xfig, sprintf('%s.png',fig_fname), 'png')
        
        close(xfig);
        
    end
end

%% ============= SYNCH: BLOCKED CONTINUOUS PLOT: RAW TARGET RANDOM EFFECT 
%*************** ANOVA per block
%*************** evidence = target 
fig_rect_sync = [0 0 500 800];
x_tick_sync   = 1:dur_sync;
x_lim_sync    = [0 dur_sync+1];
    
for xfill = 1%0:1
    for xsel = 1:length(sel_conds)
        xsel_conds = sel_conds{xsel};
        
        xfig = figure;
        set(xfig, 'Position', fig_rect_sync)
        
        clear xmean xstd fity tt pp
        
        %*************** mean line plots
        fitx = linspace(1, dur_sync, dur_sync*10);
        
        for xcond = xsel_conds
            xmean{xcond} = timecourse_bcorr.targets.sync.mean{xcond}(x_tick_sync);
            xstd{xcond}  = timecourse_bcorr.targets.sync.se{xcond}(x_tick_sync); %#ok<*AGROW>
            
            fity{xcond}  = interp1(x_tick_sync, xmean{xcond}, fitx,'spline');
            
            plot(fitx, fity{xcond}, '-','Color', xcond_color{xcond}, 'LineWidth', line_w); hold on;
        end
        
        %*************** legend
        clear xlegend
        for i=1:length(xsel_conds); xlegend{i} = condition_names{xsel_conds(i)}; end
        lg             = legend(xlegend);
        lg.Location    = 'SouthEast';
        legend(xlegend,'AutoUpdate','off')
        
        %*************** stats stars
        n = 0;
        for i = 1:(length(xsel_conds)-1)
            xcol = xsel_conds(i);
            for j = i+1:length(xsel_conds)
                xrow    = xsel_conds(j);
                
                for ss = 1:2, tt{ss} = []; pp{ss} = []; end
                
                xheight = (y_lim(2)-0.05) - (n * 0.015);
                text(-0.5, xheight, sprintf('%d vs. %d', xrow, xcol));
                plot(x_lim_sync, [xheight xheight], ':', 'Color', [0.75 0.75 0.75])
                
                for xblk = 1:(dur_sync/n_tr_blks)
                    
                    xtr = (xblk-1) * n_tr_blks;
                    plot([xtr xtr]+0.5, y_lim, '--', 'Color', [0.75 0.75 0.75])
                    
                    xpvalue = timecourse_bcorr.targets.sync.t_test.blk{xblk}.pvalue(xrow, xcol);
                    xstats  = timecourse_bcorr.targets.sync.t_test.blk{xblk}.stats{xrow, xcol};
                    
                    if (xpvalue <= (0.1/n_stat_blks)) && (xpvalue > (xblk_alpha))
                        xsig = '+';
                    elseif (xpvalue <= xblk_alpha) && (xpvalue > (xblk_alpha/5))
                        xsig = '*';
                    elseif (xpvalue <= (xblk_alpha/5)) && (xpvalue > (xblk_alpha/50))
                        xsig = '**';
                    elseif xpvalue <= (xblk_alpha/50)%0.001
                        xsig = '***';
                    end
                    
                    if xpvalue <= (0.1/n_stat_blks)
                        text(xtr + (n_tr_blks/2)+0.3, xheight, xsig, 'FontSize', 20);
                        xx = xtr + 0.7 + (0.4 * n);
                        
                        if length(xsel_conds) < 3
                            h = text(xx, y_lim(end) - 0.3, sprintf('T=%4.4f, P=%4.4f', ...
                                xstats.tstat,  xpvalue), 'FontSize', 8);
                            set(h,'Rotation',90);
                            
                            if xstats.tstat>0, ss = 1; else, ss = 2; end
                            tt{ss} = horzcat(tt{ss}, xstats.tstat);
                            pp{ss} = horzcat(pp{ss}, xpvalue);
                        end
                    end
                end
                n = n + 1;
                
                if length(xsel_conds) < 3
                    if ~isempty(tt{1})
                        text(1, 0.2 - (n * 0.05), ...
                            sprintf('Ts(%s)>%4.4f, Ps<%4.4f', num2str(xstats.df), min(tt{1}), max(pp{1})), 'FontSize', 10);
                    end
                    if  ~isempty(tt{2})
                        text(5, 0.2 - (n * 0.05), ...
                            sprintf('Ts(%s)<%4.4f, Ps<%4.4f', num2str(xstats.df), max(tt{2}), max(pp{2})), 'FontSize', 10);
                    end
                end
            end
        end
        
        %*************** anova
        if (xsel == 1) || (xsel == 2)
            xheight = y_lim(1) + 0.15;
            text(-0.5, xheight, sprintf('ANOVA %d', length(xsel_conds)));
            
            for xblk = 1:(dur_sync/n_tr_blks)
                
                xpvalue = timecourse_bcorr.targets.sync.anova{xsel}.blk{xblk}.pvalue;
                xstats  = timecourse_bcorr.targets.sync.anova{xsel}.blk{xblk}.stats;
                xtr     = (xblk-1) * n_tr_blks;
                
                if (xpvalue <= (0.1/n_stat_blks)) && (xpvalue > (xblk_alpha))
                    xsig = '+';
                elseif (xpvalue <= xblk_alpha) && (xpvalue > (xblk_alpha/5))
                    xsig = '*';
                elseif (xpvalue <= (xblk_alpha/5)) && (xpvalue > (xblk_alpha/50))
                    xsig = '**';
                elseif xpvalue <= (xblk_alpha/50)%0.001
                    xsig = '***';
                end
                
                if xpvalue <= (0.1/n_stat_blks)
                    text(xtr+(n_tr_blks/2)+0.3, xheight, xsig, 'FontSize', 20);
                    h = text(xtr+(n_tr_blks/2)+0.3, -0.1, xstats, 'FontSize', 7);
                    set(h,'Rotation',90);
                end
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
                h.FaceAlpha = .2;
                h.EdgeAlpha = .2;
                h.EdgeColor = xcond_color{xcond};
            end
        end
        
        %*************** baseline 0
        plot(x_lim_sync, [0 0], ':', 'Color', 'r', 'LineWidth', 1.5)%baseline

        %*************** real onset lines
        plot([dur_stim dur_stim]+1, y_lim, '-','Color', xonset_color)
        plot([(dur_stim + dur_manip) (dur_stim + dur_manip)]+1, y_lim, '-', 'Color', xonset_color)
        
        set(gca,'xlim', x_lim_sync, 'ylim', y_lim);
        set(gca,'XTick', x_tick_sync, 'YTick', y_lim(1):0.1:y_lim(end))
        set(gca,'XTickLabel', x_ticklable)
        
        title(sprintf('%s classifier: targets (N=%s, p<%1.4f, TR/blk=%s)', ...
            args.level,num2str(length(xsub_groups)), xblk_alpha, num2str(n_tr_blks)));
        xlabel('Volume (tr)');
        ylabel('classifier evidence (target)');
        
        %*************** real onset
        h1 = text(1.5, y_lim(1)+0.01, '(1) stim onset', 'Color', xonset_color);
        h2 = text(dur_stim + 1.5, y_lim(1)+0.01, '(2) operation onset', 'Color', xonset_color);
        h3 = text(dur_stim + dur_manip + 1.5, y_lim(1)+0.01, '(3) fixation onset', 'Color', xonset_color);
        set(h1,'Rotation', 90); set(h2,'Rotation', 90); set(h3,'Rotation', 90);
        
        %% *************** save fig
        if xfill
            fig_fname = fullfile(dirs.mvpa.group.out{xph}, ...
                sprintf('grp_plot_wn_sync_targets_sel%d_blk%d_%s', xsel, n_tr_blks, base_name));
        else
            fig_fname = fullfile(dirs.mvpa.group.out{xph}, ...
                sprintf('grp_plot_wn_sync_targets_line_sel%d_blk%d_%s', xsel, n_tr_blks, base_name));
        end
        
        savefig(xfig, sprintf('%s.fig', fig_fname));
        set(gcf,'PaperUnits','inches','PaperPosition', fig_rect_sync/100)
        print('-dpng', '-r300', sprintf('%s.png',fig_fname))
        saveas(xfig, sprintf('%s.png',fig_fname), 'png')
        
        close(xfig);
        
    end
end

end