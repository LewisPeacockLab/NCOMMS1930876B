function[grp_rsa] = grp_analysis_timecourse_rsa(grp_rsa, args, dirs)
%*************** RSA timecourse: test on study

%% ============= SETUP DIRECTORY
output_dir      = dirs.rsa.group.parse;
xsubj_grp       = args.filtered_subs;
n_subjs         = length(xsubj_grp);

%% ============= SETUP PARAMETERS
n_runs{1}       = 5;
n_runs{2}       = 6;
        
xparam          = args.index{1}.param;
n_category      = xparam.n_category;
n_subcategory   = xparam.n_subcategory;
n_targs         = 3;%1_targ, 2_nontarg, 3_baseline
n_item          = xparam.n_item;

category_names  = {'face','fruit','scene'};
subcate_names   = xparam.subcategory_name;
conds_names     = {'maintain','replace (category)','replace (subcategory)','suppress','clear'};
n_condition     = length(conds_names);
n_trs           = args.tc_tr_disp;
cate_members    = 1:n_category;
targ_names      = {'Target','RelatedNontarget','Nontarget'};
xmtcmp          = 'tukey-kramer';

xmatrix = args.index{1}.matrix;
xheader = args.index{1}.header;

for xcate = 1:n_category
    for xsubcate = 1:n_subcategory
        for xitem = 1:n_item
            it_runs{xcate}{xsubcate}{xitem} = ...
                unique(getDATA(xmatrix', xheader, ...
                {'category','subcategory','item','stimulus'}, ...
                {xcate, xsubcate, xitem, 1}, findCol(xheader,{'run'})));
            
            nn_runs{xcate}{xsubcate}(xitem) = length(it_runs{xcate}{xsubcate}{xitem});
            
            %*************** subcategory name
            args.item_names{xcate}{xsubcate}{xitem} = ...
                xparam.item_name{xcate}{xsubcate}{xitem}(1:end-4); %#ok<*NASGU>
        end
    end
end

spm_p_thresh    = args.spm_p_thresh;
spm_v_extent    = args.spm_v_extent;
spm_thresh_type = args.spm_thresh_type;

rsa_mask        = args.rsa_mask;

basename = sprintf('rsa_%s_spmT_%s_thresh%s_ext%s_%s_mask', ...
        args.level, spm_thresh_type, num2str(spm_p_thresh), num2str(spm_v_extent), rsa_mask);
    
%% *************** plot parameter
xparam = args.index{2}.param;

for xcond = 1:n_condition
    xtarg_color{1} = [238, 20, 91]/255;%#ok<*AGROW> % targ
    xtarg_color{2} = [160, 46, 83]/255;% nontarg
    xtarg_color{3} = [144, 144, 144]/255;% baseline
end
xcolor{1}   = [232 14 138]/255;%same
xcolor{2}   = [101 47 142]/255;%differ
xcolor{3}   = [0 166 156]/255;%related

xcond_color  = args.cond_color;
xcate_color  = args.cate_color;
xonset_color = args.onset_color;

dur_stim     = xparam.dur_stim;
dur_manip    = xparam.dur_manipulation;
dur_fix      = max(xparam.dur_iti);
dur_sync     = args.dur_sync;

n_tr_blks    = args.n_tr_blks;%trs in one stats block: 5=2.3s/3=1.38s
n_stat_blks  = n_trs/n_tr_blks;
xalpha       = 0.05/n_stat_blks;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%============= N Voxels/beta
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xph = 1;
for it = 1:n_subjs
    xsub  = xsubj_grp(it);
    
    for xcate = 1:n_category
        clear xsubj
        xsubj = grp_rsa{xph}.subj{xsub}.item_beta{xcate};
        
        clear t_mean t_nvox
        for xsubcate = 1:n_subcategory
            clear tt_mean tt_nvox
            for xitem = 1:n_item
                xbeta = get_mat(xsubj, 'pattern', ...
                    sprintf('spm_beta_%s', args.item_names{xcate}{xsubcate}{xitem}));
                
                tt_mean(xitem) = mean(xbeta); %#ok<*AGROW>
                tt_nvox(xitem) = length(xbeta);
            end
            
            t_mean(xsubcate) = mean(tt_mean);
            t_nvox(xsubcate) = mean(tt_nvox);
        end
        
        beta_mean(it, xcate) = mean(t_mean);
        nvox_mean(it, xcate) = mean(t_nvox);
    end
end

xtable = array2table(beta_mean, 'VariableNames', {'face','fruit','scene'});
writetable(xtable, fullfile(dirs.rsa.group.parse, sprintf('ROI_betaweight_%s.csv', basename)));

%% ============= figure
fig_rect   = [0 0 800 400];

xfig = figure;
set(xfig, 'Position', fig_rect)

for i = 1:2
    clear xmean xse b
    
    subplot(1, 2, i);
    
    if i==1
        it_mean = beta_mean;
        y_lim   = [0 1];
    else
        it_mean = nvox_mean;
        y_lim   = [0 350];
    end
    
    xmean = mean(it_mean);
    xse   = std(it_mean)/sqrt(n_subjs);
    
    for xcate = 1:n_category
        b{xcate} = bar(xcate, xmean(xcate)); hold on
        set(b{xcate}, 'facecolor', args.cate_color{xcate});
        
        %*************** mean/se
        text(xcate-0.1, xmean(xcate), ...
            sprintf(' M = %4.4f\n SE = %4.4f', xmean(xcate), xse(xcate)), 'FontSize', 10);
    end
    
    %*************** legend
    lg          = legend(category_names);
    lg.Location = 'NorthWest';%'SouthEast';
    lg.FontSize = 8;
    legend(category_names,'AutoUpdate','off')
    
    errorbar(1:n_category, xmean, xse,'k.')
    
    %*************** set
    set(gca,'xlim', [0 4],'ylim', y_lim);
    set(gca,'XTick', 1:n_category, 'XTickLabel', category_names, 'FontSize', 10);
    
    if i==1
        ylabel('item beta weight');
        title(sprintf('%s: item beta weight in RSA', args.phase_name{xph}));
    else
        ylabel('n voxels');
        title('item RSA selected voxels');
    end
    
end

%*************** save fig
fig_fname = fullfile(dirs.rsa.group.parse, sprintf('ROI_beta_nvoxs_%s', basename));

savefig(xfig, sprintf('%s.fig', fig_fname));
set(gcf,'PaperUnits','inches','PaperPosition', fig_rect/100)
print('-dpng', '-r300', sprintf('%s.png',fig_fname))
saveas(xfig, sprintf('%s.png',fig_fname), 'png')

close(xfig);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%============= TTEST
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ************* LOCALIZER within/between items
% grp_rsa{xph}.item.cate{xcate}.subcate{xsubcate}.ttest.mean(xsub, xcol)
% xcol: 1_within, 2_between
clear xcat_mean

xph = 1;

xfile    = fullfile(output_dir, sprintf('RSA_item_stats_%s_n%s.txt', basename, num2str(n_subjs)));
xout_txt = fopen(xfile, 'w+');

fprintf(xout_txt,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
fprintf(xout_txt,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
fprintf(xout_txt,'* LOCALIZER RSA\n');
fprintf(xout_txt,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
fprintf(xout_txt,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n\n');

fprintf(xout_txt,'======================================================\n');
fprintf(xout_txt,'* ITEM RSA: within vs. between items\n');
fprintf(xout_txt,'======================================================\n\n');

for i = 1:2, xcat_mean{i} = []; end
%************** paired ttest
for xcate = 1:n_category
    for xsubcate = 1:n_subcategory
        clear xmean xp xstats
        
        xmean = grp_rsa{xph}.item.ttest.cate{xcate}.subcate{xsubcate}.mean;
        for i = 1:2, xcat_mean{i} = horzcat(xcat_mean{i}, xmean(:,i)); end
        
        [~, xp, ~, xstats] = ttest(xmean(:,1), xmean(:,2));
        
        fprintf(xout_txt,'*************************************\n');
        fprintf(xout_txt,'*** category: %s, subcategory: %s\n', ...
            category_names{xcate}, subcate_names{xcate}{xsubcate});
        
        xout_stats = sprintf('T(%s)=%4.4f, P=%4.4f', num2str(xstats.df), xstats.tstat, xp);
        
        fprintf(xout_txt,'%s\n\n', xout_stats);
        fprintf(xout_txt,'within vs. between:\n mean: %4.4f/ %4.4f\n se: %4.4f/ %4.4f\n\n', ...
            mean(xmean), std(xmean)/sqrt(n_subjs));
        
        grp_rsa{xph}.item.ttest.cate{xcate}.subcate{xsubcate}.ttest.xp     = xp;
        grp_rsa{xph}.item.ttest.cate{xcate}.subcate{xsubcate}.ttest.xstats = xout_stats;
        
    end
end

%************** averaged mean/se
clear xmean xp xstats
fprintf(xout_txt,'*************************************\n');
fprintf(xout_txt,'*************************************\n');
fprintf(xout_txt,'*** average across 9 subcategories\n');

for i=1:2, xmean(:,i) = mean(xcat_mean{i},2); end

[~, xp, ~, xstats] = ttest(xmean(:,1), xmean(:,2));

xout_stats = sprintf('T(%s)=%4.4f, P=%4.4f', num2str(xstats.df), xstats.tstat, xp);

fprintf(xout_txt,'%s\n\n', xout_stats);
fprintf(xout_txt,'within vs. between:\n mean: %4.4f/ %4.4f\n se: %4.4f/ %4.4f\n\n', ...
    mean(xmean), std(xmean)/sqrt(n_subjs));

%% ============= ITEM ACROSS CATEGORY: targ vs. non-targ
% xtarg:  1_target, 2_related_nontarg 2_nontarg

fprintf(xout_txt,'======================================================\n');
fprintf(xout_txt,'* ITEM RSA: targ vs. related_nontarg vs. nontarg items\n');
fprintf(xout_txt,'======================================================\n\n');

for xtarg = 1:n_targs
    tcat_corr = [];
    
    for xcate = 1:n_category
        tt_corr = [];
        for xsubcate = 1:n_subcategory
            tt_corr = vertcat(tt_corr, grp_rsa{xph}.item.cate{xcate}.subcate{xsubcate}.targ{xtarg});
        end
        
        xrsa_corr = real(mean(tt_corr));
        grp_rsa{xph}.item_rsa.cate{xcate}.corr(:,xtarg) = xrsa_corr;
        
        tcat_corr = vertcat(tcat_corr, xrsa_corr);
    end
    
    xcat_corr = mean(tcat_corr);
    %*************** concatenated: grouped by targets
    grp_rsa{xph}.item_rsa.target.corr(:,xtarg) = xcat_corr;
end

%*************** ANOVA
for xcate = 1:n_category
    
    fprintf(xout_txt,'*************************************\n');
    fprintf(xout_txt,'*** category: %s\n\n', category_names{xcate});
    
    for xtarg = 1:n_targs
        xrsa_corr = grp_rsa{xph}.item_rsa.cate{xcate}.corr(:, xtarg);
        fprintf(xout_txt,'%s: mean: %4.4f, se: %4.4f\n', ...
            targ_names{xtarg}, mean(xrsa_corr), std(xrsa_corr)/sqrt(n_subjs));
    end

    %*************** ANOVA
    % xtarg:  1_target, 2_related_nontarg 2_nontarg
    
    clear xtable_rsa it_matrix xpvalue xpvalue_table
    
    xtable_rsa = array2table(grp_rsa{xph}.item_rsa.cate{xcate}.corr, 'VariableNames', targ_names);
    
    xanova     = {'targets'};
    xmeasures  = table((1:n_targs)','VariableNames', xanova);
    xrepmeas   = fitrm(xtable_rsa, sprintf('%s-%s~1', targ_names{1}, targ_names{end}),'WithinDesign',xmeasures);
    xanova_out = ranova(xrepmeas);
    
    xdiff_anova = sprintf('one-way ANOVA: F(%s, %s)=%2.4f, p=%1.4f', ...
        num2str(xanova_out.DF(1)), num2str(xanova_out.DF(2)),...
        xanova_out.F(1), xanova_out.pValue(1));
    
    fprintf(xout_txt,'\n*** %s\n\n', xdiff_anova);
    
    %*************** multiple comparison
    xpvalue_table = multcompare(xrepmeas, xanova,'ComparisonType', xmtcmp);%'bonferroni'
    xpvalue       = nan(n_targs-1, n_targs-1);
    
    for xcol = 1:(n_targs-1)
        for xrow = xcol+1:n_targs
            xunit = (xrow-1) + ((n_targs-1) * (xcol-1));
            xpvalue(xrow, xcol) = xpvalue_table.pValue(xunit);
        end
    end
    
    grp_rsa{xph}.item_rsa.cate{xcate}.anova.txout      = xdiff_anova;
    grp_rsa{xph}.item_rsa.cate{xcate}.anova.table      = xanova_out;
    grp_rsa{xph}.item_rsa.cate{xcate}.anova.multcomp   = xpvalue_table;
    grp_rsa{xph}.item_rsa.cate{xcate}.anova.multcomp_p = xpvalue;
    
    %*************** ttest
    clear xrsa_corr
    xrsa_corr = grp_rsa{xph}.item_rsa.cate{xcate}.corr;
    
    for xcol = 1:(n_targs-1)
        for xrow = xcol+1:n_targs
            
            xcol_corr = xrsa_corr(:, xcol);
            xrow_corr = xrsa_corr(:, xrow);
            
            [~, xp, ~, xstats] = ttest(xcol_corr, xrow_corr);
            
            fprintf(xout_txt,'****** %s vs. %s\n', targ_names{xcol}, targ_names{xrow});
            fprintf(xout_txt,'mean: %4.4f vs. %4.4f\n se: %4.4f vs. %4.4f\n\n', ...
                mean(xcol_corr), mean(xrow_corr), ...
                std(xcol_corr)/sqrt(n_subjs), std(xrow_corr)/sqrt(n_subjs));
            
            xout_stats = sprintf('T(%s)=%4.4f, P=%4.4f', num2str(xstats.df), xstats.tstat, xp);
            
            fprintf(xout_txt,'%s\n\n', xout_stats);
            
            grp_rsa{xph}.item_rsa.ttest.cate{xcate}.ttest.xp(xrow,xcol)     = xp;
            grp_rsa{xph}.item_rsa.ttest.cate{xcate}.ttest.xstats{xrow,xcol} = xout_stats;
        end
    end
end

%*************** concat
fprintf(xout_txt,'*************************************\n');
fprintf(xout_txt,'*** concatenated for targets\n');
fprintf(xout_txt,'*************************************\n\n');

for xtarg = 1:n_targs
    xrsa_corr = grp_rsa{xph}.item_rsa.target.corr(:, xtarg);
    fprintf(xout_txt,'%s: mean: %4.4f, se: %4.4f\n', ...
        targ_names{xtarg}, mean(xrsa_corr), std(xrsa_corr)/sqrt(n_subjs));
end
        
%*************** ANOVA
% xtarg:  1_target, 2_related_nontarg 2_nontarg

clear xtable_rsa it_matrix xpvalue xpvalue_table

xtable_rsa = array2table(grp_rsa{xph}.item_rsa.target.corr, 'VariableNames', targ_names);

xanova     = {'targets'};
xmeasures  = table((1:n_targs)','VariableNames', xanova);
xrepmeas   = fitrm(xtable_rsa, sprintf('%s-%s~1', targ_names{1}, targ_names{end}),'WithinDesign',xmeasures);
xanova_out = ranova(xrepmeas);

xdiff_anova = sprintf('one-way ANOVA: F(%s, %s)=%2.4f, p=%1.4f', ...
    num2str(xanova_out.DF(1)), num2str(xanova_out.DF(2)),...
    xanova_out.F(1), xanova_out.pValue(1));

fprintf(xout_txt,'\n*** %s\n\n', xdiff_anova);

%*************** multiple comparison
xpvalue_table = multcompare(xrepmeas, xanova,'ComparisonType', xmtcmp);%'bonferroni'
xpvalue       = nan(n_targs-1, n_targs-1);

for xcol = 1:(n_targs-1)
    for xrow = xcol+1:n_targs
        xunit = (xrow-1) + ((n_targs-1) * (xcol-1));
        xpvalue(xrow, xcol) = xpvalue_table.pValue(xunit);
    end
end

grp_rsa{xph}.item_rsa.target.anova.txout      = xdiff_anova;
grp_rsa{xph}.item_rsa.target.anova.table      = xanova_out;
grp_rsa{xph}.item_rsa.target.anova.multcomp   = xpvalue_table;
grp_rsa{xph}.item_rsa.target.anova.multcomp_p = xpvalue;

%*************** ttest
for xcol = 1:(n_targs-1)
    for xrow = xcol+1:n_targs
        
        xcol_corr = grp_rsa{xph}.item_rsa.target.corr(:, xcol);
        xrow_corr = grp_rsa{xph}.item_rsa.target.corr(:, xrow);
        
        [~, xp, ~, xstats] = ttest(xcol_corr, xrow_corr);
        
        fprintf(xout_txt,'****** %s vs. %s\n', targ_names{xcol}, targ_names{xrow});
        fprintf(xout_txt,'mean: %4.4f vs. %4.4f\n se: %4.4f vs. %4.4f\n\n', ...
            mean(xcol_corr), mean(xrow_corr), ...
            std(xcol_corr)/sqrt(n_subjs), std(xrow_corr)/sqrt(n_subjs));
        
        xout_stats = sprintf('T(%s)=%4.4f, P=%4.4f', num2str(xstats.df), xstats.tstat, xp);
        
        fprintf(xout_txt,'%s\n\n', xout_stats);
        
        grp_rsa{xph}.item_rsa.target.ttest.xp(xrow,xcol)     = xp;
        grp_rsa{xph}.item_rsa.target.ttest.xstats{xrow,xcol} = xout_stats;
    end
end

%% ============= STUDY PERCEPTION
% args.percept_win
% xtarg:  1_target, 2_related_nontarg 2_nontarg
% xlevel: 1_category, 2_subcategory
%*************** mean: target/ non-target
xpercept_win = args.percept_win;%
xph      = 2;
it_conds = [1 2 4 5];

for xtarg = 1:n_targs
    clear xrsa_corr
    tcat_corr = [];
    
    for xcate = 1:n_category
        t_corr = [];
        for xsubcate = 1:n_subcategory
            tt_corr = [];
            for xcond = it_conds
                ttt_corr = [];
                for xtr = xpercept_win
                    ttt_corr = vertcat(ttt_corr, grp_rsa{xph}.timecourse.cond{xcond}.cate{xcate}.subcate{xsubcate}.targ{xtarg}.tr{xtr}.corr_w);
                end
                
                tt_corr = vertcat(tt_corr, mean(ttt_corr));
            end
            
            t_corr = vertcat(t_corr, mean(tt_corr));
        end
        
        xrsa_corr = mean(t_corr);
        tcat_corr = vertcat(tcat_corr, xrsa_corr);
         
        grp_rsa{xph}.percept.cate{xcate}.corr(:,xtarg) = xrsa_corr;
    end
    
    xcat_corr = mean(tcat_corr);
    grp_rsa{xph}.percept.corr(:,xtarg) = xcat_corr;
end

fprintf(xout_txt,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
fprintf(xout_txt,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
fprintf(xout_txt,'* STUDY RSA window: %d\n', length(xpercept_win));
fprintf(xout_txt,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
fprintf(xout_txt,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n\n');

fprintf(xout_txt,'======================================================\n');
fprintf(xout_txt,'* ITEM RSA: targ vs. related_nontarg vs. nontarg items\n');
fprintf(xout_txt,'======================================================\n\n');

%*************** ANOVA
for xcate = 1:n_category
    
    fprintf(xout_txt,'*************************************\n');
    fprintf(xout_txt,'*** category: %s\n\n', category_names{xcate});
    
    for xtarg = 1:n_targs
        xrsa_corr = grp_rsa{xph}.percept.cate{xcate}.corr(:, xtarg);
        fprintf(xout_txt,'%s: mean: %4.4f, se: %4.4f\n', ...
            targ_names{xtarg}, mean(xrsa_corr), std(xrsa_corr)/sqrt(n_subjs));
    end

    %*************** ANOVA
    % xtarg:  1_target, 2_related_nontarg 2_nontarg
    
    clear xtable_rsa it_matrix xpvalue xpvalue_table
    xtable_rsa = array2table(real(grp_rsa{xph}.percept.cate{xcate}.corr), 'VariableNames', targ_names);
    
    xanova     = {'targets'};
    xmeasures  = table((1:n_targs)','VariableNames', xanova);
    xrepmeas   = fitrm(xtable_rsa, sprintf('%s-%s~1', targ_names{1}, targ_names{end}),'WithinDesign',xmeasures);
    xanova_out = ranova(xrepmeas);
    
    xdiff_anova = sprintf('one-way ANOVA: F(%s, %s)=%2.4f, p=%1.4f', ...
        num2str(xanova_out.DF(1)), num2str(xanova_out.DF(2)),...
        xanova_out.F(1), xanova_out.pValue(1));
    
    fprintf(xout_txt,'\n*** %s\n\n', xdiff_anova);
    
    %*************** multiple comparison
    xpvalue_table = multcompare(xrepmeas, xanova,'ComparisonType', xmtcmp);%'bonferroni'
    xpvalue       = nan(n_targs-1, n_targs-1);
    
    for xcol = 1:(n_targs-1)
        for xrow = xcol+1:n_targs
            xunit = (xrow-1) + ((n_targs-1) * (xcol-1));
            xpvalue(xrow, xcol) = xpvalue_table.pValue(xunit);
        end
    end
    
    grp_rsa{xph}.percept.cate{xcate}.anova.txout      = xdiff_anova;
    grp_rsa{xph}.percept.cate{xcate}.anova.table      = xanova_out;
    grp_rsa{xph}.percept.cate{xcate}.anova.multcomp   = xpvalue_table;
    grp_rsa{xph}.percept.cate{xcate}.anova.multcomp_p = xpvalue;
    
    %*************** ttest
    clear xrsa_corr
    xrsa_corr = grp_rsa{xph}.percept.cate{xcate}.corr;
    
    for xcol = 1:(n_targs-1)
        for xrow = xcol+1:n_targs
            
            xcol_corr = xrsa_corr(:, xcol);
            xrow_corr = xrsa_corr(:, xrow);
            
            [~, xp, ~, xstats] = ttest(xcol_corr, xrow_corr);
            
            fprintf(xout_txt,'****** %s vs. %s\n', targ_names{xcol}, targ_names{xrow});
            fprintf(xout_txt,'mean: %4.4f vs. %4.4f\n se: %4.4f vs. %4.4f\n\n', ...
                mean(xcol_corr), mean(xrow_corr), ...
                std(xcol_corr)/sqrt(n_subjs), std(xrow_corr)/sqrt(n_subjs));
            
            xout_stats = sprintf('T(%s)=%4.4f, P=%4.4f', num2str(xstats.df), xstats.tstat, xp);
            
            fprintf(xout_txt,'%s\n\n', xout_stats);
            
            grp_rsa{xph}.percept.cate{xcate}.ttest.xp(xrow,xcol)     = xp;
            grp_rsa{xph}.percept.cate{xcate}.ttest.xstats{xrow,xcol} = xout_stats;
        end
    end
end

%*************** concat
fprintf(xout_txt,'*************************************\n');
fprintf(xout_txt,'*** concatenated for targets\n');
fprintf(xout_txt,'*************************************\n\n');

for xtarg = 1:n_targs
    xrsa_corr = grp_rsa{xph}.percept.corr(:, xtarg);
    fprintf(xout_txt,'%s: mean: %4.4f, se: %4.4f\n', ...
        targ_names{xtarg}, mean(xrsa_corr), std(xrsa_corr)/sqrt(n_subjs));
end
        
%*************** ANOVA
% xtarg:  1_target, 2_related_nontarg 2_nontarg

clear xtable_rsa it_matrix xpvalue xpvalue_table

xtable_rsa = array2table(real(grp_rsa{xph}.percept.corr), 'VariableNames', targ_names);

xanova     = {'targets'};
xmeasures  = table((1:n_targs)','VariableNames', xanova);
xrepmeas   = fitrm(xtable_rsa, sprintf('%s-%s~1', targ_names{1}, targ_names{end}),'WithinDesign',xmeasures);
xanova_out = ranova(xrepmeas);

xdiff_anova = sprintf('one-way ANOVA: F(%s, %s)=%2.4f, p=%1.4f', ...
    num2str(xanova_out.DF(1)), num2str(xanova_out.DF(2)),...
    xanova_out.F(1), xanova_out.pValue(1));

fprintf(xout_txt,'\n*** %s\n\n', xdiff_anova);

%*************** multiple comparison
xpvalue_table = multcompare(xrepmeas, xanova,'ComparisonType', xmtcmp);%'bonferroni'
xpvalue       = nan(n_targs-1, n_targs-1);

for xcol = 1:(n_targs-1)
        for xrow = xcol+1:n_targs
            xunit = (xrow-1) + ((n_targs-1) * (xcol-1));
            xpvalue(xrow, xcol) = xpvalue_table.pValue(xunit);
        end
end
    
grp_rsa{xph}.percept.anova.txout      = xdiff_anova;
grp_rsa{xph}.percept.anova.table      = xanova_out;
grp_rsa{xph}.percept.anova.multcomp   = xpvalue_table;
grp_rsa{xph}.percept.anova.multcomp_p = xpvalue;

%*************** ttest
for xcol = 1:(n_targs-1)
    for xrow = xcol+1:n_targs
        
        xcol_corr = grp_rsa{xph}.percept.corr(:, xcol);
        xrow_corr = grp_rsa{xph}.percept.corr(:, xrow);
        
        [~, xp, ~, xstats] = ttest(xcol_corr, xrow_corr);
        
        fprintf(xout_txt,'****** %s vs. %s\n', targ_names{xcol}, targ_names{xrow});
        fprintf(xout_txt,'mean: %4.4f vs. %4.4f\n se: %4.4f vs. %4.4f\n\n', ...
            mean(xcol_corr), mean(xrow_corr), ...
            std(xcol_corr)/sqrt(n_subjs), std(xrow_corr)/sqrt(n_subjs));
        
        xout_stats = sprintf('T(%s)=%4.4f, P=%4.4f', num2str(xstats.df), xstats.tstat, xp);
        
        fprintf(xout_txt,'%s\n\n', xout_stats);
        
        grp_rsa{xph}.percept.ttest.xp(xrow,xcol)     = xp;
        grp_rsa{xph}.percept.ttest.xstats{xrow,xcol} = xout_stats;
    end
end

%% ============= LOCALIZER VS. STUDY PERCEPTION
fprintf(xout_txt,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
fprintf(xout_txt,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
fprintf(xout_txt,'* LOCALIZER VS. STUDY RSA\n');
fprintf(xout_txt,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
fprintf(xout_txt,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n\n');

fprintf(xout_txt,'======================================================\n');
fprintf(xout_txt,'* 2-way ANOVA: phase * target\n');
fprintf(xout_txt,'======================================================\n\n');

clear xfactor xarray
xfactorial_vars = {'measure','subject','phase','target'};
xfactor         = {'phase','target'};

xarray = zeros(2 * n_targs * n_subjs, length(xfactorial_vars));

for xph = 1:2
    for xtarg = 1:n_targs
        clear xrsa_corr
        xunit = (1:n_subjs) + (n_subjs * (xtarg-1)) + (n_subjs * n_targs * (xph-1));
        
        if xph==1
            xrsa_corr = grp_rsa{xph}.item_rsa.target.corr(:, xtarg);
        else
            xrsa_corr = grp_rsa{xph}.percept.corr(:, xtarg);
        end
        
        xarray(xunit, findCol(xfactorial_vars, {'measure'})) = xrsa_corr;
        xarray(xunit, findCol(xfactorial_vars, {'subject'})) = (1:n_subjs);
        xarray(xunit, findCol(xfactorial_vars, {'phase'}))   = xph;
        xarray(xunit, findCol(xfactorial_vars, {'target'}))  = xtarg;
    end
end

Y  = xarray(:, findCol(xfactorial_vars, {'measure'}));
S  = xarray(:, findCol(xfactorial_vars, {'subject'}));
F1 = xarray(:, findCol(xfactorial_vars, {'phase'}));
F2 = xarray(:, findCol(xfactorial_vars, {'target'}));

xstats = rm_anova2(Y, S, F1, F2, xfactor);
xtable = cell2table(xstats);

grp_rsa{xph}.percept.two_way.table = xtable;

for xmain = 1:length(xfactor)
    xunit     = xmain + 1;
    xtx_stats = sprintf('%s: F(%s,%s)=%4.4f, p=%4.4f', ...
        xfactor{xmain}, num2str(xtable.xstats3{xunit}), ...
        num2str(xtable.xstats3{length(xtable.xstats1)}), ...
        xtable.xstats5{xunit}, xtable.xstats6{xunit});
    
    fprintf(xout_txt,'Main effect: %s\n', xtx_stats);
    
    grp_rsa{xph}.percept.two_way.CondxPair.stats.main{xmain}  = xtx_stats;
    grp_rsa{xph}.percept.two_way.CondxPair.pvalue.main{xmain} = xtable.xstats6{xunit};
end

xunit     = 4;
xtx_stats = sprintf('%s x %s: F(%s,%s)=%4.4f, p=%4.4f', ...
    xfactor{1}, xfactor{2}, num2str(xtable.xstats3{xunit}), ...
    num2str(xtable.xstats3{length(xtable.xstats1)}), ...
    xtable.xstats5{xunit}, xtable.xstats6{xunit});

fprintf(xout_txt,'Interaction: %s\n', xtx_stats);

grp_rsa{xph}.percept.two_way.CondxPair.stats.interact  = xtx_stats;
grp_rsa{xph}.percept.two_way.CondxPair.pvalue.interact = xtable.xstats6{xunit};

%% *************** figure
clear xrsa_corr xpvalue xanova
xrsa_corr = [];
for xcate = 1:n_category
    xrsa_corr{1}{xcate}   = grp_rsa{1}.item_rsa.cate{xcate}.corr;
    xpvalue{1}{xcate} = grp_rsa{1}.item_rsa.cate{xcate}.anova.multcomp_p;
    xanova{1}{xcate}  = grp_rsa{1}.item_rsa.cate{xcate}.anova.txout;
    
    xrsa_corr{2}{xcate}   = grp_rsa{2}.percept.cate{xcate}.corr;
    xpvalue{2}{xcate} = grp_rsa{2}.percept.cate{xcate}.anova.multcomp_p;
    xanova{2}{xcate}  = grp_rsa{2}.percept.cate{xcate}.anova.txout;
end

xrsa_corr{1}{4}   = grp_rsa{1}.item_rsa.target.corr;
xpvalue{1}{4} = grp_rsa{1}.item_rsa.target.anova.multcomp_p;
xanova{1}{4}  = grp_rsa{1}.item_rsa.target.anova.txout;

xrsa_corr{2}{4} = grp_rsa{2}.percept.corr;
xpvalue{2}{4} = grp_rsa{2}.percept.anova.multcomp_p;
xanova{2}{4}  = grp_rsa{2}.percept.anova.txout;

%*************** bar plot
xtkalpha   = 0.05;
if strcmp(args.mask_name, 'hippocampus')
    y_lim  = [-0.05 0.45];
elseif strcmp(args.mask_name, 'PPAsphere_10_bin')
    y_lim  = [-0.05 0.55];
else
    y_lim  = [0 0.45];
end
xfont_size = 8;
fig_rect   = [0 0 2000 1000];

xfig = figure;
set(xfig, 'Position', fig_rect)

for xph = 1:2
    for xcate = 1:(n_category+1)
        clear xmean xse b
        
        subplot(2, 4, xcate + (4 * (xph-1)));
        
        xcorrs = real(xrsa_corr{xph}{xcate});
        xcorrs(find(isnan(xcorrs)))=0;
        xmean = mean(xcorrs);
        xse   = std(xcorrs)/sqrt(n_subjs);
        
        for xtarg = 1:n_targs
            b{xtarg} = bar(xtarg, xmean(xtarg)); hold on
            set(b{xtarg}, 'facecolor', xtarg_color{xtarg});
            
            %*************** mean/se
            text(xtarg-0.1, xmean(xtarg) + 0.03, ...
                sprintf(' M = %4.4f\n SE = %4.4f', xmean(xtarg), xse(xtarg)), 'FontSize', xfont_size);
        end
        
        %*************** legend
        if (xph==1) && (xcate == 1)
            lg          = legend(targ_names);
            lg.Location = 'Best';%'SouthEast';
            lg.FontSize = 8;   
            legend(targ_names,'AutoUpdate','off')
        end
        
        errorbar(1:n_targs, xmean, xse,'k.')
        
        %*************** ttest
        n = 1;
        for xcol = 1:(n_targs-1)
            for xrow = xcol+1:n_targs
                xp = xpvalue{xph}{xcate}(xrow, xcol);
                
                if xp <= 0.1
                    yy_sig = y_lim(2) - (0.02 * n);
                    n = n + 1;
                    
                    if (xp <= 0.1) && (xp > xtkalpha)
                        xsig = '+';
                    elseif (xp <= xtkalpha) && (xp > (xtkalpha/5))
                        xsig = '*';
                    elseif (xp <= (xtkalpha/5)) && (xp > (xtkalpha/50))
                        xsig = '**';
                    elseif xp <= (xtkalpha/50)%0.001
                        xsig = '***';
                    end
                    
                    txt_x = xcol + (xrow-xcol)/2;
                    plot([xcol xrow],[yy_sig yy_sig]-0.005,'k-');
                    text(txt_x, yy_sig, xsig, 'FontSize', 20);
                    text(txt_x + 0.2, yy_sig+0.001, sprintf('p=%4.3f', xp), 'FontSize', xfont_size);
                end
            end
        end
        
        %*************** one-way anova
        text(0.1, y_lim(1) + 0.02, xanova{xph}{xcate}, 'FontSize', 8);
        
        %*************** two-way anova
        if xph == 2
            text(0.1, 0.32, grp_rsa{xph}.percept.two_way.CondxPair.stats.interact, 'FontSize', xfont_size);
        end
        
        %*************** set
        set(gca,'xlim', [0 4],'ylim', y_lim);
        set(gca,'XTick', 1:n_targs, 'XTickLabel', targ_names, 'FontSize', xfont_size);
        
        if xcate~=4, xname = category_names{xcate};
        else,        xname = 'concat'; end
        
        title(sprintf('%s/ %s: item RSA decoding: %s', ...
            args.phase_name{xph}, xname, xmtcmp));
        ylabel('RSA');
    end
end

%*************** save fig
fig_fname = fullfile(output_dir, sprintf('bar_RSA_decoding_%s_n%s_win%d', ...
    basename, num2str(n_subjs), length(xpercept_win)));

savefig(xfig, sprintf('%s.fig', fig_fname));
set(gcf,'PaperUnits','inches','PaperPosition', fig_rect/100)
print('-dpng', '-r300', sprintf('%s.png',fig_fname))
saveas(xfig, sprintf('%s.png',fig_fname), 'png')

% close(xfig);

%% ************** save in csv
xarray = grp_rsa{2}.percept.corr;%{'Target','RelatedNontarget','Nontarget'}
for xtarg = 1:3
    xunit = (1:n_subjs) + (n_subjs * (xtarg-1));
    
    xprecept_rsa(xunit, 1) = 1:n_subjs;
    xprecept_rsa(xunit, 2) = xtarg;
    xprecept_rsa(xunit, 3) = xarray(:, xtarg);
end

xtable = array2table(xprecept_rsa, 'VariableNames', {'subject','target','decoding_r'});
writetable(xtable, fullfile(dirs.rsa.group.pattern, 'group_rsa_study_decoding.csv'));

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%============= LOCALIZER PLOT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%************** CONFUSION MATRIX WITHIN SUBCATEGORY
xph = 1;

for xcate = 1:n_category
    for xsubcate = 1:n_subcategory
        grp_rsa{xph}.item.cate{xcate}.subcate{xsubcate}.heatmap   = [];
        grp_rsa{xph}.item.cate{xcate}.subcate{xsubcate}.heatmap_w = [];
    end
end

for xcate = 1:n_category
    for xsubcate = 1:n_subcategory
        for xitem = 1:n_item
            for xrun = 1:nn_runs{xcate}{xsubcate}(xitem)
                
                xcell = xrun + sum(nn_runs{xcate}{xsubcate}(1:(xitem-1)));
                
                for yitem = 1:n_item
                    for yrun = 1:nn_runs{xcate}{xsubcate}(yitem)
                        
                        ycell = yrun + sum(nn_runs{xcate}{xsubcate}(1:(yitem-1)));
                        
                        grp_rsa{xph}.item.cate{xcate}.subcate{xsubcate}.heatmap(xcell, ycell) = ...
                            mean(grp_rsa{xph}.item.cate{xcate}.subcate{xsubcate}.corr{xcell, ycell});
                        
                        grp_rsa{xph}.item.cate{xcate}.subcate{xsubcate}.heatmap_w(xcell, ycell) = ...
                            mean(grp_rsa{xph}.item.cate{xcate}.subcate{xsubcate}.corr_w{xcell, ycell});
                        
                    end
                end
            end
        end
    end
end

xmin = []; xmax = []; xmin_w = []; xmax_w = [];

for xcate = 1:n_category
    for xsubcate = 1:n_subcategory
        xx     = grp_rsa{xph}.item.cate{xcate}.subcate{xsubcate}.heatmap;
        xx_w   = grp_rsa{xph}.item.cate{xcate}.subcate{xsubcate}.heatmap_w;
        
        xmin   = horzcat(xmin, min(xx(~isinf(xx))));
        xmin_w = horzcat(xmin_w, min(xx_w(~isinf(xx_w))));
        
        xmax   = horzcat(xmax, max(xx(~isinf(xx))));
        xmax_w = horzcat(xmax_w, max(xx_w(~isinf(xx_w))));
    end
end

%% --------- figure
clear xticks

for xcate = 1:n_category
    for xsubcate = 1:n_subcategory
        for xitem = 1:n_item
            for xrun = 1:nn_runs{xcate}{xsubcate}(xitem)
                
                xticks{xcate}{xsubcate}(xitem) = ...
                    (nn_runs{xcate}{xsubcate}(xitem)/2) + ...
                    sum(nn_runs{xcate}{xsubcate}(1:xitem-1));
            end
        end
    end
end

rect_w      = 1200;
fig_rect    = [0 0 rect_w (rect_w/1.7)];
bar_lim     = [min(xmin) max(xmax)];%0.9
bar_lim_w   = [min(xmin_w) max(xmax_w)];%0.9

for xcate = 1:n_category
    heatmap_fig = figure;
    set(heatmap_fig, 'Position', fig_rect)
    
    for xsubcate = 1:n_subcategory
        %*************** similarity: pearson correlation
        subplot(2, n_subcategory, xsubcate)
        imagesc(grp_rsa{xph}.item.cate{xcate}.subcate{xsubcate}.heatmap);
        colormap('jet'); colorbar; caxis(bar_lim); hold on
        
        %*************** boundary
        for xitem = 1:(n_item-1)
            plot([0, sum(nn_runs{xcate}{xsubcate})]+0.5,...
                [sum(nn_runs{xcate}{xsubcate}(1:xitem)), sum(nn_runs{xcate}{xsubcate}(1:xitem))]+0.5, ...
                'k--','LineWidth',1);
            
            plot([sum(nn_runs{xcate}{xsubcate}(1:xitem)), sum(nn_runs{xcate}{xsubcate}(1:xitem))]+0.5,...
                [0, sum(nn_runs{xcate}{xsubcate})]+0.5, ...
                'k--','LineWidth',1);
        end
        
        title(sprintf('Localizer RSA: %s-%s (N=%s)', ...
            category_names{xcate}, subcate_names{xcate}{xsubcate}, ...
            num2str(n_subjs)),'FontSize',15,'FontWeight','bold');
        xlabel('item', 'FontSize', 15);
        ylabel('item', 'FontSize', 15);
        
        set(gca,'YTick', xticks{xcate}{xsubcate}, 'YTickLabel', 1:n_item, 'FontSize', 10);
        set(gca,'XTick', xticks{xcate}{xsubcate}, 'XTickLabel', 1:n_item, 'FontSize', 10);
        set(gca,'xaxisLocation','top')
        
        %*************** similarity: pearson correlation: weighted
        subplot(2, n_subcategory, xsubcate + n_subcategory)
        imagesc(grp_rsa{xph}.item.cate{xcate}.subcate{xsubcate}.heatmap_w);
        colormap('jet'); colorbar; caxis(bar_lim_w); hold on
        
        %*************** boundary
        for xitem = 1:(n_item-1)
            plot([0, sum(nn_runs{xcate}{xsubcate})]+0.5,...
                [sum(nn_runs{xcate}{xsubcate}(1:xitem)), sum(nn_runs{xcate}{xsubcate}(1:xitem))]+0.5, ...
                'k--','LineWidth',1);
            
            plot([sum(nn_runs{xcate}{xsubcate}(1:xitem)), sum(nn_runs{xcate}{xsubcate}(1:xitem))]+0.5,...
                [0, sum(nn_runs{xcate}{xsubcate})]+0.5, ...
                'k--','LineWidth',1);
        end
        
        title(sprintf('Localizer RSA w: %s-%s (N=%s)', ...
            category_names{xcate}, subcate_names{xcate}{xsubcate}, ...
            num2str(n_subjs)),'FontSize',15,'FontWeight','bold');
        xlabel('item', 'FontSize', 15);
        ylabel('item', 'FontSize', 15);
        
        set(gca,'YTick', xticks{xcate}{xsubcate}, 'YTickLabel', 1:n_item, 'FontSize', 10);
        set(gca,'XTick', xticks{xcate}{xsubcate}, 'XTickLabel', 1:n_item, 'FontSize', 10);
        set(gca,'xaxisLocation','top')
    end
    
    %% *************** save fig
    fig_fname = fullfile(output_dir, sprintf('plot_similarity_matrix_%s_%s_%s_n%s', ...
        category_names{xcate}, args.phase_name{xph}, basename, num2str(n_subjs)));
    
    savefig(heatmap_fig, sprintf('%s.fig', fig_fname));
    set(gcf,'PaperUnits','inches','PaperPosition', fig_rect/100)
    print('-dpng', '-r300', sprintf('%s.png',fig_fname))
    saveas(heatmap_fig, sprintf('%s.png',fig_fname), 'png')
    
    close(heatmap_fig);
end

%% ============= TEMPLATE CONFUSION MATRIX WITHIN LOCALIZER: MEAN
xph = 1;

for xmaskcate = 1:n_category
    for xcate = 1:n_category
        for xsubcate = 1:n_subcategory
            for xitem = 1:n_item
                
                xcell = xitem + (n_item * (xsubcate-1)) + (n_item * n_subcategory * (xcate-1));
                
                for ycate = 1:n_category
                    for ysubcate = 1:n_subcategory
                        for yitem = 1:n_item
                            
                            ycell = yitem + (n_item * (ysubcate-1)) + (n_item * n_subcategory * (ycate-1));
                            
                            grp_rsa{xph}.template.maskcate{xmaskcate}.heatmap(xcell, ycell) = ...
                                mean(grp_rsa{xph}.template.maskcate{xmaskcate}.corr{xcell, ycell});
                            
                            grp_rsa{xph}.template.maskcate{xmaskcate}.heatmap_w(xcell, ycell) = ...
                                mean(grp_rsa{xph}.template.maskcate{xmaskcate}.corr_w{xcell, ycell});
                        end
                        
                    end
                end
            end
        end
    end
end

xmin = []; xmax = []; xmin_w = []; xmax_w = [];

for xmaskcate = 1:n_category
    
    xx     = grp_rsa{xph}.template.maskcate{xmaskcate}.heatmap;
    xx_w   = grp_rsa{xph}.template.maskcate{xmaskcate}.heatmap_w;
    
    xmin   = horzcat(xmin, min(xx(~isinf(xx))));
    xmin_w = horzcat(xmin_w, min(xx_w(~isinf(xx_w))));
    
    xmax   = horzcat(xmax, max(xx(~isinf(xx))));
    xmax_w = horzcat(xmax_w, max(xx_w(~isinf(xx_w))));
end

for xcate = 1:n_category
    for xsubcate = 1:n_subcategory
        xnames{xsubcate + n_subcategory * (xcate-1)} = ...
            subcate_names{xcate}{xsubcate};
    end
end

%% --------- figure
clear xticks

all_items   = n_item * n_category * n_subcategory;
xticks      = (n_item/2):n_item:all_items;

fig_rect    = [0 0 rect_w (rect_w/2.4)];
bar_lim     = [min(xmin) max(xmax)];%0.9
bar_lim_w   = [min(xmin_w) max(xmax_w)];%0.9

for xmaskcate = 1:n_category
    heatmap_fig = figure;
    set(heatmap_fig, 'Position', fig_rect)
    
    %*************** similarity: pearson correlation
    subplot(1, 2, 1)
    imagesc(grp_rsa{xph}.template.maskcate{xmaskcate}.heatmap);
    colormap('jet'); colorbar; caxis(bar_lim); hold on
    
    %*************** boundary
    for xcate = 1:n_category
        for xsubcate = 1:(n_subcategory-1)
            xx = n_item * xsubcate + (n_item * n_subcategory * (xcate-1));
            plot([0, all_items]+0.5, [xx xx]+0.5, 'k--','LineWidth', 1);
            plot([xx xx]+0.5, [0, all_items]+0.5, 'k--','LineWidth', 1);
        end
        
        xx = (n_item * n_subcategory * xcate);
        plot([0, all_items]+0.5, [xx, xx]+0.5, 'k-','LineWidth', 1);
        plot([xx, xx]+0.5, [0, all_items]+0.5, 'k-','LineWidth', 1);
        
    end
    
    title(sprintf('Localizer RSA: ROI: %s (N=%s)', ...
        category_names{xmaskcate}, num2str(n_subjs)),...
        'FontSize',15,'FontWeight','bold');
    xlabel('item', 'FontSize', 15);
    ylabel('item', 'FontSize', 15);
    xtickangle(45)
    
    set(gca,'YTick', xticks, 'YTickLabel', xnames, 'FontSize', 10);
    set(gca,'XTick', xticks, 'XTickLabel', xnames, 'FontSize', 10);
    set(gca,'xaxisLocation','top')
    
    %*************** similarity: pearson correlation: weighted
    subplot(1, 2, 2)
    imagesc(grp_rsa{xph}.template.maskcate{xmaskcate}.heatmap_w);
    colormap('jet'); colorbar; caxis(bar_lim_w); hold on
    
    %*************** boundary
    for xcate = 1:n_category
        for xsubcate = 1:(n_subcategory-1)
            xx = n_item * xsubcate + (n_item * n_subcategory * (xcate-1));
            plot([0, all_items]+0.5, [xx xx]+0.5, 'k--','LineWidth',2);
            plot([xx xx]+0.5, [0, all_items]+0.5, 'k--','LineWidth',2);
        end
        
        xx = (n_item * n_subcategory * xcate);
        plot([0, all_items]+0.5, [xx, xx]+0.5, 'k-','LineWidth',2);
        plot([xx, xx]+0.5, [0, all_items]+0.5, 'k-','LineWidth',2);
        
    end
    
    title(sprintf('Localizer RSA w: ROI: %s (N=%s)', ...
        category_names{xmaskcate}, num2str(n_subjs)),...
        'FontSize',15,'FontWeight','bold');
    xlabel('item', 'FontSize', 15);
    ylabel('item', 'FontSize', 15);
    xtickangle(45)
    
    set(gca,'YTick', xticks, 'YTickLabel', xnames, 'FontSize', 10);
    set(gca,'XTick', xticks, 'XTickLabel', xnames, 'FontSize', 10);
    set(gca,'xaxisLocation','top')
    
    %% *************** save fig
    fig_fname = fullfile(output_dir, sprintf('plot_template_similarity_matrix_%s_%s_%s_n%s', ...
        category_names{xmaskcate}, args.phase_name{xph}, basename, num2str(n_subjs)));
    
    savefig(heatmap_fig, sprintf('%s.fig', fig_fname));
    set(gcf,'PaperUnits','inches','PaperPosition', fig_rect/100)
    print('-dpng', '-r300', sprintf('%s.png',fig_fname))
    saveas(heatmap_fig, sprintf('%s.png',fig_fname), 'png')
    
    close(heatmap_fig);
    
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%============= PREDICTING NEW ITEM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%*************** for replace category/subcategory
%*************** repCat (2): targ:1_item, 2_new_item, 3_related_item, 4_related new_item, 5_others
%*************** repSub (3): targ:1_item, 2_new_item, 3_related, 4_others
%*************** grp_rsa{xph}.subj{xsub}.timecourse.cond{xcond}.item_targ{xtarg}.tr{xtr}.corr_w
xph       = 2;
it_conds  = 2:3;
item_peak = args.item_peak;%for replace: 1_before, 2_after switch
n_peak    = item_peak(2, 2) - item_peak(2, 1) + 1;

for xcond = it_conds
    
    if xcond==2, n_targs = 5; else, n_targs = 4; end
    
    grp_rsa{xph}.timecourse.item_targ.cond{xcond}.mean_w = zeros(n_targs, n_trs);
    grp_rsa{xph}.timecourse.item_targ.cond{xcond}.se_w   = zeros(n_targs, n_trs);
    
    for xtarg = 1:n_targs
        for xtr = 1:n_trs
            tt = [];
            for xsub = xsubj_grp
                tt = horzcat(tt, ...
                    mean(grp_rsa{xph}.subj{xsub}.timecourse.cond{xcond}.item_targ{xtarg}.tr{xtr}.corr_w));
            end
            
            grp_rsa{xph}.timecourse.item_targ.cond{xcond}.targ{xtarg}.tr{xtr}.corr_w = tt;
            grp_rsa{xph}.timecourse.item_targ.cond{xcond}.mean_w(xtarg, xtr) = mean(tt);
            grp_rsa{xph}.timecourse.item_targ.cond{xcond}.se_w(xtarg, xtr)   = std(tt)/sqrt(n_subjs);
        end
    end
end

%% *************** baseline correction
clear xbase_mean
for xcond = it_conds
    
    if xcond==2, n_targs = 5; else, n_targs = 4; end
    
    for xtarg = 1:n_targs
        
        for it = 1:n_subjs
            xbase = [];
            for xtr = (1:args.baseline_trs)
                xbase = horzcat(xbase, ...
                    grp_rsa{xph}.timecourse.item_targ.cond{xcond}.targ{xtarg}.tr{xtr}.corr_w(it));
            end
            
            xbase_mean{it}{xcond}(xtarg) = mean(xbase);
            
            for xtr = 1:n_trs
                grp_rsa{xph}.timecourse.item_targ_bs.cond{xcond}.targ{xtarg}.tr{xtr}.corr_w(it) = ...
                    grp_rsa{xph}.timecourse.item_targ.cond{xcond}.targ{xtarg}.tr{xtr}.corr_w(it) - mean(xbase);
                
            end
        end
        
        for xtr = 1:n_trs
            tt = grp_rsa{xph}.timecourse.item_targ_bs.cond{xcond}.targ{xtarg}.tr{xtr}.corr_w;
            grp_rsa{xph}.timecourse.item_targ_bs.cond{xcond}.mean_w(xtarg, xtr) = mean(tt);
            grp_rsa{xph}.timecourse.item_targ_bs.cond{xcond}.se_w(xtarg, xtr)   = std(tt)/sqrt(n_subjs);
        end
    end
end

%% *************** STATS: RANDOM
for xcond = it_conds
    
    if xcond==2, n_targs = 5; else, n_targs = 4; end
    
    for xblk = 1:n_stat_blks
        
        grp_rsa{xph}.timecourse.item_targ.cond{xcond}.blk{xblk}.pvalue_w = nan(n_targs, n_targs);
        grp_rsa{xph}.timecourse.item_targ_bs.cond{xcond}.blk{xblk}.pvalue_w = nan(n_targs, n_targs);
        
        for xcol = 1:(n_targs-1)
            for xrow = (xcol+1):n_targs
                clear xpat
                
                it_trs = (1:n_tr_blks) + (n_tr_blks * (xblk-1));
                
                %#########################################################
                %*************** similarity: pearson correlation
                xpat = grp_rsa{xph}.timecourse.item_targ.cond{xcond};
                xtarg_col = []; xtarg_row = [];
                
                for xtr = it_trs
                    xtarg_col = vertcat(xtarg_col, xpat.targ{xcol}.tr{xtr}.corr_w);
                    xtarg_row = vertcat(xtarg_row, xpat.targ{xrow}.tr{xtr}.corr_w);
                end
                
                %*************** ttest
                [~, xpvalue, ~, xstats] = ttest(mean(xtarg_col), mean(xtarg_row), 'Alpha', xalpha);
                
                grp_rsa{xph}.timecourse.item_targ.cond{xcond}.blk{xblk}.pvalue_w(xrow, xcol) = xpvalue;
                grp_rsa{xph}.timecourse.item_targ.cond{xcond}.blk{xblk}.stats_w{xrow, xcol}  = xstats;
                
                %#########################################################
                %*************** baseline corrected
                %*************** similarity: pearson correlation
                xpat = grp_rsa{xph}.timecourse.item_targ_bs.cond{xcond};
                xtarg_col = []; xtarg_row = [];
                
                for xtr = it_trs
                    xtarg_col = vertcat(xtarg_col, xpat.targ{xcol}.tr{xtr}.corr_w);
                    xtarg_row = vertcat(xtarg_row, xpat.targ{xrow}.tr{xtr}.corr_w);
                end
                
                %*************** ttest
                [~, xpvalue, ~, xstats] = ttest(mean(xtarg_col), mean(xtarg_row), 'Alpha', xalpha);
                
                grp_rsa{xph}.timecourse.item_targ_bs.cond{xcond}.blk{xblk}.pvalue_w(xrow, xcol) = xpvalue;
                grp_rsa{xph}.timecourse.item_targ_bs.cond{xcond}.blk{xblk}.stats_w{xrow, xcol}  = xstats;
            end
        end
    end
end

%% *************** timecourse plot
xrep_color{1}    = [238, 20, 91]/255;% 1_item;
xrep_color{2}    = [0, 188, 182]/255;% 2_new_item
xrep_color{3}    = [160, 46, 83]/255;% 3_related_item
xrep_color{4}{2} = [34, 127, 126]/255;% 4_related_new_item
xrep_color{4}{3} = [144, 144, 144]/255;% 4_others
xrep_color{5}    = [144, 144, 144]/255;% 5_others

legend_names{2}  = {'item','new item','related','new related','others'};
legend_names{3}  = {'item','new item','related','others'};

bin_names        = {'before','after'};

f_size       = 8;
x_tick       = 1:n_trs;
x_ticklable  = x_tick;
x_lim        = [x_tick(1) x_tick(end)];
y_lim        = [-0.3 0.3];
y_tick       = y_lim(1):0.05:y_lim(end);

fig_rect     = [0 0 1400 1000];

xfig = figure;
set(xfig, 'Position', fig_rect)

for xbs = 1:2
    for it = 1:length(it_conds)
        xcond = it_conds(it);
        clear xx
        
        if xbs==1
            xx_corr = grp_rsa{xph}.timecourse.item_targ.cond{xcond};
        else
            xx_corr = grp_rsa{xph}.timecourse.item_targ_bs.cond{xcond};
        end
        
        if xcond==2, n_targs = 5; else, n_targs = 4; end
        
        subplot(2, 2, it + (2 * (xbs-1)))
        
        %*************** mean line plots
        clear xmean xse fity xpvalue
        
        fitx = linspace(1, n_trs, n_trs*10);
        
        for xtarg = 1:n_targs
            xmean{xtarg} = xx_corr.mean_w(xtarg, x_tick);
            xse{xtarg}   = xx_corr.se_w(xtarg, x_tick);
            
            fity{xtarg}  = interp1(x_tick, xmean{xtarg}, fitx,'spline');
            
            if xtarg==4
                xx_color = xrep_color{xtarg}{xcond};
            else
                xx_color = xrep_color{xtarg};
            end
            
            plot(fitx, fity{xtarg}, '-','Color', xx_color, 'LineWidth', 2); hold on;
            
        end
        
        %*************** legend
        xlegend        = legend_names{xcond};
        lg             = legend(xlegend);
        lg.Location    = 'southeast';%'bestoutside';
        lg.FontSize    = f_size;
        legend(xlegend,'AutoUpdate','off')
        grid on
        
        %*************** std error-bar filling
        for xtarg = 1:n_targs
            clear xerr fit_err
            xerr(1,:) = xmean{xtarg} - xse{xtarg};
            xerr(2,:) = xmean{xtarg} + xse{xtarg};
            
            for i = 1:2, fit_err(:,i) = interp1(x_tick, xerr(i,:), fitx,'spline'); end
            
            in_between = [fit_err(:,1)', fliplr(fit_err(:,2)')];
            
            if xtarg==4
                xx_color = xrep_color{xtarg}{xcond};
            else
                xx_color = xrep_color{xtarg};
            end
            
            h = fill([fitx fliplr(fitx)], in_between, xx_color);
            h.FaceAlpha = 0.4;
            h.EdgeAlpha = 0.4;
            h.EdgeColor = xx_color;
            
        end
        
        %*************** stats stars
        n = 0;
        for xcol = 1:(n_targs-1)
            for xrow = (xcol+1):n_targs
                yy_sig = (y_lim(2)-0.02) - (n * 0.02);
                text(-3.5, yy_sig, sprintf('%s vs. %s', xlegend{xcol}, ...
                    xlegend{xrow}), 'FontSize', f_size/1.5);
                
                for xblk = 1:n_stat_blks
                    
                    xpvalue = xx_corr.blk{xblk}.pvalue_w(xrow, xcol);
                    xx_sig  = (n_tr_blks/2) + (n_tr_blks * (xblk-1)) + 0.5 - 0.08;
                    
                    if xpvalue <= xalpha
                        text(xx_sig, yy_sig, '*', 'FontSize', f_size*1.5);
                    end
                    
                    xx = n_tr_blks + (n_tr_blks * (xblk-1)) + 0.5;
                    plot([xx xx], y_lim, '--', 'Color', 'b')
                    
                end
                
                n = n + 1;
            end
        end
        
        %*************** real onset lines
        plot([dur_stim dur_stim]+1, y_lim, '-', 'Color', xonset_color)
        plot([(dur_stim + dur_manip) (dur_stim + dur_manip)]+1, y_lim, '-', 'Color', xonset_color)
        
        h1 = text(1.5, y_lim(1)+0.05, '(1) stim onset', 'Color', xonset_color, 'FontSize', f_size);
        h2 = text(dur_stim + 1.5, y_lim(1)+0.01, '(2) operation onset', 'Color', xonset_color, 'FontSize', f_size);
        h3 = text(dur_stim + dur_manip + 1.5, y_lim(1)+0.05, '(3) fixation onset', 'Color', xonset_color, 'FontSize', f_size);
        set(h1,'Rotation', 90); set(h2,'Rotation', 90); set(h3,'Rotation', 90);
        
        %*************** set gca
        set(gca,'xlim', x_lim, 'ylim', y_lim);
        set(gca,'XTick', x_tick,'YTick', y_tick)
        set(gca,'XTickLabel', x_ticklable);
        set(gca,'FontSize', f_size/1.5)
        
        if xbs==1
            title(sprintf('RSAw item : %s (N=%s, p<%1.4f)', ...
                conds_names{xcond}, num2str(n_subjs), xalpha),...
                'FontSize', f_size*1.5,'FontWeight','bold');
        else
            title(sprintf('RSAw item corrected: %s (N=%s, p<%1.4f)', ...
                conds_names{xcond}, num2str(n_subjs), xalpha),...
                'FontSize', f_size*1.5,'FontWeight','bold');
        end
        xlabel('Volume (tr)','FontSize', f_size);
        ylabel('Similarity (Z Pearson R)','FontSize', f_size);
    end
end

%*************** save fig
fig_fname = fullfile(output_dir, sprintf('plot_tc_replace_predicted_item_%s_n%s', ...
    basename, num2str(n_subjs)));

savefig(xfig, sprintf('%s.fig', fig_fname));
set(gcf,'PaperUnits','inches','PaperPosition', fig_rect/100)
print('-dpng', '-r300', sprintf('%s.png',fig_fname))
saveas(xfig, sprintf('%s.png',fig_fname), 'png')

close(xfig);

%% *************** REPLACE BIN: PREDICTING ITEM
% xbin: item_peak: for replace: 1_before, 2_after switch
% grp_rsa{xph}.bin_mean.item_targ.cond{xcond}.bin{xbin}.corr_w(xsub, xtarg)

%*************** trial_basis: target/ non-target
for xcond = it_conds
    clear grp_mean_corr
    
    if xcond==2, n_targs = 5; else, n_targs = 4; end
    
    for xtarg = 1:n_targs
        for it_sub = 1:n_subjs
            xsub = xsubj_grp(it_sub);
            clear xr xrsa xtrial_corr % for 2 bins
            
            xrsa     = grp_rsa{xph}.subj{xsub}.timecourse.cond{xcond}.item_targ{xtarg};
            n_trials = length(xrsa.tr{1}.corr_w);
            
            rm_trials = [];
            
            for xtrial = 1:n_trials
                
                for xbin = 1:size(item_peak,1)
                    
                    xbin_corr    = [];% for each bin
                    for xtr = item_peak(xbin,1):item_peak(xbin,2)
                        if xtrial <= length(xrsa.tr{xtr}.corr_w)
                            xbin_corr     = vertcat(xbin_corr, xrsa.tr{xtr}.corr_w(xtrial));
                            
                        end
                    end
                    
                    xtrial_corr(xtrial, xbin)    = mean(xbin_corr);
                    xtrial_corr_bs(xtrial, xbin) = mean(xbin_corr - xbase_mean{it_sub}{xcond}(xtarg));
                    
                    if isempty(xbin_corr)
                        rm_trials = horzcat(rm_trials, xtrial);
                    end
                end
            end
            
            xtrial_corr(rm_trials, :)    = [];
            xtrial_corr_bs(rm_trials, :) = [];
            
            %*************** RANDOM: subj mean level
            for xbin = 1:size(item_peak,1)
                grp_mean_corr{xbin}(it_sub, xtarg)    = mean(xtrial_corr(:,xbin));
                grp_mean_corr_bs{xbin}(it_sub, xtarg) = mean(xtrial_corr_bs(:,xbin));
            end
            
        end
    end
    
    grp_rsa{xph}.bin_mean.item_targ.cond{xcond}.bin    = grp_mean_corr;
    grp_rsa{xph}.bin_mean.item_targ_bs.cond{xcond}.bin = grp_mean_corr_bs;
    
end

%% *************** ttest
for xcond = it_conds
    
    if xcond==2, n_targs = 5; else, n_targs = 4; end
    
    for xbin = 1:2
        
        grp_rsa{xph}.bin_mean.item_targ.cond{xcond}.pvalue.bin{xbin}    = nan(n_targs-1, n_targs-1);
        grp_rsa{xph}.bin_mean.item_targ_bs.cond{xcond}.pvalue.bin{xbin} = nan(n_targs-1, n_targs-1);
        
        for xcol = 1:(n_targs-1)
            for xrow = (xcol+1):n_targs
                
                xtarg_col = grp_rsa{xph}.bin_mean.item_targ.cond{xcond}.bin{xbin}(:,xcol);
                xtarg_row = grp_rsa{xph}.bin_mean.item_targ.cond{xcond}.bin{xbin}(:,xrow);
                
                %*************** ttest
                [~, xpvalue, ~, xstats] = ttest(xtarg_col, xtarg_row, 'Alpha', 0.05);
                
                grp_rsa{xph}.bin_mean.item_targ.cond{xcond}.pvalue.bin{xbin}(xrow, xcol) = xpvalue;
                grp_rsa{xph}.bin_mean.item_targ.cond{xcond}.stats.bin{xbin}{xrow, xcol}  = xstats;
                
                %%%%%%%%%%%%%%%% baseline corrected
                xtarg_col = grp_rsa{xph}.bin_mean.item_targ_bs.cond{xcond}.bin{xbin}(:,xcol);
                xtarg_row = grp_rsa{xph}.bin_mean.item_targ_bs.cond{xcond}.bin{xbin}(:,xrow);
                
                %*************** ttest
                [~, xpvalue, ~, xstats] = ttest(xtarg_col, xtarg_row, 'Alpha', 0.05);
                
                grp_rsa{xph}.bin_mean.item_targ_bs.cond{xcond}.pvalue.bin{xbin}(xrow, xcol) = xpvalue;
                grp_rsa{xph}.bin_mean.item_targ_bs.cond{xcond}.stats.bin{xbin}{xrow, xcol}  = xstats;
                
            end
        end
    end
end

%% *************** plot
%*************** repCat (2): targ:1_item, 2_new_item, 3_related_item, 4_related_new_item, 5_others
%*************** repSub (3): targ:1_item, 2_new_item, 3_related, 4_others

y_lim    = [-0.2 0.5];
fig_rect = [0 0 1400 1000];

xfig = figure;
set(xfig, 'Position', fig_rect)

for i = 1:length(it_conds)
    
    xcond = it_conds(i);
    
    if xcond==2, n_targs = 5; else, n_targs = 4; end
    
    for xbin = 1:size(item_peak,1)
        clear b xrsa_corr
        
        subplot(length(it_conds), 2, xbin + 2 * (i-1))
        
        xrsa_corr = grp_rsa{xph}.bin_mean.item_targ_bs.cond{xcond}.bin{xbin};
        xmean = mean(xrsa_corr);
        xse   = std(xrsa_corr)/sqrt(n_subjs);
        
        for xtarg = 1:n_targs
            b{xtarg} = bar(xtarg, xmean(xtarg)); hold on
            
            if xtarg==4
                set(b{xtarg}, 'facecolor', xrep_color{xtarg}{xcond});
            else
                set(b{xtarg}, 'facecolor', xrep_color{xtarg});
            end
        end
        
        for xtarg = 1:n_targs
            errorbar(xtarg, xmean(xtarg), xse(xtarg),'k.')
        end
        
        %*************** ttest
        n = 1;
        for xcol = 1:(n_targs-1)
            for xrow = (xcol+1):n_targs
                xpvalue = grp_rsa{xph}.bin_mean.item_targ_bs.cond{xcond}.pvalue.bin{xbin}(xrow, xcol);
                xstats  = grp_rsa{xph}.bin_mean.item_targ_bs.cond{xcond}.stats.bin{xbin}{xrow, xcol};
                
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
        
        set(gca,'xlim',[1 n_targs]+[-0.5 0.5],'ylim', y_lim);
        set(gca,'XTick', 1:n_targs, 'YTick', y_lim(1):0.1:y_lim(end));
        set(gca,'XTickLabel', legend_names{xcond})
        
        title(sprintf('%s: %s replacing: %s-%s', ...
            conds_names{xcond}, bin_names{xbin}, ...
            num2str(item_peak(xbin,1)), num2str(item_peak(xbin,2))));
        ylabel('RSA');
    end
end

%*************** save fig
fig_fname = fullfile(output_dir, sprintf('plot_replace_predicted_item_%s_n%s', ...
    basename, num2str(n_subjs)));

savefig(xfig, sprintf('%s.fig', fig_fname));
set(gcf,'PaperUnits','inches','PaperPosition', fig_rect/100)
print('-dpng', '-r300', sprintf('%s.png',fig_fname))
saveas(xfig, sprintf('%s.png',fig_fname), 'png')

close(xfig);

%% *************** plot: only target vs. new target
%*************** repCat (2): targ:1_item, 2_new_item, 3_related_item, 4_related_new_item, 5_others
%*************** repSub (3): targ:1_item, 2_new_item, 3_related, 4_others

y_lim    = [-0.2 0.5];
fig_rect = [0 0 1400 1000];
n_targs   = 2;

xfig = figure;
set(xfig, 'Position', fig_rect)

for xbin = 1:size(item_peak,1)
    clear b xrsa_corr
    
    subplot(2, 2, xbin)
    
    xrsa_corr = [];
    for xcond = it_conds
        xrsa_corr = horzcat(xrsa_corr, ...
            grp_rsa{xph}.bin_mean.item_targ_bs.cond{xcond}.bin{xbin}(:,1:2));
    end
    
    xmean = mean(xrsa_corr);
    xse   = std(xrsa_corr)/sqrt(n_subjs);
    
    for xtarg = 1:4
        b = bar(xtarg, xmean(xtarg)); hold on
        
        if xtarg < 3
            b.FaceColor = xrep_color{mod(xtarg-1,2)+1};
        else
            b.FaceColor = 'w';
            b.EdgeColor = xrep_color{mod(xtarg-1,2)+1};
            b.LineWidth = 3;
        end
        
    end
    
    for xtarg = 1:4
        errorbar(xtarg, xmean(xtarg), xse(xtarg),'k.')
    end
    
    %*************** ttest
    n = 1;
    for xcol = 1:(4-1)
        for xrow = (xcol+1):4
            
            %%%%%%%%%%%%%%%% baseline corrected
            xtarg_col = xrsa_corr(:,xcol);
            xtarg_row = xrsa_corr(:,xrow);
            
            %*************** ttest
            [~, xpvalue, ~, xstats] = ttest(xtarg_col, xtarg_row, 'Alpha', 0.05);
            
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
    xlegend        = {'repCat: item','repCat: new item','repSubcat: item','repSubcat: new item'};
    lg             = legend(xlegend);
    lg.Location    = 'southeastoutside';%'SouthEast';
    legend(xlegend,'AutoUpdate','off')
    
    set(gca,'xlim',[1 4]+[-0.5 0.5],'ylim', y_lim);
    set(gca,'XTick', 1:4, 'YTick', y_lim(1):0.1:y_lim(end));
    set(gca,'XTickLabel', xlegend)
    
    title(sprintf('%s: %s replacing: %s-%s', ...
        conds_names{xcond}, bin_names{xbin}, ...
        num2str(item_peak(xbin,1)), num2str(item_peak(xbin,2))));
    ylabel('RSA');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%*************** timecourse
x_tick       = 1:n_trs;
x_ticklable  = x_tick;
x_lim        = [x_tick(1) x_tick(end)];
y_lim        = [-0.3 0.3];
y_tick       = y_lim(1):0.05:y_lim(end);

for xtarg = 1:2
    clear xmean xse fity xpvalue fitx
    
    subplot(2, 2, 2 + xtarg)
    
    %*************** mean line plots
    fitx  = linspace(1, n_trs, n_trs*10);
    for xcond = it_conds
        xmean = grp_rsa{xph}.timecourse.item_targ_bs.cond{xcond}.mean_w(xtarg, x_tick);
        xse   = grp_rsa{xph}.timecourse.item_targ_bs.cond{xcond}.se_w(xtarg, x_tick);
        
        fity  = interp1(x_tick, xmean, fitx,'spline');
        
        if xcond == 2
            plot(fitx, fity, '-','Color', xrep_color{xtarg}, 'LineWidth', 2); hold on;
        else
            plot(fitx, fity, '--','Color', xrep_color{xtarg}, 'LineWidth', 2); hold on;
        end
    end
    
    %*************** legend
    xlegend = {'replace (category)','replace (subcategory)'};
    lg             = legend(xlegend);
    lg.Location    = 'southeast';%'bestoutside';
    legend(xlegend,'AutoUpdate','off')
    grid on
    
    %*************** std error-bar filling
    for xcond = it_conds
        clear xerr fit_err
        xmean = grp_rsa{xph}.timecourse.item_targ_bs.cond{xcond}.mean_w(xtarg, x_tick);
        xse   = grp_rsa{xph}.timecourse.item_targ_bs.cond{xcond}.se_w(xtarg, x_tick);
        
        xerr(1,:) = xmean - xse;
        xerr(2,:) = xmean + xse;
        
        for i = 1:2, fit_err(:,i) = interp1(x_tick, xerr(i,:), fitx,'spline'); end
        
        in_between = [fit_err(:,1)', fliplr(fit_err(:,2)')];
        
        h = fill([fitx fliplr(fitx)], in_between, xrep_color{xtarg});
        h.FaceAlpha = 0.2;
        h.EdgeAlpha = 0.2;
        h.EdgeColor = xrep_color{xtarg};
    end
    
    %*************** stats stars
    for xblk = 1:n_stat_blks
        
        it_trs = (1:n_tr_blks) + (n_tr_blks * (xblk-1));
        
        yy_sig = (y_lim(2)-0.02);
        xpat   = grp_rsa{xph}.timecourse.item_targ_bs;
        
        xtarg_col = []; xtarg_row = [];
        
        for xtr = it_trs
            xtarg_col = vertcat(xtarg_col, xpat.cond{2}.targ{xtarg}.tr{xtr}.corr_w);
            xtarg_row = vertcat(xtarg_row, xpat.cond{3}.targ{xtarg}.tr{xtr}.corr_w);
        end
        
        %*************** ttest
        [~, xpvalue, ~, xstats] = ttest(mean(xtarg_col), mean(xtarg_row), 'Alpha', xalpha);
        
        xx_sig = (n_tr_blks/2) + (n_tr_blks * (xblk-1)) + 0.5 - 0.08;
        
        if xpvalue <= xalpha
            text(xx_sig, yy_sig, '*', 'FontSize', f_size*1.5);
        end
        
        xx = n_tr_blks + (n_tr_blks * (xblk-1)) + 0.5;
        plot([xx xx], y_lim, '--', 'Color', 'b')
    end
    
    %*************** real onset lines
    plot([dur_stim dur_stim]+1, y_lim, '-', 'Color', xonset_color)
    plot([(dur_stim + dur_manip) (dur_stim + dur_manip)]+1, y_lim, '-', 'Color', xonset_color)
    
    h1 = text(1.5, y_lim(1)+0.05, '(1) stim onset', 'Color', xonset_color, 'FontSize', f_size);
    h2 = text(dur_stim + 1.5, y_lim(1)+0.01, '(2) operation onset', 'Color', xonset_color, 'FontSize', f_size);
    h3 = text(dur_stim + dur_manip + 1.5, y_lim(1)+0.05, '(3) fixation onset', 'Color', xonset_color, 'FontSize', f_size);
    set(h1,'Rotation', 90); set(h2,'Rotation', 90); set(h3,'Rotation', 90);
    
    %*************** set gca
    set(gca,'xlim', x_lim, 'ylim', y_lim);
    set(gca,'XTick', x_tick,'YTick', y_tick)
    set(gca,'XTickLabel', x_ticklable);
    set(gca,'FontSize', f_size/1.5)
    
    title(sprintf('RSAw item corrected: %s (N=%s, p<%1.4f)', ...
        legend_names{2}{xtarg}, num2str(n_subjs), xalpha),...
        'FontSize', f_size*1.5,'FontWeight','bold');
    
    xlabel('Volume (tr)','FontSize', f_size);
    ylabel('Similarity (Z Pearson R)','FontSize', f_size);
    
end

%*************** save fig
fig_fname = fullfile(output_dir, sprintf('plot_comb_replace_predicted_item_%s_n%s', ...
    basename, num2str(n_subjs)));

savefig(xfig, sprintf('%s.fig', fig_fname));
set(gcf,'PaperUnits','inches','PaperPosition', fig_rect/100)
print('-dpng', '-r300', sprintf('%s.png',fig_fname))
saveas(xfig, sprintf('%s.png',fig_fname), 'png')

close(xfig);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%============= STUDY TIMECOURSE: PER TARGET
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xph = 2;
n_targs = 3;
%*************** mean: target/ non-target
for xcond = 1:n_condition
    
    grp_rsa{xph}.timecourse.cond{xcond}.mean_w = zeros(n_targs, n_trs);
    grp_rsa{xph}.timecourse.cond{xcond}.se_w   = zeros(n_targs, n_trs);
    
    for xtarg = 1:n_targs
        for xtr = 1:n_trs
            
            grp_rsa{xph}.timecourse.cond{xcond}.mean_w(xtarg, xtr) = ...
                mean(grp_rsa{xph}.timecourse.cond{xcond}.targ{xtarg}.tr{xtr}.corr_w);
            grp_rsa{xph}.timecourse.cond{xcond}.se_w(xtarg, xtr) = ...
                std(grp_rsa{xph}.timecourse.cond{xcond}.targ{xtarg}.tr{xtr}.corr_w)/sqrt(n_subjs);
            
            %*************** same N, N+1 category/subcategory
            for xlevel = 1:2
                
                if xlevel==1, it_sames = 1:2; else, it_sames = 1:3; end
                
                for xsame = it_sames % 1_same, 2_differ, 3_related
                    grp_rsa{xph}.timecourse.sameseq{xlevel}{xsame}.cond{xcond}.mean_w(xtarg, xtr) = ...
                        mean(grp_rsa{xph}.timecourse.sameseq{xlevel}{xsame}.cond{xcond}.targ{xtarg}.tr{xtr}.corr_w);
                    grp_rsa{xph}.timecourse.sameseq{xlevel}{xsame}.cond{xcond}.se_w(xtarg, xtr) = ...
                        std(grp_rsa{xph}.timecourse.sameseq{xlevel}{xsame}.cond{xcond}.targ{xtarg}.tr{xtr}.corr_w)/sqrt(n_subjs);
                end
            end
        end
        
        %*************** sync TRs
        for xtr = 1:dur_sync
            
            grp_rsa{xph}.timecourse.sync.cond{xcond}.mean_w(xtarg, xtr) = ...
                mean(grp_rsa{xph}.timecourse.sync.cond{xcond}.targ{xtarg}.tr{xtr}.corr_w);
            grp_rsa{xph}.timecourse.sync.cond{xcond}.se_w(xtarg, xtr) = ...
                std(grp_rsa{xph}.timecourse.sync.cond{xcond}.targ{xtarg}.tr{xtr}.corr_w)/sqrt(n_subjs);
            
            %*************** same N, N+1 category/subcategory
            for xlevel = 1:2
                
                if xlevel==1, it_sames = 1:2; else, it_sames = 1:3; end
                
                for xsame = it_sames % 1_same, 2_differ, 3_related
                    grp_rsa{xph}.timecourse.sameseq{xlevel}{xsame}.sync.cond{xcond}.mean_w(xtarg, xtr) = ...
                        mean(grp_rsa{xph}.timecourse.sameseq{xlevel}{xsame}.sync.cond{xcond}.targ{xtarg}.tr{xtr}.corr_w);
                    grp_rsa{xph}.timecourse.sameseq{xlevel}{xsame}.sync.cond{xcond}.se_w(xtarg, xtr) = ...
                        std(grp_rsa{xph}.timecourse.sameseq{xlevel}{xsame}.sync.cond{xcond}.targ{xtarg}.tr{xtr}.corr_w)/sqrt(n_subjs);
                    
                    %*************** next
                    grp_rsa{xph}.timecourse.sameseq{xlevel}{xsame}.sync.next.cond{xcond}.mean_w(xtarg, xtr) = ...
                        mean(grp_rsa{xph}.timecourse.sameseq{xlevel}{xsame}.sync.next.cond{xcond}.targ{xtarg}.tr{xtr}.corr_w);
                    grp_rsa{xph}.timecourse.sameseq{xlevel}{xsame}.sync.next.cond{xcond}.se_w(xtarg, xtr) = ...
                        std(grp_rsa{xph}.timecourse.sameseq{xlevel}{xsame}.sync.next.cond{xcond}.targ{xtarg}.tr{xtr}.corr_w)/sqrt(n_subjs);
                end
            end
        end
    end
end

%% ============= STATS: RANDOM
for xcond = 1:n_condition
    for xblk = 1:n_stat_blks
        
        grp_rsa{xph}.timecourse.cond{xcond}.blk{xblk}.pvalue_w = nan(n_targs, n_targs);
        
        for xcol = 1:(n_targs-1)
            for xrow = (xcol+1):n_targs
                clear xpat
                
                it_trs = (1:n_tr_blks) + (n_tr_blks * (xblk-1));
                
                xpat = grp_rsa{xph}.timecourse.cond{xcond};
                
                for xx = 2%:2%1_raw, 2_weight
                    %#########################################################
                    %*************** similarity: pearson correlation
                    xtarg_col = []; xtarg_row = [];
                    
                    for xtr = it_trs
                        if xx == 1
                            xtarg_col = vertcat(xtarg_col, xpat.targ{xcol}.tr{xtr}.corr);
                            xtarg_row = vertcat(xtarg_row, xpat.targ{xrow}.tr{xtr}.corr);
                        elseif xx == 2
                            xtarg_col = vertcat(xtarg_col, xpat.targ{xcol}.tr{xtr}.corr_w);
                            xtarg_row = vertcat(xtarg_row, xpat.targ{xrow}.tr{xtr}.corr_w);
                        end
                    end
                    
                    %*************** ttest
                    [~, xpvalue, ~, xstats] = ttest(mean(xtarg_col), mean(xtarg_row), 'Alpha', xalpha);
                    
                    if xx == 1
                        grp_rsa{xph}.timecourse.cond{xcond}.blk{xblk}.pvalue(xrow, xcol)   = xpvalue;
                        grp_rsa{xph}.timecourse.cond{xcond}.blk{xblk}.stats{xrow, xcol}    = xstats;
                    elseif xx == 2
                        grp_rsa{xph}.timecourse.cond{xcond}.blk{xblk}.pvalue_w(xrow, xcol) = xpvalue;
                        grp_rsa{xph}.timecourse.cond{xcond}.blk{xblk}.stats_w{xrow, xcol}  = xstats;
                    end
                end
            end
        end
    end
end

%% ============= SAME N, N+1 CATEGORY: STATS: RANDOM:
for xcond = 1:n_condition
    for xblk = 1:n_stat_blks
        for xlevel = 1:2
            
            if xlevel==1, it_sames = 1:2; else, it_sames = 1:3; end
            
            for xsame = it_sames % 1_same, 2_differ, 3_related
                
                grp_rsa{xph}.timecourse.sameseq{xlevel}{xsame}.cond{xcond}.blk{xblk}.pvalue_w = nan(n_targs, n_targs);
                
                for xcol = 1:(n_targs-1)
                    for xrow = (xcol+1):n_targs
                        clear xpat
                        
                        it_trs = (1:n_tr_blks) + (n_tr_blks * (xblk-1));
                        
                        xpat = grp_rsa{xph}.timecourse.sameseq{xlevel}{xsame}.cond{xcond};
                        
                        %#########################################################
                        %*************** similarity: pearson correlation
                        xtarg_col = []; xtarg_row = [];
                        
                        for xtr = it_trs
                            xtarg_col = vertcat(xtarg_col, xpat.targ{xcol}.tr{xtr}.corr_w);
                            xtarg_row = vertcat(xtarg_row, xpat.targ{xrow}.tr{xtr}.corr_w);
                        end
                        
                        %*************** ttest
                        [~, xpvalue, ~, xstats] = ttest(mean(xtarg_col), mean(xtarg_row), 'Alpha', xalpha);
                        
                        grp_rsa{xph}.timecourse.sameseq{xlevel}{xsame}.cond{xcond}.blk{xblk}.pvalue_w(xrow, xcol) = xpvalue;
                        grp_rsa{xph}.timecourse.sameseq{xlevel}{xsame}.cond{xcond}.blk{xblk}.stats_w{xrow, xcol}  = xstats;
                        
                    end
                end
            end
        end
    end
end

%% ============= STATS: B/W CONDS WITH ONLY-TARGET
clear xpat
xpat = grp_rsa{xph}.timecourse;

for xblk = 1:n_stat_blks
    
    grp_rsa{xph}.timecourse.targs.blk{xblk}.pvalue_w = nan(n_condition, n_condition);
    
    for xcol = 1:(n_condition-1)
        for xrow = (xcol+1):n_condition
            it_trs = (1:n_tr_blks) + (n_tr_blks * (xblk-1));
            
            for xx = 2%1:2
                %*************** similarity: pearson correlation\
                xtarg_col = []; xtarg_row = [];
                
                for xtr = it_trs
                    if xx == 1
                        xtarg_col = vertcat(xtarg_col, xpat.cond{xcol}.targ{1}.tr{xtr}.corr);
                        xtarg_row = vertcat(xtarg_row, xpat.cond{xrow}.targ{1}.tr{xtr}.corr);
                    elseif xx == 2
                        xtarg_col = vertcat(xtarg_col, xpat.cond{xcol}.targ{1}.tr{xtr}.corr_w);
                        xtarg_row = vertcat(xtarg_row, xpat.cond{xrow}.targ{1}.tr{xtr}.corr_w);
                    end
                end
                
                %*************** ttest
                [~, xpvalue, ~, xstats] = ttest(mean(xtarg_col), mean(xtarg_row), 'Alpha', xalpha);
                
                if xx == 1
                    grp_rsa{xph}.timecourse.targs.blk{xblk}.pvalue(xrow, xcol)   = xpvalue;
                    grp_rsa{xph}.timecourse.targs.blk{xblk}.stats{xrow, xcol}    = xstats;
                elseif xx == 2
                    grp_rsa{xph}.timecourse.targs.blk{xblk}.pvalue_w(xrow, xcol) = xpvalue;
                    grp_rsa{xph}.timecourse.targs.blk{xblk}.stats_w{xrow, xcol}  = xstats;
                end
            end
        end
    end
end

%% ============= SAME N, N+1 CATEGORY: STATS: B/W CONDS WITH ONLY-TARGET
for xlevel = 1:2
    
    if xlevel==1, it_sames = 1:2; else, it_sames = 1:3; end
    
    for xsame = it_sames % 1_same, 2_differ, 3_related
        clear xpat
        xpat = grp_rsa{xph}.timecourse.sameseq{xlevel}{xsame};
        
        for xblk = 1:n_stat_blks
            
            grp_rsa{xph}.timecourse.sameseq{xlevel}{xsame}.targs.blk{xblk}.pvalue_w = nan(n_condition, n_condition);
            
            for xcol = 1:(n_condition-1)
                for xrow = (xcol+1):n_condition
                    it_trs = (1:n_tr_blks) + (n_tr_blks * (xblk-1));
                    
                    %*************** similarity: pearson correlation
                    xtarg_col = []; xtarg_row = [];
                    
                    for xtr = it_trs
                        xtarg_col = vertcat(xtarg_col, xpat.cond{xcol}.targ{1}.tr{xtr}.corr_w);
                        xtarg_row = vertcat(xtarg_row, xpat.cond{xrow}.targ{1}.tr{xtr}.corr_w);
                    end
                    
                    %*************** ttest
                    [~, xpvalue, ~, xstats] = ttest(mean(xtarg_col), mean(xtarg_row), 'Alpha', xalpha);
                    
                    grp_rsa{xph}.timecourse.sameseq{xlevel}{xsame}.targs.blk{xblk}.pvalue_w(xrow, xcol) = xpvalue;
                    grp_rsa{xph}.timecourse.sameseq{xlevel}{xsame}.targs.blk{xblk}.stats_w{xrow, xcol}  = xstats;
                end
            end
        end
    end
end

%% ============= SYNC STATS: B/W CONDS WITH ONLY-TARGET
clear xpat
xpat = grp_rsa{xph}.timecourse.sync;

for xblk = 1:(dur_sync/n_tr_blks)
    
    grp_rsa{xph}.timecourse.sync.targs.blk{xblk}.pvalue_w = nan(n_condition, n_condition);
    
    for xcol = 1:(n_condition-1)
        for xrow = (xcol+1):n_condition
            it_trs = (1:n_tr_blks) + (n_tr_blks * (xblk-1));
            
            for xx = 2%1:2
                %*************** similarity: pearson correlation
                xtarg_col = []; xtarg_row = [];
                
                for xtr = it_trs
                    if xx == 1
                        xtarg_col = vertcat(xtarg_col, xpat.cond{xcol}.targ{1}.tr{xtr}.corr);
                        xtarg_row = vertcat(xtarg_row, xpat.cond{xrow}.targ{1}.tr{xtr}.corr);
                    elseif xx == 2
                        xtarg_col = vertcat(xtarg_col, xpat.cond{xcol}.targ{1}.tr{xtr}.corr_w);
                        xtarg_row = vertcat(xtarg_row, xpat.cond{xrow}.targ{1}.tr{xtr}.corr_w);
                    end
                end
                
                %*************** ttest
                [~, xpvalue, ~, xstats] = ttest(mean(xtarg_col), mean(xtarg_row), 'Alpha', xalpha);
                
                if xx == 1
                    grp_rsa{xph}.timecourse.sync.targs.blk{xblk}.pvalue(xrow, xcol)   = xpvalue;
                    grp_rsa{xph}.timecourse.sync.targs.blk{xblk}.stats{xrow, xcol}    = xstats;
                elseif xx == 2
                    grp_rsa{xph}.timecourse.sync.targs.blk{xblk}.pvalue_w(xrow, xcol) = xpvalue;
                    grp_rsa{xph}.timecourse.sync.targs.blk{xblk}.stats_w{xrow, xcol}  = xstats;
                end
            end
        end
    end
end

%% ============= SAME N, N+1 CATEGORY: SYNC STATS: B/W CONDS WITH ONLY-TARGET
for xlevel = 1:2
    
    if xlevel==1, it_sames = 1:2; else, it_sames = 1:3; end
    
    for xsame = it_sames
        clear xpat
        
        xpat = grp_rsa{xph}.timecourse.sameseq{xlevel}{xsame}.sync;
        
        for xblk = 1:(dur_sync/n_tr_blks)
            
            grp_rsa{xph}.timecourse.sameseq{xlevel}{xsame}.sync.targs.blk{xblk}.pvalue_w = nan(n_condition, n_condition);
            
            for xcol = 1:(n_condition-1)
                for xrow = (xcol+1):n_condition
                    it_trs = (1:n_tr_blks) + (n_tr_blks * (xblk-1));
                    
                    for xx = 2
                        %*************** similarity: pearson correlation\
                        xtarg_col = []; xtarg_row = [];
                        
                        for xtr = it_trs
                            xtarg_col = vertcat(xtarg_col, xpat.cond{xcol}.targ{1}.tr{xtr}.corr_w);
                            xtarg_row = vertcat(xtarg_row, xpat.cond{xrow}.targ{1}.tr{xtr}.corr_w);
                        end
                        
                        %*************** ttest
                        [~, xpvalue, ~, xstats] = ttest(mean(xtarg_col), mean(xtarg_row), 'Alpha', xalpha);
                        
                        grp_rsa{xph}.timecourse.sameseq{xlevel}{xsame}.sync.targs.blk{xblk}.pvalue_w(xrow, xcol) = xpvalue;
                        grp_rsa{xph}.timecourse.sameseq{xlevel}{xsame}.sync.targs.blk{xblk}.stats_w{xrow, xcol}  = xstats;
                    end
                end
            end
        end
    end
end

%% ============= SAME N, N+1 CATEGORY: NEXT: SYNC STATS: B/W CONDS WITH ONLY-TARGET
for xlevel = 1:2
    
    if xlevel==1, it_sames = 1:2; else, it_sames = 1:3; end
    
    for xsame = it_sames
        clear xpat
        
        xpat = grp_rsa{xph}.timecourse.sameseq{xlevel}{xsame}.sync.next;
        
        for xblk = 1:(dur_sync/n_tr_blks)
            
            grp_rsa{xph}.timecourse.sameseq{xlevel}{xsame}.sync.next.targs.blk{xblk}.pvalue_w = nan(n_condition, n_condition);
            
            for xcol = 1:(n_condition-1)
                for xrow = (xcol+1):n_condition
                    it_trs = (1:n_tr_blks) + (n_tr_blks * (xblk-1));
                    
                    for xx = 2
                        %*************** similarity: pearson correlation\
                        xtarg_col = []; xtarg_row = [];
                        
                        for xtr = it_trs
                            xtarg_col = vertcat(xtarg_col, xpat.cond{xcol}.targ{1}.tr{xtr}.corr_w);
                            xtarg_row = vertcat(xtarg_row, xpat.cond{xrow}.targ{1}.tr{xtr}.corr_w);
                        end
                        
                        %*************** ttest
                        [~, xpvalue, ~, xstats] = ttest(mean(xtarg_col), mean(xtarg_row), 'Alpha', xalpha);
                        
                        grp_rsa{xph}.timecourse.sameseq{xlevel}{xsame}.sync.next.targs.blk{xblk}.pvalue_w(xrow, xcol) = xpvalue;
                        grp_rsa{xph}.timecourse.sameseq{xlevel}{xsame}.sync.next.targs.blk{xblk}.stats_w{xrow, xcol}  = xstats;
                    end
                end
            end
        end
    end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%============= STUDY BIN: N VS. N+1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% xbin:   1_n trial (args.bin_trs{1}), 2_n+1 trial (args.bin_trs{2}), 3_n+1 baseline
% xtarg:  1_target, 2_related_nontarg 2_nontarg
% xlevel: 1_category, 2_subcategory
% xsame:  1_same, 2_differ, 3_related

xph    = 2;
xlevel = 1;% 1_category, 2_subcategory

if xlevel==1, it_sames = 1:2; else, it_sames = 1:3; end

for xcond = 1:n_condition
    for xtarg = 1:n_targs
        for xsub = xsubj_grp
            grp_rsa{xph}.sub{xsub}.bin_dist.cond{xcond}.targ{xtarg}.corr_w = [];
        end
    end
end

%*************** trial_basis: target/ non-target
for xcond = 1:n_condition
    for xtarg = 1:n_targs
        for xsame = it_sames
            
            grp_trial_corr = [];
            grp_mean_corr  = [];
            grp_xr         = [];
            
            for it_sub = 1:n_subjs
                xsub = xsubj_grp(it_sub);
                clear xr xrsa xtrial_corr % for 2 bins
                
                for xbin = 1:length(args.bin_trs)
                    if xbin==1
                        xrsa{xbin} = grp_rsa{xph}.subj{xsub}.timecourse.sameseq{xlevel}{xsame}.cond{xcond}.targ{xtarg};
                    else
                        xrsa{xbin} = grp_rsa{xph}.subj{xsub}.timecourse.sameseq{xlevel}{xsame}.sync.next.cond{xcond}.targ{xtarg};
                    end
                end
                
                n_trials = length(xrsa{1}.tr{1}.corr_w);
                
                rm_trials = [];
                
                for xtrial = 1:n_trials
                    
                    for xbin = 1:3
                        
                        xbin_corr = [];% for each bin
                        
                        for xtr = args.bin_trs{xbin}
                            if xtrial <= length(xrsa{xbin}.tr{xtr}.corr_w)
                                xbin_corr = vertcat(xbin_corr, xrsa{xbin}.tr{xtr}.corr_w(xtrial));
                                
                                if xbin==1
                                    %*************** DISTRIBUTION
                                    grp_rsa{xph}.sub{xsub}.bin_dist.cond{xcond}.targ{xtarg}.corr_w = ...
                                        vertcat(grp_rsa{xph}.sub{xsub}.bin_dist.cond{xcond}.targ{xtarg}.corr_w,...
                                        xrsa{xbin}.tr{xtr}.corr_w(xtrial));
                                end
                            end
                        end
                        
                        xtrial_corr(xtrial, xbin) = mean(xbin_corr);
                        
                        if isempty(xbin_corr)
                            rm_trials = horzcat(rm_trials, xtrial);
                        end
                    end
                end
                
                xtrial_corr(rm_trials, :) = [];
                
                %*************** FIXED: trial level
                grp_trial_corr = vertcat(grp_trial_corr, xtrial_corr);
                grp_rsa{xph}.sub{xsub}.bin_trial.sameseq{xlevel}.cond{xcond}.targ{xtarg}.same{xsame}.corr_w = xtrial_corr;
                
                %*************** RANDOM: subj mean level
                grp_mean_corr(it_sub, :) = mean(xtrial_corr);
                
                %*************** CORRELATION RANDOM
                xr = corrcoef(xtrial_corr(:,1), xtrial_corr(:,2));
                grp_xr(it_sub) = xr(1,2);
                
            end
            
            grp_rsa{xph}.bin_trial.sameseq{xlevel}.cond{xcond}.targ{xtarg}.same{xsame}.corr_w = grp_trial_corr;
            grp_rsa{xph}.bin_mean.sameseq{xlevel}.cond{xcond}.targ{xtarg}.same{xsame}.corr_w  = grp_mean_corr;
            
            grp_rsa{xph}.bin_trial.sameseq{xlevel}.cond{xcond}.targ{xtarg}.same{xsame}.r = grp_xr;
            
        end
    end
end

%% ============= REPLACE N+1 SAME/DIFF BASED ON NEW-ITEM
% grp_rsa{xph}.bin_mean.sameseq{xlevel}.cond{xcond}.newitem{2}.targ{xtarg}.same{2}.corr_w

xlevel = 1;%1_cate
xcond  = 2;%2_repCat
xsame  = 2;%2_diff
xtarg  = 1;

clear grp_mean_corr

for it_sub = 1:n_subjs
    xsub = xsubj_grp(it_sub);
    clear xrsa
    
    xrsa = grp_rsa{xph}.subj{xsub}.timecourse.sameseq{xlevel}{xsame}.sync.next.cond{xcond};
    
    for xnew_same = 1:2
        clear xbin_corr bs_corr
        
        for xbin = 2:3
            
%             %*************** baseline correction
%             for xtr = 1:args.baseline_trs
%                 bs_corr(xtr) = xrsa.newitem{xnew_same}.targ{xtarg}.tr{xtr}.corr_w;
%             end
%             
            xbin_corr = [];% for each bin
            
            for xtr = args.bin_trs{xbin}
                xbin_corr = vertcat(xbin_corr, ...
                    mean(xrsa.newitem{xnew_same}.targ{xtarg}.tr{xtr}.corr_w));
            end
            
            %*************** RANDOM: subj mean level
            grp_mean_corr{xnew_same}(it_sub, xbin) = mean(xbin_corr);
            
            if xbin==3
                n_newsame{xnew_same}(it_sub) = length(xrsa.newitem{xnew_same}.targ{xtarg}.tr{xtr}.corr_w);
            end
        end
        
        grp_rsa{xph}.bin_mean.sameseq{xlevel}.cond{xcond}.newitem{xnew_same}.targ{xtarg}.same{xsame}.corr_w  = grp_mean_corr{xnew_same};
    end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%============= STUDY BIN PLOT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
level_names = {'category','subcategory'};
same_names  = {'same','differ','related'};

f_size       = 10;

%% ============= DISTRIBUTION
xtarg      = 1;
n_bin      = 5;

x_lim      = [0 n_condition+1];
y_lim      = [0 1/(n_bin/2)];

min_x      = 0;
max_x      = 1;

bin_width  = (max_x - min_x)/n_bin;
x_bins     = min_x:bin_width:max_x;
x_tick     = 1:n_condition;
y_tick     = y_lim(1):0.05:y_lim(2);

point_size = 8;

rect_w     = 600;
rect_h     = 500;
fig_rect   = [0 0 rect_w*n_bin rect_h];

clear y_prob_mean y_prob_se sub_prob
y_prob_mean = zeros(n_bin, n_condition);
y_prob_se   = zeros(n_bin, n_condition);

for xcond = 1:n_condition
    for it_sub = 1:n_subjs
        xsub = xsubj_grp(it_sub);
        
        xx_evidence = [];
        xx_evidence = grp_rsa{xph}.sub{xsub}.bin_dist.cond{xcond}.targ{xtarg}.corr_w;
        sub_prob{xcond}(it_sub, 1:n_bin) = histcounts(xx_evidence, x_bins)/length(xx_evidence);
    end
    
    y_prob_mean(:, xcond) = mean(sub_prob{xcond});
    y_prob_se(:, xcond)   = std(sub_prob{xcond})/sqrt(n_subjs);
end

xfig = figure;
set(xfig, 'Position', fig_rect)

for xbin = 1:n_bin
    
    subplot(1, n_bin, xbin)
    
    for xcond = 1:n_condition
        b = bar(xcond, y_prob_mean(xbin, xcond)); hold on;
        b.FaceColor = xcond_color{xcond};
    end
    
    %*************** legend
    if xbin==n_bin
        xlegend        = conds_names;
        lg             = legend(xlegend);
        lg.Location    = 'NorthEast';%'SouthEast';
        legend(xlegend,'AutoUpdate','off')
    end
    
    for xcond = 1:n_condition
        errorbar(xcond, y_prob_mean(xbin, xcond), y_prob_se(xbin, xcond),...
            '-k','LineWidth',1);
    end
    
    set(gca,'XLim', x_lim, 'YLim',y_lim, 'XTick',x_tick, 'YTick',y_tick);
    h=set(gca,'XTickLabel', conds_names, 'FontSize', f_size);
    xtickangle(45)
    
    title(sprintf('%s-%s TR, bin %s: %4.4f-%4.4f RSA', ...
        num2str(args.bin_trs{1}(1)), num2str(args.bin_trs{1}(end)),...
        num2str(xbin), x_bins(xbin), x_bins(xbin+1)));
    ylabel('proportion');
end

%*************** save fig
fig_fname = fullfile(output_dir, sprintf('plot_dist_bin%s_%s_n%s', ...
    num2str(n_bin), basename, num2str(n_subjs)));

%     savefig(xfig, sprintf('%s.fig', fig_fname));
set(gcf,'PaperUnits','inches','PaperPosition', fig_rect/100)
print('-dpng', '-r300', sprintf('%s.png',fig_fname))
saveas(xfig, sprintf('%s.png',fig_fname), 'png')

close(xfig);

%% ============= RANDOM: NEXT-TRIAL BIN: TARGET-ONLY
clear xmean xse xpvalue

xlevel     = 1;%1_cate
it_bins    = [2 3];
xtarg      = 1;
xbin_alpha = 0.05;

if xlevel==1, it_sames = 1:2; else, it_sames = 1:3; end

clear xpvalue xrsa xstats it_rsa

for it_bin = 1:length(it_bins)
    xbin = it_bins(it_bin);
    
    for xcond = 1:n_condition
        %*************** mean/se
        for xsame = it_sames
            
            if (xcond == 2) && (xsame == 2)
                xnew_same = 2;%only when new-item and next-item is different
                it_rsa = grp_rsa{xph}.bin_mean.sameseq{xlevel}.cond{xcond}.newitem{xnew_same}.targ{xtarg};
            else
                it_rsa = grp_rsa{xph}.bin_mean.sameseq{xlevel}.cond{xcond}.targ{xtarg};
            end
            
            xmean{it_bin}(xsame, xcond) = mean(it_rsa.same{xsame}.corr_w(:, xbin));
            xse{it_bin}(xsame, xcond)   = std(it_rsa.same{xsame}.corr_w(:, xbin))/sqrt(n_subjs);
            
            xrsa.same{xsame}.corr_w = it_rsa.same{xsame}.corr_w;
            
        end
        
        %*************** ttest
        [~, xpvalue{it_bin}(xcond), ~, xstats{it_bin}{xcond}] = ttest(xrsa.same{1}.corr_w(:, xbin), xrsa.same{2}.corr_w(:, xbin), 'Alpha', xbin_alpha);
        
    end
end

%%*************** figure
f_size       = 10;
x_ticklable  = same_names(it_sames);
x_lim        = [1 length(it_sames)]+[-0.5 0.5];
y_lim        = [-0.05 0.25];
x_tick       = it_sames;
y_tick       = y_lim(1):0.05:y_lim(end);

rect_w       = 350;
rect_h       = 450;
fig_rect     = [0 0 rect_w*n_condition rect_h*2];

xfig = figure;
set(xfig, 'Position', fig_rect)

for it_bin = 1:length(it_bins)
    xbin = it_bins(it_bin);
    
    for xcond = 1:n_condition
        clear b
        
        subplot(length(it_bins), n_condition, ...
            xcond + n_condition * mod(it_bin, 2))
        
        for xsame = it_sames
            b{xsame} = bar(xsame, xmean{it_bin}(xsame, xcond)); hold on
            set(b{xsame}, 'facecolor', xcolor{xsame});
        end
        
        if xcond==n_condition
            %*************** legend
            xlegend        = same_names(it_sames);
            lg             = legend(xlegend);
            lg.Location    = 'Best';%'SouthEast';
            legend(xlegend,'AutoUpdate','off')
        end
        
        for xsame = it_sames
            errorbar(xsame, xmean{it_bin}(xsame, xcond), xse{it_bin}(xsame, xcond),'k.')
        end
        
        %*************** ttest
        xp = xpvalue{it_bin}(xcond);
        xs = xstats{it_bin}{xcond};
        if xp <= 0.1%xbin_alpha
            yy_sig = y_lim(2) - 0.01;
            
            if (xp <= 0.1) && (xp > xbin_alpha)
                xsig = '+';
            elseif (xp <= xbin_alpha) && (xp > (xbin_alpha/5))
                xsig = '*';
            elseif (xp <= (xbin_alpha/5)) && (xp > (xbin_alpha/50))
                xsig = '**';
            elseif xp <= (xbin_alpha/50)%0.001
                xsig = '***';
            end
            
            plot([1 2],[yy_sig yy_sig]-0.005,'k-');
            text(1.5, yy_sig, xsig, 'FontSize', f_size+10);
            text(1.5, yy_sig-0.015, sprintf('p=%4.4f', xp), 'FontSize', 10);
            text(1.5, yy_sig-0.03, sprintf('T(%s)=%4.4f', num2str(xs.df), xs.tstat), 'FontSize', 10);
        end
        
        %*************** set
        set(gca,'xlim', x_lim,'ylim', y_lim);
        set(gca,'XTick', x_tick, 'YTick', y_tick);
        set(gca,'XTickLabel', x_ticklable, 'FontSize', f_size);
        
        title(sprintf('next item RSA: %s/ %s:%s TR', conds_names{xcond}, ...
            num2str(args.bin_trs{xbin}(1)), num2str(args.bin_trs{xbin}(end))));
        ylabel('RSA');
    end
end

%*************** save fig
fig_fname = fullfile(output_dir, sprintf('plot_bin_next_sameseq_%s_n%s', ...
    basename, num2str(n_subjs)));

savefig(xfig, sprintf('%s.fig', fig_fname));
set(gcf,'PaperUnits','inches','PaperPosition', fig_rect/100)
print('-dpng', '-r300', sprintf('%s.png',fig_fname))
saveas(xfig, sprintf('%s.png',fig_fname), 'png')

close(xfig);

%% *************** baseline (different category) stats
xmtcmp       = 'tukey-kramer';
xconds_names = {'maintain','replace','suppress','clear'};

xlevel       = 1;%1_cate
xbin         = 2;
xtarg        = 1;
xsame        = 1:2;
it_conds     = [1 2 4 5];
y_lim        = [-0.05 0.3];

fig_rect     = [0 0 900 450];

xfig = figure;
set(xfig, 'Position', fig_rect)

for xsame = 1:2
    clear xmean xse xrsa_corr it_matrix xpvalue
    
    subplot(1, 2, xsame)
    
    for it = 1:length(it_conds)
        xcond = it_conds(it);
        %*************** mean/se
        if (xcond == 2) && (xsame == 2)
            xnew_same = 2;%only when new-item and next-item is different
            it_rsa = grp_rsa{xph}.bin_mean.sameseq{xlevel}.cond{xcond}.newitem{xnew_same}.targ{xtarg};
        else
            it_rsa = grp_rsa{xph}.bin_mean.sameseq{xlevel}.cond{xcond}.targ{xtarg};
        end
        
        xmean(it)   = mean(it_rsa.same{xsame}.corr_w(:, xbin));
        xse(it)     = std(it_rsa.same{xsame}.corr_w(:, xbin))/sqrt(n_subjs);
        it_matrix(:,it) = it_rsa.same{xsame}.corr_w(:, xbin);
    end
    
    %*************** anova
    xtable_acc = array2table(it_matrix, 'VariableNames', xconds_names);
    
    xanova     = {'condition'};
    xmeasures  = table((1:length(it_conds))','VariableNames', xanova);
    xrepmeas   = fitrm(xtable_acc,'maintain-clear~1','WithinDesign',xmeasures);
    xanova_out = ranova(xrepmeas);
    
    xdiff_anova = sprintf('one-way ANOVA: F(%s, %s)=%2.4f, p=%1.4f', ...
        num2str(xanova_out.DF(1)), num2str(xanova_out.DF(2)),...
        xanova_out.F(1), xanova_out.pValue(1));
    
    xpvalue_table = multcompare(xrepmeas,'condition','ComparisonType',xmtcmp);%'bonferroni'
    
    %*************** ttest
    xpvalue = nan(length(it_conds)-1, length(it_conds)-1);
    
    for xrow = 1:(length(it_conds)-1)
        for xcol = xrow+1:length(it_conds)
            xunit = (xcol-1) + ((length(it_conds)-1) * (xrow-1));
            xpvalue(xrow, xcol) = xpvalue_table.pValue(xunit);
        end
    end
    
    %*************** plot
    clear b
    for it = 1:length(it_conds)
        xcond = it_conds(it);
        b{it} = bar(it, xmean(it)); hold on
        set(b{it}, 'facecolor', xcond_color{xcond});
    end
    
    if xcond==n_condition
        %*************** legend
        xlegend        = xconds_names;
        lg             = legend(xlegend);
        lg.Location    = 'SouthEast';%'SouthEast';
        legend(xlegend,'AutoUpdate','off')
    end
    
    for it = 1:length(it_conds)
        errorbar(it, xmean(it), xse(it),'k.')
    end
    
    %*************** ttest
    n = 1;
    for xrow = 1:(length(it_conds)-1)
        for xcol = xrow+1:length(it_conds)
            xp = xpvalue(xrow, xcol);
            if xp <= 0.1%xbin_alpha
                yy_sig = y_lim(2) - (0.015 * n);
                n = n+1;
                if (xp <= 0.1) && (xp > xbin_alpha)
                    xsig = '+';
                elseif (xp <= xbin_alpha) && (xp > (xbin_alpha/5))
                    xsig = '*';
                elseif (xp <= (xbin_alpha/5)) && (xp > (xbin_alpha/50))
                    xsig = '**';
                elseif xp <= (xbin_alpha/50)%0.001
                    xsig = '***';
                end
                
                plot([xrow xcol],[yy_sig yy_sig],'k-');
                text(xrow + (xcol-xrow)/2, yy_sig, xsig, 'FontSize', f_size+10);
                text(xcol, yy_sig, sprintf('p=%4.4f', xp), 'FontSize', 10);
            end
        end
    end
    
    %*************** set
    set(gca,'xlim', [0 length(it_conds)+1],'ylim', y_lim);
    set(gca,'XTick', 1:length(it_conds), 'YTick', y_lim(1):0.05:y_lim(2));
    set(gca,'XTickLabel', xconds_names, 'FontSize', f_size);
    
    title(sprintf('next item RSA: %s', same_names{xsame}));
    ylabel('RSA');
end

%*************** save fig
fig_fname = fullfile(output_dir, sprintf('plot_bin_cond_next_sameseq_%s_n%s', ...
    basename, num2str(n_subjs)));

savefig(xfig, sprintf('%s.fig', fig_fname));
set(gcf,'PaperUnits','inches','PaperPosition', fig_rect/100)
print('-dpng', '-r300', sprintf('%s.png',fig_fname))
saveas(xfig, sprintf('%s.png',fig_fname), 'png')

close(xfig);

%% ============= NEXT-TRIAL BIN: SAME-DIFFERENCE: TARGET-ONLY
% for replace(category): when n~=n+1, using only new-item~=n+1 item
clear xmean xse xpvalue xdiff xstats xrsa it_rsa
xlevel = 1;
xtarg  = 1;
xdiff_table = zeros(((n_condition-1) * n_subjs), 3);

for it_bin = 1:length(it_bins)
    xbin = it_bins(it_bin);
    
    for xcond = 1:n_condition
        %*************** mean/se
        for xsame = it_sames
            
            if (xcond == 2) && (xsame == 2)
                xnew_same = 2;%only when new-item and next-item is different
                it_rsa = grp_rsa{xph}.bin_mean.sameseq{xlevel}.cond{xcond}.newitem{xnew_same}.targ{xtarg};
            else
                it_rsa = grp_rsa{xph}.bin_mean.sameseq{xlevel}.cond{xcond}.targ{xtarg};
            end
            
            xrsa.same{xsame}.corr_w = it_rsa.same{xsame}.corr_w;
            
        end
        
        xdiff{it_bin}{xcond} = xrsa.same{1}.corr_w(:, xbin) - xrsa.same{2}.corr_w(:, xbin);
        
        %*************** mean/se
        xmean{it_bin}(xcond) = mean(xdiff{it_bin}{xcond});
        xse{it_bin}(xcond)   = std(xdiff{it_bin}{xcond})/sqrt(n_subjs);
        
        %*************** table
        if it_bin == 1
            if xcond~=3
                xit = find(ismember(it_conds, xcond));
                xunit = (1:n_subjs) + (n_subjs * (xit-1));
                xdiff_table(xunit, 1) = 1:n_subjs;
                xdiff_table(xunit, 2) = xcond;
                xdiff_table(xunit, 3) = xdiff{it_bin}{xcond};
            end
        end
    end
    
    %*************** ttest
    xpvalue{it_bin} = nan(n_condition-1, n_condition-1);
    
    for xrow = 1:(n_condition-1)
        for xcol = xrow+1:n_condition
            
            [~, xpvalue{it_bin}(xrow, xcol)] = ttest(xdiff{it_bin}{xrow}, xdiff{it_bin}{xcol}, 'Alpha', xbin_alpha);
        end
    end
end

%% *************** save table for R violin plot
table_header = {'subject','condition','nextencoding'};

xtable = array2table(xdiff_table, 'VariableNames', table_header);
writetable(xtable, fullfile(dirs.rsa.group.pattern, 'group_rsa_next_item.csv'));

%% *************** figure
x_ticklable  = conds_names;
x_lim        = [1 n_condition]+[-0.5 0.5];
y_lim        = [-0.06 0.07];
x_tick       = 1:n_condition;
y_tick       = y_lim(1):0.01:y_lim(end);

fig_rect     = [0 0 700 1200];

xfig = figure;
set(xfig, 'Position', fig_rect)

for it_bin = 1:length(it_bins)
    xbin = it_bins(it_bin);
    
    subplot(length(it_bins), 1, mod(it_bin, 2)+1)
    
    clear b
    for xcond = 1:n_condition
        b{xcond} = bar(xcond, xmean{it_bin}(xcond)); hold on
        set(b{xcond}, 'facecolor', xcond_color{xcond});
    end
    
    %*************** legend
    xlegend        = conds_names;
    lg             = legend(xlegend);
    lg.Location    = 'BestOutside';%'SouthEast';
    legend(xlegend,'AutoUpdate','off')
    lg.FontSize    = f_size-2;
    
    %*************** errorbar
    errorbar(1:n_condition, xmean{it_bin}, xse{it_bin},'k.')
    
    %*************** ttest
    n = 1;
    for xrow = 1:(n_condition-1)
        for xcol = xrow+1:n_condition
            
            xp = xpvalue{it_bin}(xrow, xcol);
            
            if xp <= xbin_alpha
                yy_sig = y_lim(2) - (0.005 * n);
                n = n+1;
                
                if (xp <= xbin_alpha) && (xp > (xbin_alpha/5))
                    xsig = '*';
                elseif (xp <= (xbin_alpha/5)) && (xp > (xbin_alpha/50))
                    xsig = '**';
                elseif xp <= (xbin_alpha/50)%0.001
                    xsig = '***';
                end
                
                plot([xrow xcol],[yy_sig yy_sig]-0.001,'k-');
                text(xrow + 0.5, yy_sig, xsig, 'FontSize', f_size+5);
                text(xrow + 0.8, yy_sig+0.0015, sprintf('p=%4.4f', xp), 'FontSize', f_size);
            end
        end
    end
    
    %*************** set
    set(gca,'xlim', x_lim,'ylim', y_lim);
    set(gca,'XTick', x_tick, 'YTick', y_tick);
    set(gca,'XTickLabel', x_ticklable, 'FontSize', f_size);
    xtickangle(45)
    
    title(sprintf('item N+1: same-differ/ %s:%s TR', num2str(args.bin_trs{xbin}(1)), num2str(args.bin_trs{xbin}(end))));
    ylabel('RSA');
    
end

%*************** save fig
fig_fname = fullfile(output_dir, sprintf('plot_bin_next_diff_sameseq_%s_n%s', ...
    basename, num2str(n_subjs)));

savefig(xfig, sprintf('%s.fig', fig_fname));
set(gcf,'PaperUnits','inches','PaperPosition', fig_rect/100)
print('-dpng', '-r300', sprintf('%s.png',fig_fname))
saveas(xfig, sprintf('%s.png',fig_fname), 'png')

close(xfig);

%% *************** SAME-DIFFERENCE ANOVA: [1 2 4 5] CONDS
clear it_matrix xpvalue

it_bin = 1;
xbin   = it_bins(it_bin);%N+1 peak

it_conds     = [1 2 4 5];
xmtcmp       = 'tukey-kramer';
xconds_names = {'maintain','replace_category','suppress','clear'};

for it = 1:length(it_conds)
    xcond = it_conds(it);
    it_matrix(:, it) = xdiff{it_bin}{xcond};
end

xtable_acc = array2table(it_matrix, 'VariableNames', xconds_names);

xanova     = {'condition'};
xmeasures  = table((1:length(it_conds))','VariableNames', xanova);
xrepmeas   = fitrm(xtable_acc,'maintain-clear~1','WithinDesign',xmeasures);
xanova_out = ranova(xrepmeas);

xdiff_anova = sprintf('one-way ANOVA: F(%s, %s)=%2.4f, p=%1.4f', ...
    num2str(xanova_out.DF(1)), num2str(xanova_out.DF(2)),...
    xanova_out.F(1), xanova_out.pValue(1));

xpvalue_table = multcompare(xrepmeas,'condition','ComparisonType',xmtcmp);%'bonferroni'

%*************** ttest
xpvalue = nan(length(it_conds)-1, length(it_conds)-1);

for xrow = 1:(length(it_conds)-1)
    for xcol = xrow+1:length(it_conds)
        xunit = (xcol-1) + ((length(it_conds)-1) * (xrow-1));
        xpvalue(xrow, xcol) = xpvalue_table.pValue(xunit);
    end
end

%% *************** figure
x_ticklable  = conds_names(it_conds);
x_lim        = [1 length(it_conds)]+[-0.5 0.5];
y_lim        = [-0.06 0.07];
x_tick       = 1:length(it_conds);
y_tick       = y_lim(1):0.01:y_lim(end);

fig_rect     = [0 0 700 600];

xfig = figure;
set(xfig, 'Position', fig_rect)

clear b
for it = 1:length(it_conds)
    xcond = it_conds(it);
    b{it} = bar(it, xmean{it_bin}(xcond)); hold on
    set(b{it}, 'facecolor', xcond_color{xcond});
    
    text(it-0.2, 0.01, sprintf('M=%1.3f\nSE=%1.3f', ...
        xmean{it_bin}(xcond), xse{it_bin}(xcond)), 'FontSize', f_size);
end

%*************** legended
xlegend        = {'maintain','replace (category)','suppress','clear'};%conds_names(it_conds);
lg             = legend('maintain','replace (category)','suppress','clear');
lg.Location    = 'BestOutside';%'SouthEast';
% legend(xlegend,'AutoUpdate','off')
lg.FontSize    = f_size-2;

%*************** errorbar
errorbar(1:length(it_conds), xmean{it_bin}(it_conds), xse{it_bin}(it_conds),'k.')

%*************** ttest
n = 1;
for xrow = 1:(length(it_conds)-1)
    for xcol = xrow+1:length(it_conds)
        
        xp = xpvalue(xrow, xcol);
        
        if xp <= xbin_alpha
            if (xp <= xbin_alpha) && (xp > (xbin_alpha/5))
                xsig = '*';
            elseif (xp <= (xbin_alpha/5)) && (xp > (xbin_alpha/50))
                xsig = '**';
            elseif xp <= (xbin_alpha/50)%0.001
                xsig = '***';
            end
        else
            xsig = '-';
        end
        n = n+1;
        yy_sig = y_lim(2) - (0.005 * n);
        
        plot([xrow xcol],[yy_sig yy_sig]-0.001,'k-');
        text(xrow + 0.5, yy_sig, xsig, 'FontSize', f_size+5);
        text(xrow + 0.8, yy_sig+0.0015, sprintf('p=%4.4f', xp), 'FontSize', f_size);
    end
end

text(2, y_lim(1)+0.01, xdiff_anova);

%*************** set
set(gca,'xlim', x_lim,'ylim', y_lim);
set(gca,'XTick', x_tick, 'YTick', y_tick);
set(gca,'XTickLabel', x_ticklable, 'FontSize', f_size);
xtickangle(45)

title(sprintf('item N+1 RSA: same-differ/ %s, %s:%s TR', xmtcmp, num2str(args.bin_trs{xbin}(1)), num2str(args.bin_trs{xbin}(end))));
ylabel('RSA');

%*************** save fig
fig_fname = fullfile(output_dir, sprintf('plot_bin_next_diff_4conds_sameseq_%s_n%s', ...
    basename, num2str(n_subjs)));

savefig(xfig, sprintf('%s.fig', fig_fname));
set(gcf,'PaperUnits','inches','PaperPosition', fig_rect/100)
print('-dpng', '-r300', sprintf('%s.png',fig_fname))
saveas(xfig, sprintf('%s.png',fig_fname), 'png')

% close(xfig);

%% ============= REPLACE N+1 SAME/DIFF BASED ON NEW-ITEM

xph    = 2;
xlevel = 1;%1_cate
xcond  = 2;%2_repCat
xsame  = 2;%2_diff
xtarg  = 1;

clear xx xmean xse xpvalue
for it = 1:length(it_bins)
    xbin = it_bins(it);
    for xnew_same = 1:2
        xx{it}{xnew_same} = grp_rsa{xph}.bin_mean.sameseq{xlevel}.cond{xcond}.newitem{xnew_same}.targ{xtarg}.same{xsame}.corr_w(:,xbin);
        
        xmean{it}(xnew_same) = mean(xx{it}{xnew_same});
        xse{it}(xnew_same)   = std(xx{it}{xnew_same})/sqrt(n_subjs);
    end
    
    %*************** ttest
    [~, xpvalue{it}, ~, xstats{it}] = ttest(xx{it}{1}, xx{it}{2}, 'Alpha', xbin_alpha);
    
    xx{it}{3} = grp_rsa{xph}.bin_mean.sameseq{xlevel}.cond{xcond}.newitem{1}.targ{xtarg}.same{xsame}.corr_w(:,xbin) - ...
        grp_rsa{xph}.bin_mean.sameseq{xlevel}.cond{xcond}.newitem{2}.targ{xtarg}.same{xsame}.corr_w(:, xbin);
    
    xmean{it}(3) = mean(xx{it}{3});
    xse{it}(3)   = std(xx{it}{3})/sqrt(n_subjs);
    
end

%*************** ttest
[~, xpvalue{3}, ~, xstats{3}] = ttest(xx{1}{3}, xx{2}{3}, 'Alpha', xbin_alpha);

%% *************** plot
it_bins  = [2 3];
x_lim    = [0.5 3.5];
y_lim    = [-0.05 0.4];
y_tick   = y_lim(1):0.05:y_lim(end);
fig_rect = [0 0 1000 500];

xfig = figure;
set(xfig, 'Position', fig_rect)

for it = 1:length(it_bins)
    xbin = it_bins(it);
    subplot(1, 3, it)
    
    for xnew_same = 1:3
        b{xnew_same} = bar(xnew_same, xmean{it}(xnew_same)); hold on
    end
    
    for xnew_same = 1:2
        set(b{xnew_same}, 'facecolor', xcolor{xnew_same});
    end
    set(b{3}, 'facecolor', 'w');
    
    %*************** legend
    xlegend        = same_names(it_sames);
    xlegend{3}     = 'same-diff';
    lg             = legend(xlegend);
    lg.Location    = 'NorthEast';%'SouthEast';
    legend(xlegend,'AutoUpdate','off')
    
    for xnew_same = 1:3
        errorbar(xnew_same, xmean{it}(xnew_same), xse{it}(xnew_same),'k.')
    end
    
    %*************** ttest
    xp = xpvalue{it}; xs = xstats{it};
    if xp <= 0.1%xbin_alpha
        yy_sig = 0.3 - 0.01;
        
        if (xp <= 0.1) && (xp > xbin_alpha)
            xsig = '+';
        elseif (xp <= xbin_alpha) && (xp > (xbin_alpha/5))
            xsig = '*';
        elseif (xp <= (xbin_alpha/5)) && (xp > (xbin_alpha/50))
            xsig = '**';
        elseif xp <= (xbin_alpha/50)%0.001
            xsig = '***';
        end
        
        plot([1 2],[yy_sig yy_sig]-0.005,'k-');
        text(1.5, yy_sig, xsig, 'FontSize', f_size+10);
        text(1.5, yy_sig-0.02, ...
            sprintf('T(%s)=%4.4f\np=%4.4f', num2str(xs.df), xs.tstat, xp), 'FontSize', 10);
    end
    
    %*************** set
    set(gca,'xlim', x_lim,'ylim', y_lim);
    
    xtrs = args.bin_trs{xbin};
    title(sprintf('%s, %s-%s trs', ...
        conds_names{xcond}, num2str(xtrs(1)), num2str(xtrs(end))));
    ylabel('RSA');
end

subplot(1, 3, 3)

for xnew_same = 1:2
    b{xnew_same} = bar(xnew_same, mean(n_newsame{xnew_same})); hold on
    set(b{xnew_same}, 'facecolor', xcolor{xnew_same});
    text(xnew_same, 2, sprintf('n=%4.4f', mean(n_newsame{xnew_same})));
end

%*************** legend
xlegend        = same_names(it_sames);
lg             = legend(xlegend);
lg.Location    = 'NorthEast';%'SouthEast';
legend(xlegend,'AutoUpdate','off')

for xnew_same = 1:2
    errorbar(xnew_same, mean(n_newsame{xnew_same}), ...
        std(n_newsame{xnew_same})/sqrt(n_subjs),'k.')
end

%*************** set
set(gca,'xlim', [0.5 2.5],'ylim', [0 30]);

title(sprintf('n trials: %s', conds_names{xcond}));
ylabel('n trials');

%*************** save fig
fig_fname = fullfile(output_dir, sprintf('plot_repCat_next_sameseq_%s_n%s', ...
    basename, num2str(n_subjs)));

savefig(xfig, sprintf('%s.fig', fig_fname));
set(gcf,'PaperUnits','inches','PaperPosition', fig_rect/100)
print('-dpng', '-r300', sprintf('%s.png',fig_fname))
saveas(xfig, sprintf('%s.png',fig_fname), 'png')

% close(xfig);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ============= RANDOM: NEXT-TRIAL BIN: TARGET-ONLY
% correlation b/w n-trial & n+1 RSA
clear xmean xse xpvalue

xlevel     = 1;
xbin       = 2;
xtarg      = 1;
xbin_alpha = 0.05;

if xlevel==1, it_sames = 1:2; else, it_sames = 1:3; end

for xcond = 1:n_condition
    clear xrsa
    xrsa = grp_rsa{xph}.bin_trial.sameseq{xlevel}.cond{xcond}.targ{xtarg};
    
    %*************** ttest
    [~, xpvalue(xcond)] = ttest(xrsa.same{1}.r, xrsa.same{2}.r, 'Alpha', xbin_alpha);
    
    %*************** mean/se
    for xsame = it_sames
        xmean(xsame, xcond) = mean(xrsa.same{xsame}.r);
        xse(xsame, xcond)   = std(xrsa.same{xsame}.r)/sqrt(n_subjs);
    end
end

%% *************** figure
f_size       = 10;
x_ticklable  = same_names(it_sames);
x_lim        = [1 length(it_sames)]+[-0.5 0.5];
y_lim        = [-0.1 0.3];
x_tick       = it_sames;
y_tick       = y_lim(1):0.05:y_lim(end);

rect_w       = 350;
rect_h       = 450;
fig_rect     = [0 0 rect_w*n_condition rect_h];

xfig = figure;
set(xfig, 'Position', fig_rect)

for xcond = 1:n_condition
    clear b
    
    subplot(1,n_condition, xcond)
    
    for xsame = it_sames
        b{xsame} = bar(xsame, xmean(xsame, xcond)); hold on
        set(b{xsame}, 'facecolor', xcolor{xsame});
    end
    
    if xcond==n_condition
        %*************** legend
        xlegend        = same_names(it_sames);
        lg             = legend(xlegend);
        lg.Location    = 'Best';%'SouthEast';
        legend(xlegend,'AutoUpdate','off')
    end
    
    for xsame = it_sames
        errorbar(xsame, xmean(xsame, xcond), xse(xsame, xcond),'k.')
    end
    
    %*************** ttest
    if xpvalue(xcond) <= 0.1%xbin_alpha
        yy_sig = y_lim(2) - 0.01;
        
        if (xpvalue(xcond) <= 0.1) && (xpvalue(xcond) > xbin_alpha)
            xsig = '+';
        elseif (xpvalue(xcond) <= xbin_alpha) && (xpvalue(xcond) > (xbin_alpha/5))
            xsig = '*';
        elseif (xpvalue(xcond) <= (xbin_alpha/5)) && (xpvalue(xcond) > (xbin_alpha/50))
            xsig = '**';
        elseif xpvalue(xcond) <= (xbin_alpha/50)%0.001
            xsig = '***';
        end
        
        plot([1 2],[yy_sig yy_sig]-0.005,'k-');
        text(1.5, yy_sig, xsig, 'FontSize', f_size+10);
        text(1.5, yy_sig-0.015, sprintf('p=%4.4f', xpvalue(xcond)), 'FontSize', 10);
    end
    
    %*************** set
    set(gca,'xlim', x_lim,'ylim', y_lim);
    set(gca,'XTick', x_tick, 'YTick', y_tick);
    set(gca,'XTickLabel', x_ticklable, 'FontSize', f_size);
    
    title(sprintf('N & N+1: %s', conds_names{xcond}));
    ylabel('R (within subject)');
end

%*************** save fig
fig_fname = fullfile(output_dir, sprintf('plot_bin_sameseq_corr_bar_rand_%s_n%s', ...
    basename, num2str(n_subjs)));

savefig(xfig, sprintf('%s.fig', fig_fname));
set(gcf,'PaperUnits','inches','PaperPosition', fig_rect/100)
print('-dpng', '-r300', sprintf('%s.png',fig_fname))
saveas(xfig, sprintf('%s.png',fig_fname), 'png')

% close(xfig);

%% ============= FIXED: CORR: CURRENT VS. NEXT-TRIAL BIN: TARGET-ONLY
%% *************** bootstrap
n_iteration = 1000;
xtarg = 1;

for xboot = 1:n_iteration
    
    if ~(mod(xboot,100))
        fprintf('%s.',num2str(xboot));
    end
    
    %************** SAMPLING
    xsample = sort(datasample(xsubj_grp, n_subjs, 'Replace', true));
    
    for xcond = 1:n_condition
        for xsame = it_sames
            clear xrsa xx yy xr
            xx = []; yy = [];
            for xsub = xsample
                xrsa = grp_rsa{xph}.sub{xsub}.bin_trial.sameseq{xlevel}.cond{xcond}.targ{xtarg};
                
                %*************** corr
                xx = vertcat(xx, xrsa.same{xsame}.corr_w(:, 1));
                yy = vertcat(yy, xrsa.same{xsame}.corr_w(:, 2));
            end
            
            xr = corrcoef(xx, yy);
            
            xboot_corr{xcond}{xsame}(xboot) = xr(1,2);
        end
    end
end

fprintf('\n');

%% *************** plot
x_lim        = [-1 1.5];
y_lim        = [-0.05 0.4];
y_tick       = y_lim(1):0.05:y_lim(end);

xfig = figure;
set(xfig, 'Position', fig_rect)

for xcond = 1:n_condition
    clear xrsa
    
    subplot(1, n_condition, xcond); hold on;
    
    %*************** legend
    for xsame = it_sames
        plot([0 0], [0 0], 'Color', xcolor{xsame}, 'LineWidth', 2)
        hold on
    end
    
    xlegend        = same_names(it_sames);
    lg             = legend(xlegend);
    lg.Location    = 'SouthEast';%'SouthEast';
    legend(xlegend,'AutoUpdate','off')
    
    xrsa = grp_rsa{xph}.bin_trial.sameseq{xlevel}.cond{xcond}.targ{xtarg};
    
    for xsame = it_sames
        clear xx yy xr xpvalue
        
        %*************** corr
        xx = xrsa.same{xsame}.corr_w(:, 1);
        yy = xrsa.same{xsame}.corr_w(:, 2);
        
        [xr, xpvalue] = corrcoef(xx, yy);
        
        %*************** plot
        P = polyfit(xx, yy, 1);%linear fit
        y_fit = P(1) * xx + P(2);
        
        %         scatter(xx, yy, 20, xcolor{xsame});
        plot(xx, y_fit, 'Color', xcolor{xsame}, 'LineWidth', 2)
        
        text(x_lim(1)+0.05, y_lim(2) - (0.02 * xsame), ...
            sprintf('%s: r2 = %1.4f, p = %1.4f', same_names{xsame}, xr(1,2)^2, xpvalue(1,2)),...
            'FontSize', f_size, 'FontWeight', 'bold');
    end
    
    %*************** bootstrap
    xboot_pvalue = 1-size(find(xboot_corr{xcond}{1} > xboot_corr{xcond}{2}), 2)/n_iteration;
    text(x_lim(1)+0.05, y_lim(2) - 0.1, ...
        sprintf('bootstrap(1000)\nsame>differ: p = %1.4f', xboot_pvalue),...
        'FontSize', f_size, 'FontWeight', 'bold');
    
    %*************** set
    set(gca,'xlim', x_lim, 'ylim', y_lim);
    set(gca,'YTick', y_tick)
    
    title(conds_names{xcond});
    xlabel('N tiral item'); ylabel('N+1 trial item');
    
end

%*************** save fig
fig_fname = fullfile(output_dir, sprintf('plot_bin_sameseq_corr_fixed_%s_n%s', ...
    basename, num2str(n_subjs)));

savefig(xfig, sprintf('%s.fig', fig_fname));
set(gcf,'PaperUnits','inches','PaperPosition', fig_rect/100)
print('-dpng', '-r300', sprintf('%s.png',fig_fname))
saveas(xfig, sprintf('%s.png',fig_fname), 'png')

close(xfig);

grp_rsa{xph}.sub{xsub}.bin_corr.sameseq{xlevel}.boot_corr = xboot_corr;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%============= STUDY PLOT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ============= TIMECOURSE: PER TARGET
f_size       = 8;
x_tick       = 1:n_trs;
x_ticklable  = x_tick;
x_lim        = [x_tick(1) x_tick(end)];
y_lim        = [-0.3 0.3];
y_tick       = y_lim(1):0.05:y_lim(end);

rect_w       = 350;
rect_h       = 450;
fig_rect     = [0 0 rect_w*2 rect_h*n_condition];
n_targs       = 3;

xfig = figure;
set(xfig, 'Position', fig_rect)

for xcond = 1:n_condition
    it_subplt = 1:2;
    
    subplot(n_condition, 1, xcond)
    
    %*************** mean line plots
    clear xmean xse fity xpvalue
    
    fitx = linspace(1, n_trs, n_trs*10);
    
    for xtarg = 1:n_targs
        xmean{xtarg} = grp_rsa{xph}.timecourse.cond{xcond}.mean_w(xtarg, x_tick);
        xse{xtarg}   = grp_rsa{xph}.timecourse.cond{xcond}.se_w(xtarg, x_tick);
        
        fity{xtarg}  = interp1(x_tick, xmean{xtarg}, fitx,'spline');
        
        plot(fitx, fity{xtarg}, '-','Color', xtarg_color{xtarg}, 'LineWidth', 2); hold on;
    end
    
    %*************** legend
    xlegend = {'target','nontarget','baseline'};
    lg             = legend(xlegend);
    lg.Location    = 'southeast';%'bestoutside';
    lg.FontSize    = f_size;
    legend(xlegend,'AutoUpdate','off')
    grid on
    
    %*************** std error-bar filling
    for xtarg = 1:n_targs
        clear xerr fit_err
        xerr(1,:) = xmean{xtarg} - xse{xtarg};
        xerr(2,:) = xmean{xtarg} + xse{xtarg};
        
        for i = 1:2, fit_err(:,i) = interp1(x_tick, xerr(i,:), fitx,'spline'); end
        
        in_between = [fit_err(:,1)', fliplr(fit_err(:,2)')];
        h = fill([fitx fliplr(fitx)], in_between, xtarg_color{xtarg});
        set(h,'facealpha', .4)
        
        plot(fitx, fit_err(:,1), '-', 'Color', xtarg_color{xtarg})
        plot(fitx, fit_err(:,2), '-', 'Color', xtarg_color{xtarg})
    end
    
    %*************** stats stars
    n = 0;
    for xcol = 1:(n_targs-1)
        for xrow = (xcol+1):n_targs
            yy_sig = (y_lim(2)-0.02) - (n * 0.02);
            text(-3.5, yy_sig, sprintf('%s vs. %s', xlegend{xcol}, ...
                xlegend{xrow}), 'FontSize', f_size/1.5);
            
            for xblk = 1:n_stat_blks
                
                xpvalue = grp_rsa{xph}.timecourse.cond{xcond}.blk{xblk}.pvalue_w(xrow, xcol);
                xx_sig  = (n_tr_blks/2) + (n_tr_blks * (xblk-1)) + 0.5 - 0.08;
                
                if xpvalue <= xalpha
                    text(xx_sig, yy_sig, '*', 'FontSize', f_size*1.5);
                end
                
                xx = n_tr_blks + (n_tr_blks * (xblk-1)) + 0.5;
                plot([xx xx], y_lim, '--', 'Color', 'b')
                
            end
            
            n = n + 1;
        end
    end
    
    %*************** real onset lines
    plot([dur_stim dur_stim]+1, y_lim, '-', 'Color', xonset_color)
    plot([(dur_stim + dur_manip) (dur_stim + dur_manip)]+1, y_lim, '-', 'Color', xonset_color)
    
    h1 = text(1.5, y_lim(1)+0.05, '(1) stim onset', 'Color', xonset_color, 'FontSize', f_size);
    h2 = text(dur_stim + 1.5, y_lim(1)+0.01, '(2) operation onset', 'Color', xonset_color, 'FontSize', f_size);
    h3 = text(dur_stim + dur_manip + 1.5, y_lim(1)+0.05, '(3) fixation onset', 'Color', xonset_color, 'FontSize', f_size);
    set(h1,'Rotation', 90); set(h2,'Rotation', 90); set(h3,'Rotation', 90);
    
    %*************** set gca
    set(gca,'xlim', x_lim, 'ylim', y_lim);
    set(gca,'XTick', x_tick,'YTick', y_tick)
    set(gca,'XTickLabel', x_ticklable);
    set(gca,'FontSize', f_size/1.5)
    
    title(sprintf('RSAw item : %s (N=%s, p<%1.4f)', ...
        conds_names{xcond}, num2str(n_subjs), xalpha),...
        'FontSize', f_size*1.5,'FontWeight','bold');
    xlabel('Volume (tr)','FontSize', f_size);
    ylabel('Similarity (Z Pearson R)','FontSize', f_size);
    
end

%*************** save fig
fig_fname = fullfile(output_dir, sprintf('plot_tc_blk%s_%s_%s_%s_n%s', ...
    num2str(n_tr_blks), args.level, args.phase_name{xph}, basename, num2str(n_subjs)));

savefig(xfig, sprintf('%s.fig', fig_fname));
set(gcf,'PaperUnits','inches','PaperPosition', fig_rect/100)
print('-dpng', '-r300', sprintf('%s.png',fig_fname))
saveas(xfig, sprintf('%s.png',fig_fname), 'png')

close(xfig);

%% ============= TIMECOURSE: TARGET-Only
xfill           = 1;
rect_w          = 250;

cat_x_tick      = horzcat(1:n_trs, (1:dur_sync) + n_trs);
cat_x_lim       = [1 (dur_sync + n_trs)];
cat_x_ticklable = horzcat(1:n_trs, (1:dur_sync) + 21);

sel_conds{1}  = 1:n_condition;
sel_conds{2}  = [1 2 3];
sel_conds{3}  = [1 4];
sel_conds{4}  = [1 2];
sel_conds{5}  = [2 3];
sel_conds{6}  = [2 4];
sel_conds{7}  = [2 5];
sel_conds{8}  = [4 5];

n_sels          = length(sel_conds);
fig_rect        = [0 0 rect_w*3 rect_h*n_sels];

xfig = figure;
set(xfig, 'Position', fig_rect)

%*************** plot selected conditions
for xsel = 1:length(sel_conds)
    
    xsel_conds = sel_conds{xsel};
    
    subplot(n_sels, 1, xsel)
    
    %*************** time window
    for xwin = 1:2 %1_tr, 2_sync
        
        if xwin==1
            it_trs       = n_trs;
            x_tick       = 1:it_trs;
            fitx         = linspace(1, it_trs, it_trs*10);
        else
            it_trs       = dur_sync;
            x_tick       = (1:dur_sync) + n_trs;
            fitx         = linspace(n_trs + 1, n_trs + dur_sync + 1, dur_sync*10);
        end
        
        %*************** mean line plots
        clear fity xmean xstd
        
        for xcond = xsel_conds
            if xwin == 1
                xmean{xcond} = grp_rsa{xph}.timecourse.cond{xcond}.mean_w(1, 1:it_trs);
                xse{xcond}   = grp_rsa{xph}.timecourse.cond{xcond}.se_w(1, 1:it_trs);
            elseif xwin == 2
                xmean{xcond} = grp_rsa{xph}.timecourse.sync.cond{xcond}.mean_w(1, 1:it_trs);
                xse{xcond}   = grp_rsa{xph}.timecourse.sync.cond{xcond}.se_w(1, 1:it_trs);
            end
            
            fity{xcond}  = interp1(x_tick, xmean{xcond}, fitx,'spline');
            
            plot(fitx, fity{xcond}, '-','Color', xcond_color{xcond}, 'LineWidth', 2); hold on;
            
        end
        
        %*************** legend
        clear xlegend
        for xcond = 1:length(xsel_conds)
            it_cond = xsel_conds(xcond);
            xlegend{xcond} = conds_names{it_cond};
        end
        lg             = legend(xlegend);
        lg.Location    = 'southeast';%'bestoutside';'best';%
        lg.FontSize    = f_size/1.5;
        legend(xlegend,'AutoUpdate','off')
        grid on
        
        if xfill
            %*************** std error-bar filling
            for xcond = xsel_conds %#ok<*UNRCH>
                clear xerr fit_err
                xerr(1,:) = xmean{xcond} - xse{xcond};
                xerr(2,:) = xmean{xcond} + xse{xcond};
                
                for i = 1:2, fit_err(:,i) = interp1(x_tick, xerr(i,:), fitx,'spline'); end
                
                in_between = [fit_err(:,1)', fliplr(fit_err(:,2)')];
                h = fill([fitx fliplr(fitx)], in_between, xcond_color{xcond});
                set(h,'facealpha', .4)
                
                plot(fitx, fit_err(:,1), '-', 'Color', xcond_color{xcond})
                plot(fitx, fit_err(:,2), '-', 'Color', xcond_color{xcond})
            end
        end
        
        %*************** stats stars
        n = 0;
        
        for i = 1:(length(xsel_conds)-1)
            xcol = xsel_conds(i);
            for j = i+1:length(xsel_conds)
                xrow = xsel_conds(j);
                
                yy_sig = (y_lim(2)-0.02) - (n * 0.02);
                
                if xsel==1
                    text(-2, yy_sig, sprintf('%d vs. %d', xcol, xrow), 'FontSize', f_size);
                end
                
                for xblk = 1:it_trs/n_tr_blks
                    
                    if xwin == 1
                        xpvalue = grp_rsa{xph}.timecourse.targs.blk{xblk}.pvalue_w(xrow, xcol);
                        xx_blk = xblk;
                    else
                        xpvalue = grp_rsa{xph}.timecourse.sync.targs.blk{xblk}.pvalue_w(xrow, xcol);
                        xx_blk = xblk + (n_trs/n_tr_blks);
                    end
                    
                    xx_sig = (n_tr_blks/2) + (n_tr_blks * (xx_blk-1)) + 0.5 - 0.08;
                    
                    if xpvalue <= xalpha
                        text(xx_sig, yy_sig, '*', 'FontSize', f_size*1.5);
                    end
                    
                    xx = n_tr_blks + (n_tr_blks * (xx_blk-1)) + 0.5;
                    plot([xx xx], y_lim, '--', 'Color', 'b')
                    
                end
                
                n = n + 1;
            end
        end
    end
    
    %*************** real onset lines
    plot([dur_stim dur_stim]+1, y_lim, '-', 'Color', xonset_color)
    plot([(dur_stim + dur_manip) (dur_stim + dur_manip)]+1, y_lim, '-', 'Color', xonset_color)
    plot([(dur_stim + dur_manip + dur_fix) (dur_stim + dur_manip + dur_fix)]+1, y_lim, '--', 'Color', xonset_color, 'Linewidth', 2)
    plot([n_trs n_trs] + 1, y_lim, '-', 'Color', 'k', 'Linewidth', 3)
    
    h1 = text(1.5, y_lim(1)+0.05, '(1) stim onset', 'Color', xonset_color, 'FontSize', f_size);
    h2 = text(dur_stim + 1.5, y_lim(1)+0.05, '(2) operation onset', 'Color', xonset_color, 'FontSize', f_size);
    h3 = text(dur_stim + dur_manip + 1.5, y_lim(1)+0.05, '(3) fixation onset', 'Color', xonset_color, 'FontSize', f_size);
    h4 = text(dur_stim + dur_manip + dur_fix + 1.5, y_lim(1)+0.05, 'synch point', 'Color', 'k', 'FontSize', f_size);
    h5 = text(n_trs + 1.5, y_lim(1)+0.05, 'synch on next-trial onset', 'Color', 'k', 'FontSize', f_size);
    set(h1,'Rotation', 90); set(h2,'Rotation', 90); set(h3,'Rotation', 90);
    set(h4,'Rotation', 90); set(h5,'Rotation', 90);
    
    %*************** set gca
    set(gca,'xlim', cat_x_lim, 'ylim', y_lim);
    set(gca,'XTick', cat_x_tick,'YTick', y_tick)
    set(gca,'XTickLabel', cat_x_ticklable);
    set(gca,'FontSize', f_size/1.5)
    
    title(sprintf('RSAw item targets (N=%s, p<%1.4f)', ...
        num2str(n_subjs), xalpha),...
        'FontSize',f_size*1.5,'FontWeight','bold');
    xlabel('Volume (tr)', 'FontSize', f_size);
    ylabel('Similarity (Pearson R)', 'FontSize', f_size);
end

%*************** save fig
fig_fname = fullfile(output_dir, sprintf('plot_tc_targs_blk%s_%s_%s_fill%d_%s_n%s', ...
    num2str(n_tr_blks), args.level, args.phase_name{xph}, xfill, basename, num2str(n_subjs)));

savefig(xfig, sprintf('%s.fig', fig_fname));
set(gcf,'PaperUnits','inches','PaperPosition', fig_rect/100)
print('-dpng', '-r300', sprintf('%s.png',fig_fname))
saveas(xfig, sprintf('%s.png',fig_fname), 'png')

close(xfig);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SAME CATEGORY: N = N+1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ============= TIMECOURSE: TARGET-Only
cat_x_tick      = 1:(n_trs + (dur_sync*2));
cat_x_lim       = [1 (n_trs + (dur_sync*2))];
cat_x_ticklable = horzcat(1:n_trs, (1:dur_sync) + 21, (1:dur_sync) + 21);

%*************** plot selected conditions
for xsel = 1:length(sel_conds)
    
    xsel_conds = sel_conds{xsel};
    
    for xlevel = 1%:2
        
        if xlevel==1
            it_sames = 1:2; it_subplot = 1:2;
        else
            it_sames = 1:3; it_subplot = [1 3 2];
        end
        
        fig_rect   = [0 0 rect_w*4 rect_h*(length(it_sames))];
        xfig       = figure;
        set(xfig, 'Position', fig_rect)
        
        for xsame = it_sames
            
            subplot(length(it_sames), 1, it_subplot(xsame))
            
            xrsa = grp_rsa{xph}.timecourse.sameseq{xlevel}{xsame};
            
            %*************** time window
            for xwin = 1:3 %1_tr, 2_sync, 3_next item
                if xwin==1
                    it_trs       = n_trs;
                    x_tick       = 1:it_trs;
                    fitx         = linspace(1, n_trs, n_trs*10);
                elseif xwin==2
                    it_trs       = dur_sync;
                    x_tick       = (1:dur_sync) + n_trs;
                    fitx         = linspace(n_trs + 1, n_trs + dur_sync, dur_sync*10);
                elseif xwin==3
                    it_trs       = dur_sync;
                    x_tick       = (1:dur_sync) + (n_trs + dur_sync);
                    fitx         = linspace(n_trs + dur_sync + 1, n_trs + (dur_sync*2), dur_sync*10);
                end
                
                %*************** mean line plots
                clear fity xmean xstd
                
                for xcond = xsel_conds
                    if xwin == 1
                        xmean{xcond} = xrsa.cond{xcond}.mean_w(1, 1:it_trs);
                        xse{xcond}   = xrsa.cond{xcond}.se_w(1, 1:it_trs);
                    elseif xwin == 2
                        xmean{xcond} = xrsa.sync.cond{xcond}.mean_w(1, 1:it_trs);
                        xse{xcond}   = xrsa.sync.cond{xcond}.se_w(1, 1:it_trs);
                    elseif xwin == 3
                        xmean{xcond} = xrsa.sync.next.cond{xcond}.mean_w(1, 1:it_trs);
                        xse{xcond}   = xrsa.sync.next.cond{xcond}.se_w(1, 1:it_trs);
                    end
                    
                    fity{xcond}  = interp1(x_tick, xmean{xcond}, fitx,'spline');
                    
                    plot(fitx, fity{xcond}, '-','Color', xcond_color{xcond}, 'LineWidth', 2); hold on;
                    
                end
                
                %*************** legend
                clear xlegend
                for xcond = 1:length(xsel_conds)
                    it_cond = xsel_conds(xcond);
                    xlegend{xcond} = conds_names{it_cond};
                end
                lg             = legend(xlegend);
                lg.Location    = 'southeast';%'best';%%'bestoutside';
                lg.Position    = [0 0 5 0.1];% + [n_trs + 2 y_lim(1) + 0.1 n_trs+2 y_lim(1) + 0.1];
                lg.FontSize    = f_size/1.5;
                legend(xlegend,'AutoUpdate','off')
                grid on
                
                if xfill
                    %*************** std error-bar filling
                    for xcond = xsel_conds %#ok<*UNRCH>
                        clear xerr fit_err
                        xerr(1,:) = xmean{xcond} - xse{xcond};
                        xerr(2,:) = xmean{xcond} + xse{xcond};
                        
                        for i = 1:2, fit_err(:,i) = interp1(x_tick, xerr(i,:), fitx,'spline'); end
                        
                        in_between = [fit_err(:,1)', fliplr(fit_err(:,2)')];
                        h = fill([fitx fliplr(fitx)], in_between, xcond_color{xcond});
                        set(h,'facealpha', .4)
                        
                        plot(fitx, fit_err(:,1), '-', 'Color', xcond_color{xcond})
                        plot(fitx, fit_err(:,2), '-', 'Color', xcond_color{xcond})
                    end
                end
                
                %*************** stats stars
                n = 0;
                
                for i = 1:(length(xsel_conds)-1)
                    xcol = xsel_conds(i);
                    for j = i+1:length(xsel_conds)
                        xrow = xsel_conds(j);
                        
                        yy_sig = (y_lim(2)-0.02) - (n * 0.02);
                        
                        if xsel==1
                            text(-2, yy_sig, sprintf('%d vs. %d', xcol, xrow), 'FontSize', f_size);
                        end
                        
                        for xblk = 1:it_trs/n_tr_blks
                            
                            if xwin == 1
                                xpvalue = xrsa.targs.blk{xblk}.pvalue_w(xrow, xcol);
                                xx_blk = xblk;
                            elseif xwin == 2
                                xpvalue = xrsa.sync.targs.blk{xblk}.pvalue_w(xrow, xcol);
                                xx_blk = xblk + (n_trs/n_tr_blks);
                            elseif xwin == 3
                                xpvalue = xrsa.sync.next.targs.blk{xblk}.pvalue_w(xrow, xcol);
                                xx_blk = xblk + (n_trs/n_tr_blks) + (dur_sync/n_tr_blks);
                            end
                            
                            xx_sig = (n_tr_blks/2) + (n_tr_blks * (xx_blk-1)) + 0.5 - 0.08;
                            
                            if xpvalue <= xalpha
                                text(xx_sig, yy_sig, '*', 'FontSize', f_size*1.5);
                            end
                            
                            xx = n_tr_blks + (n_tr_blks * (xx_blk-1)) + 0.5;
                            plot([xx xx], y_lim, '--', 'Color', 'b')
                            
                        end
                        
                        n = n + 1;
                    end
                end
            end
            
            %*************** real onset lines
            plot([dur_stim dur_stim]+1, y_lim, '-', 'Color', xonset_color)
            plot([(dur_stim + dur_manip) (dur_stim + dur_manip)]+1, y_lim, '-', 'Color', xonset_color)
            plot([(dur_stim + dur_manip + dur_fix) (dur_stim + dur_manip + dur_fix)]+1, y_lim, '--', 'Color', xonset_color, 'Linewidth', 2)
            plot([n_trs n_trs]+1, y_lim, '-', 'Color', 'k', 'Linewidth', 3)
            plot([(n_trs + dur_sync) (n_trs + dur_sync)]+1, y_lim, '-', 'Color', 'k', 'Linewidth', 3)
            
            h1 = text(1.5, y_lim(1)+0.05, '(1) stim onset', 'Color', xonset_color, 'FontSize', f_size);
            h2 = text(dur_stim + 1.5, y_lim(1)+0.05, '(2) operation onset', 'Color', xonset_color, 'FontSize', f_size);
            h3 = text(dur_stim + dur_manip + 1.5, y_lim(1)+0.05, '(3) fixation onset', 'Color', xonset_color, 'FontSize', f_size);
            h4 = text(dur_stim + dur_manip + dur_fix + 1.5, y_lim(1)+0.05, 'synch point', 'Color', 'k', 'FontSize', f_size);
            h5 = text(n_trs + 1.5, y_lim(1)+0.05, 'synch on next-trial onset', 'Color', 'k', 'FontSize', f_size);
            h6 = text((n_trs + dur_sync) + 1.5, y_lim(1)+0.05, 'next item', 'Color', 'k', 'FontSize', f_size);
            set(h1,'Rotation', 90); set(h2,'Rotation', 90); set(h3,'Rotation', 90);
            set(h4,'Rotation', 90); set(h5,'Rotation', 90); set(h6,'Rotation', 90);
            
            %*************** set gca
            set(gca,'xlim', cat_x_lim, 'ylim', y_lim);
            set(gca,'XTick', cat_x_tick,'YTick', y_tick)
            set(gca,'XTickLabel', cat_x_ticklable);
            set(gca,'FontSize', f_size/1.5)
            
            title(sprintf('RSA item targets %s: %s (N=%s, p<%1.4f)', ...
                level_names{xlevel}, same_names{xsame}, num2str(n_subjs), xalpha),...
                'FontSize',f_size*1.5,'FontWeight','bold');
            xlabel('Volume (tr)', 'FontSize', f_size);
            ylabel('Similarity (Pearson R)', 'FontSize', f_size);
        end
        
        %*************** save fig
        fig_fname = fullfile(output_dir, sprintf('plot_tc_sameseq_targs_%s_sel%s_blk%s_fill%s_%s_n%s', ...
            level_names{xlevel}, num2str(xsel), num2str(n_tr_blks),  num2str(xfill), basename, num2str(n_subjs)));
        
        %             savefig(xfig, sprintf('%s.fig', fig_fname));
        set(gcf,'PaperUnits','inches','PaperPosition', fig_rect/100)
        print('-dpng', '-r300', sprintf('%s.png',fig_fname))
        saveas(xfig, sprintf('%s.png',fig_fname), 'png')
        
        close(xfig);
        
    end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SAME vs. DIFFERENT: N = N+1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ============= TIMECOURSE: TARGET-Only
%*************** plot selected conditions
xph      = 2;
xlevel   = 1;
it_sames = 1:2;
rect_w   = 250;
rect_h   = 450;

fig_rect = [0 0 rect_w*4 rect_h*(length(it_sames))];
xfig     = figure;
set(xfig, 'Position', fig_rect)

for xcond = 1:n_condition
    
    subplot(n_condition, 1, xcond)
    
    %*************** legend
    for xsame = it_sames
        plot([0 0], [0 0], '-','Color', xcolor{xsame}, 'LineWidth', 2); hold on
    end
    
    clear xlegend
    xlegend        = same_names(it_sames);
    lg             = legend(xlegend);
    lg.Location    = 'southeast';%'best';%%'bestoutside';
    lg.FontSize    = f_size/1.5;
    legend(xlegend,'AutoUpdate','off')
    grid on
    
    for xsame = it_sames
        
        xrsa = grp_rsa{xph}.timecourse.sameseq{xlevel}{xsame};
        
        %*************** time window
        for xwin = 1:3 %1_tr, 2_sync, 3_next item
            if xwin==1
                it_trs = n_trs;
                x_tick = 1:it_trs;
                fitx   = linspace(1, n_trs, n_trs*10);
            elseif xwin==2
                it_trs = dur_sync;
                x_tick = (1:dur_sync) + n_trs;
                fitx   = linspace(n_trs + 1, n_trs + dur_sync, dur_sync*10);
            elseif xwin==3
                it_trs = dur_sync;
                x_tick = (1:dur_sync) + (n_trs + dur_sync);
                fitx   = linspace(n_trs + dur_sync + 1, n_trs + (dur_sync*2), dur_sync*10);
            end
            
            %*************** mean line plots
            clear fity xmean xstd
            
            if xwin == 1
                xmean = xrsa.cond{xcond}.mean_w(1, 1:it_trs);
                xse   = xrsa.cond{xcond}.se_w(1, 1:it_trs);
            elseif xwin == 2
                xmean = xrsa.sync.cond{xcond}.mean_w(1, 1:it_trs);
                xse   = xrsa.sync.cond{xcond}.se_w(1, 1:it_trs);
            elseif xwin == 3
                xmean = xrsa.sync.next.cond{xcond}.mean_w(1, 1:it_trs);
                xse   = xrsa.sync.next.cond{xcond}.se_w(1, 1:it_trs);
            end
            
            fity = interp1(x_tick, xmean, fitx,'spline');
            plot(fitx, fity, '-','Color', xcolor{xsame}, 'LineWidth', 2); hold on;
            
            if xfill
                %*************** std error-bar filling
                clear xerr fit_err
                xerr(1,:) = xmean - xse;
                xerr(2,:) = xmean + xse;
                
                for i = 1:2, fit_err(:,i) = interp1(x_tick, xerr(i,:), fitx,'spline'); end
                
                in_between = [fit_err(:,1)', fliplr(fit_err(:,2)')];
                h = fill([fitx fliplr(fitx)], in_between, xcolor{xsame});
                set(h,'facealpha', .4)
                
                plot(fitx, fit_err(:,1), '-', 'Color', xcolor{xsame})
                plot(fitx, fit_err(:,2), '-', 'Color', xcolor{xsame})
            end
        end
    end
    
    %*************** real onset lines
    plot([dur_stim dur_stim]+1, y_lim, '-', 'Color', xonset_color)
    plot([(dur_stim + dur_manip) (dur_stim + dur_manip)]+1, y_lim, '-', 'Color', xonset_color)
    plot([(dur_stim + dur_manip + dur_fix) (dur_stim + dur_manip + dur_fix)]+1, y_lim, '--', 'Color', xonset_color, 'Linewidth', 2)
    plot([n_trs n_trs]+1, y_lim, '-', 'Color', 'k', 'Linewidth', 3)
    plot([(n_trs + dur_sync) (n_trs + dur_sync)]+1, y_lim, '-', 'Color', 'k', 'Linewidth', 3)
    
    h1 = text(1.5, y_lim(1)+0.05, '(1) stim onset', 'Color', xonset_color, 'FontSize', f_size);
    h2 = text(dur_stim + 1.5, y_lim(1)+0.05, '(2) operation onset', 'Color', xonset_color, 'FontSize', f_size);
    h3 = text(dur_stim + dur_manip + 1.5, y_lim(1)+0.05, '(3) fixation onset', 'Color', xonset_color, 'FontSize', f_size);
    h4 = text(dur_stim + dur_manip + dur_fix + 1.5, y_lim(1)+0.05, 'synch point', 'Color', 'k', 'FontSize', f_size);
    h5 = text(n_trs + 1.5, y_lim(1)+0.05, 'synch on next-trial onset', 'Color', 'k', 'FontSize', f_size);
    h6 = text((n_trs + dur_sync) + 1.5, y_lim(1)+0.05, 'next item', 'Color', 'k', 'FontSize', f_size);
    set(h1,'Rotation', 90); set(h2,'Rotation', 90); set(h3,'Rotation', 90);
    set(h4,'Rotation', 90); set(h5,'Rotation', 90); set(h6,'Rotation', 90);
    
    %*************** set gca
    set(gca,'xlim', cat_x_lim, 'ylim', y_lim);
    set(gca,'XTick', cat_x_tick,'YTick', y_tick)
    set(gca,'XTickLabel', cat_x_ticklable);
    set(gca,'FontSize', f_size/1.5)
    
    title(sprintf('RSA item targets: %s (N=%s, p<%1.4f)', ...
        conds_names{xcond}, num2str(n_subjs), xalpha),...
        'FontSize',f_size*1.5,'FontWeight','bold');
    xlabel('Volume (tr)', 'FontSize', f_size);
    ylabel('Similarity (Pearson R)', 'FontSize', f_size);
end

%*************** save fig
fig_fname = fullfile(output_dir, sprintf('plot_tc_samediff_targs_%s_sel%s_blk%s_fill%s_%s_n%s', ...
    level_names{xlevel}, num2str(xsel), num2str(n_tr_blks),  num2str(xfill), basename, num2str(n_subjs)));

%     savefig(xfig, sprintf('%s.fig', fig_fname));
set(gcf,'PaperUnits','inches','PaperPosition', fig_rect/100)
print('-dpng', '-r300', sprintf('%s.png',fig_fname))
saveas(xfig, sprintf('%s.png',fig_fname), 'png')

% close(xfig);

%% ============= N trials
clear n_total
for xlevel = 1:2
    if xlevel==1, it_sames = 1:2; else, it_sames = 1:3; end
    
    for xsame = it_sames
        n_total{xlevel}{xsame} = [];
    end
end

for xcond = 1:n_condition
    for xlevel = 1:2
        if xlevel==1, it_sames = 1:2; else, it_sames = 1:3; end
        
        for xsame = it_sames
            mean_trials{xlevel}(xcond, xsame) = mean(grp_rsa{xph}.timecourse.sameseq{xlevel}{xsame}.cond{xcond}.n_trial);
            std_trials{xlevel}(xcond, xsame)  = std(grp_rsa{xph}.timecourse.sameseq{xlevel}{xsame}.cond{xcond}.n_trial);
            
            n_total{xlevel}{xsame} = horzcat(n_total{xlevel}{xsame},...
                grp_rsa{xph}.timecourse.sameseq{xlevel}{xsame}.cond{xcond}.n_trial);
        end
    end
end

for xlevel = 1:2
    if xlevel==1, it_sames = 1:2; else, it_sames = 1:3; end
    
    for xsame = it_sames
        mean_trials{xlevel}(n_condition+1, xsame) = mean(n_total{xlevel}{xsame});
        std_trials{xlevel}(n_condition+1, xsame)  = std(n_total{xlevel}{xsame});
    end
end

%% *************** figure
clear h
y_lim     = [0 70];
fig_rect  = [0 0 1800 600];
xfig      = figure;

set(xfig, 'Position', fig_rect)

for xlevel = 1:2%1_category, 2_subcategory
    
    if xlevel==1
        n_sames = 2;
        xlabels = {'same','different'};
    else
        n_sames = 3;
        xlabels = {'same','different','related'};
    end
    
    for xcond = 1:(n_condition + 1)
        clear xmean xstd
        
        subplot(2, n_condition+1, xcond + ((n_condition+1) * (xlevel-1)))
        
        xmean = mean_trials{xlevel}(xcond, :);
        xstd  = std_trials{xlevel}(xcond, :);
        
        for i = 1:n_sames
            h{i} = bar(i, xmean(i)); hold on
            set(h{i}, 'facecolor', xcolor{i});
            
            text(i-0.2, 10, sprintf('m=%4.2f', xmean(i)));
            text(i-0.2, 5, sprintf('sd=%4.2f', xstd(i)));
        end
        
        for i = 1:n_sames
            errorbar(i, xmean(i), xstd(i),'k.')
        end
        legend(xlabels)
        
        set(gca,'xlim',[1 n_sames]+[-0.5 0.5],'ylim', y_lim);
        set(gca,'XTick', 1:n_sames,'YTick', y_lim(1):10:y_lim(end));
        set(gca,'XTickLabel', xlabels)
        
        if xcond <= n_condition
            title(sprintf('%s, %s', conds_names{xcond}, level_names{xlevel}));
        else
            title(sprintf('concatenated, %s', level_names{xlevel}));
        end
        ylabel('n trials');
    end
end

fig_fname = fullfile(output_dir, sprintf('plot_sameseq_sync_n_trials_%s_n%s', basename, num2str(n_subjs)));

%     savefig(xfig, sprintf('%s.fig', fig_fname));
set(gcf,'PaperUnits','inches','PaperPosition', fig_rect/100)
print('-dpng', '-r300', sprintf('%s.png',fig_fname))
saveas(xfig, sprintf('%s.png',fig_fname), 'png')

% close(xfig);

%% *************** figure
clear h xlabels
y_lim     = [0 80];
fig_rect  = [0 0 1000 300];
xfig      = figure;

set(xfig, 'Position', fig_rect)

it_conds = [1 2 4 5];
xlevel   = 1;
n_sames  = 2;
xtitle   = {'same','different'};

for xsame = 1:n_sames
    
    subplot(1, n_sames, xsame)
    
    clear xmean xstd xpvalue xstats
    xmean = mean_trials{xlevel}(:, xsame);
    xstd  = std_trials{xlevel}(:, xsame);
    
    for it = 1:length(it_conds)
        xcond = it_conds(it);
        h{it} = bar(it, xmean(xcond)); hold on
        set(h{it}, 'facecolor', xcond_color{xcond});
        
        text(it-0.2, 5, sprintf('n=%4.2f', xmean(xcond)));
        xlabels{it} = conds_names{xcond};
    end
    
    for it = 1:length(it_conds)
        xcond = it_conds(it);
        errorbar(it, xmean(xcond), xstd(xcond),'k.')
    end
    
    lg             = legend(xlabels);
    lg.Location    = 'bestoutside';%'best';%%'bestoutside';
    lg.FontSize    = f_size/1.5;
    legend(xlabels,'AutoUpdate','off')
    
    %*************** ttest
    n = 1;
    for xcol = 1:(length(it_conds)-1)
        xcol_it = it_conds(xcol);
        for xrow = xcol+1:length(it_conds)
            xrow_it = it_conds(xrow);
            xcol_n = grp_rsa{xph}.timecourse.sameseq{xlevel}{xsame}.cond{xcol_it}.n_trial;
            xrow_n = grp_rsa{xph}.timecourse.sameseq{xlevel}{xsame}.cond{xrow_it}.n_trial;
            
            [~, xpvalue, ~, xstats] = ttest(xcol_n, xrow_n, 'Alpha', 0.05);
            
            if xpvalue <= .05
                yy_sig = (y_lim(2)-5) - (n * 5);
                xx_sig = xcol + ((xrow - xcol)/2);
                n = n + 1;
                
                plot([xcol xrow], [yy_sig yy_sig], '-k')
                text(xx_sig, yy_sig, '*', 'FontSize', f_size*1.5);
                text(xx_sig, yy_sig-1, sprintf('T(%s)=%4.2f, p=%4.3f', ...
                    num2str(xstats.df), xstats.tstat, xpvalue), 'FontSize', f_size);
            end
        end
    end
    
    set(gca,'xlim',[1 length(it_conds)]+[-0.5 0.5],'ylim', y_lim);
    set(gca,'XTick', 1:length(it_conds),'YTick', y_lim(1):10:y_lim(end));
    
    title(sprintf('%s, %s', xtitle{xsame}, level_names{xlevel}));
    ylabel('n trials');
    
end

fig_fname = fullfile(output_dir, sprintf('plot_sameseq_sync_n_trials_ttest_%s_n%s', basename, num2str(n_subjs)));

%     savefig(xfig, sprintf('%s.fig', fig_fname));
set(gcf,'PaperUnits','inches','PaperPosition', fig_rect/100)
print('-dpng', '-r300', sprintf('%s.png',fig_fname))
saveas(xfig, sprintf('%s.png',fig_fname), 'png')

% close(xfig);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% REACTIVATION OF N ITEM IN N+1 TRIAL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear it_targs
% it_targs{xcond}{xsame}{xnewsame}
for xcond = [1 4 5]
    it_targs{xcond}{1} = [1 3 4 6];
    it_targs{xcond}{2} = [1 3 4 5 6];
end

it_targs{2}{1}    = [1 2 3 4 6];
it_targs{2}{2}{1} = [1 2 3 4 6];
it_targs{2}{2}{2} = [1 2 3 4 5];
it_targs{3}{1}    = [1 2 3 6];
it_targs{3}{2}    = [1 2 3 5 6];

for xcond = 1:n_condition
    for xsame = 1:2
        
        for it = 1:length(xsubj_grp)
            xsub = xsubj_grp(it);
            
            if (xcond==2) && (xsame == 2)
                for xnewsame = 1:2
                    for xtarg = it_targs{xcond}{xsame}{xnewsame}
                        t_corr = [];
                        for xtr = 1:n_trs
                            t_corr = horzcat(t_corr, ...
                                mean(grp_rsa{xph}.subj{xsub}.react.cond{xcond}.same{xsame}.newsame{xnewsame}.targ{xtarg}.tr{xtr}.corr_w));
                        end
                        
                        grp_rsa{xph}.react.cond{xcond}.same{xsame}.newsame{xnewsame}.targ{xtarg}.corr_w(it, 1:n_trs) = t_corr;
                        grp_rsa{xph}.react.cond{xcond}.same{xsame}.newsame{xnewsame}.targ{xtarg}.n_trials(it) = ...
                            size(grp_rsa{xph}.subj{xsub}.react.cond{xcond}.same{xsame}.newsame{xnewsame}.targ{xtarg}.tr{1}.corr_w, 2);
                    end
                end
                
            else
                for xtarg = it_targs{xcond}{xsame}
                    t_corr = [];
                    for xtr = 1:n_trs
                        t_corr = horzcat(t_corr, ...
                            mean(grp_rsa{xph}.subj{xsub}.react.cond{xcond}.same{xsame}.targ{xtarg}.tr{xtr}.corr_w));
                    end
                    
                    grp_rsa{xph}.react.cond{xcond}.same{xsame}.targ{xtarg}.corr_w(it, 1:n_trs) = t_corr;
                    grp_rsa{xph}.react.cond{xcond}.same{xsame}.targ{xtarg}.n_trials(it) = ...
                        size(grp_rsa{xph}.subj{xsub}.react.cond{xcond}.same{xsame}.targ{xtarg}.tr{1}.corr_w, 2);
                end
            end
        end
    end
end

%% *************** figure: timecourse
% xcond_color  = args.cond_color;

clear xtarg_color it_conds
% 1_N-target, 2_N-newtarg, 3_N+1-target,
% 4_N-relative, 5_N+1-relative, 6_others
targ_names = {'N-target','N-newtarg','N+1-target','N-relative','N+1-relative','others'};
same_names = {'same','differ'};

xtarg_color{1} = [238, 20, 91]/255;% 1_item;
xtarg_color{2} = [0, 188, 182]/255;
xtarg_color{3} = [137 62 194]/255;
xtarg_color{4} = [127, 37, 66]/255;% 3_related_item
xtarg_color{5} = [84 51 105]/255;
xtarg_color{6} = [144, 144, 144]/255;% baseline

x_tick         = 1:n_trs;
x_ticklable    = x_tick;
x_lim          = [x_tick(1) x_tick(end)];
y_lim          = [-0.3 0.3];
y_tick         = y_lim(1):0.05:y_lim(end);

it_conds{1} = 1:5;
it_conds{2} = [1 2 4 5];

for xf = 1:2
    
    n_condition = length(it_conds{xf});
    fig_rect    = [0 0 500*3 300*n_condition];
    xfig        = figure;
    set(xfig, 'Position', fig_rect)
    
    for jj = 1:length(it_conds{xf})
        xcond = it_conds{xf}(jj);
        
        for xsame = 1:2
            
            if (xcond==2) && (xsame == 2)
                it_newsame = 1:2;
            else
                it_newsame = 1;
            end
            
            for xnewsame = it_newsame
                
                if (xcond==2) && (xsame == 2)
                    xsplt   = 4 + xnewsame;
                    it_targ = it_targs{xcond}{xsame}{xnewsame};
                else
                    xsplt   = xsame + 3 * (jj-1);
                    it_targ = it_targs{xcond}{xsame};
                end
                
                subplot(n_condition, 3, xsplt)
                
                clear xlegend xrsa_corr
                %*************** legend
                for it = 1:length(it_targ)
                    xtarg = it_targ(it);
                    xlegend{it} = targ_names{xtarg};
                    plot([0 0], [0 0], '-','Color', xtarg_color{xtarg}, 'LineWidth', 2); hold on
                end
                
                lg             = legend(xlegend);
                lg.Location    = 'NorthEast';%'best';%%'bestoutside';
                lg.FontSize    = f_size/1.5;
                legend(xlegend,'AutoUpdate','off')
                grid on
                
                %*************** timecourse plot
                for xtarg = it_targ
                    if (xcond==2) && (xsame == 2)
                        xrsa_corr = grp_rsa{xph}.react.cond{xcond}.same{xsame}.newsame{xnewsame}.targ{xtarg}.corr_w;
                    else
                        xrsa_corr = grp_rsa{xph}.react.cond{xcond}.same{xsame}.targ{xtarg}.corr_w;
                    end
                    
                    %*************** mean line plots
                    clear fity xmean xse
                    xmean = mean(xrsa_corr);
                    xse   = std(xrsa_corr)/sqrt(n_subjs);
                    
                    fitx  = linspace(1, n_trs, n_trs*10);
                    fity  = interp1(x_tick, xmean, fitx,'spline');
                    plot(fitx, fity, '-','Color', xtarg_color{xtarg}, 'LineWidth', 2); hold on;
                    
                    if xfill
                        %*************** std error-bar filling
                        clear xerr fit_err
                        xerr(1,:) = xmean - xse;
                        xerr(2,:) = xmean + xse;
                        
                        for i = 1:2, fit_err(:,i) = interp1(x_tick, xerr(i,:), fitx,'spline'); end
                        
                        in_between = [fit_err(:,1)', fliplr(fit_err(:,2)')];
                        h = fill([fitx fliplr(fitx)], in_between, xtarg_color{xtarg});
                        set(h,'facealpha', .4)
                        
                        plot(fitx, fit_err(:,1), '-', 'Color', xtarg_color{xtarg})
                        plot(fitx, fit_err(:,2), '-', 'Color', xtarg_color{xtarg})
                    end
                end
                
                %*************** real onset lines
                plot([dur_stim dur_stim]+1, y_lim, '-', 'Color', xonset_color)
                plot([(dur_stim + dur_manip) (dur_stim + dur_manip)]+1, y_lim, '-', 'Color', xonset_color)
                
                h1 = text(dur_stim + 1.5, y_lim(1)+0.05, '(2) operation', 'Color', xonset_color, 'FontSize', f_size);
                h2 = text(dur_stim + dur_manip + 1.5, y_lim(1)+0.05, '(3) fixation', 'Color', xonset_color, 'FontSize', f_size);
                set(h1,'Rotation', 90); set(h2,'Rotation', 90);
                
                %*************** set gca
%                 set(gca,'xlim', x_lim, 'ylim', y_lim);
                set(gca,'XTick', x_tick,'YTick', y_tick)
                set(gca,'XTickLabel', x_ticklable);
                set(gca,'FontSize', f_size/1.5)
                
                xlabel('Volume (tr)', 'FontSize', f_size);
                ylabel('Similarity (z-R)', 'FontSize', f_size);
                
                %*************** title
                if (xcond==2) && (xsame == 2)
                    title(sprintf('iRSA N in N+1: %s\nN vs. N+1: %s, N-new vs. N+1: %s', ...
                        conds_names{xcond}, same_names{xsame}, same_names{xnewsame}),...
                        'FontSize',f_size,'FontWeight','bold');
                else
                    title(sprintf('iRSA N in N+1: %s\nN vs. N+1: %s', ...
                        conds_names{xcond}, same_names{xsame}),...
                        'FontSize',f_size,'FontWeight','bold');
                end
                
            end
        end
    end
    
    %*************** save fig
    fig_fname = fullfile(output_dir, sprintf('plot_tc_react_n+1_%s_n%s_%dconds', basename, num2str(n_subjs), n_condition));
    
    savefig(xfig, sprintf('%s.fig', fig_fname));
    set(gcf,'PaperUnits','inches','PaperPosition', fig_rect/100)
    print('-dpng', '-r300', sprintf('%s.png',fig_fname))
    saveas(xfig, sprintf('%s.png',fig_fname), 'png')
    
%     close(xfig);
end

%% *************** figure: timecourse
% xcond_color  = args.cond_color;
% 1_N-target, 2_N-newtarg, 3_N+1-target,
% 4_N-relative, 5_N+1-relative, 6_others

for xf = 2
    
    n_condition = length(it_conds{xf});
    fig_rect    = [0 0 500*3 300*n_condition];
    xfig        = figure;
    set(xfig, 'Position', fig_rect)
    
    for jj = 1:length(it_conds{xf})
        xcond = it_conds{xf}(jj);
        
        for xsame = 1:2
            
            if (xcond==2) && (xsame == 2)
                it_newsame = 1:2;
            else
                it_newsame = 1;
            end
            
            for xnewsame = it_newsame
                if (xcond==2) && (xsame == 2) 
                    xsplt   = 4 + xnewsame;
                    if (xnewsame==2)
                        it_targ = [1 3];
                    else
                        it_targ = [1 3 6];
                    end
                else
                    xsplt   = xsame + 3 * (jj-1);
                    it_targ = [1 3 6];
                end
                
                subplot(n_condition, 3, xsplt)
                
                clear xlegend xrsa_corr
                %*************** legend
                for it = 1:length(it_targ)
                    xtarg = it_targ(it);
                    xlegend{it} = targ_names{xtarg};
                    plot([0 0], [0 0], '-','Color', xtarg_color{xtarg}, 'LineWidth', 2); hold on
                end
                
                lg             = legend(xlegend);
                lg.Location    = 'NorthEast';%'best';%%'bestoutside';
                lg.FontSize    = f_size/1.5;
                legend(xlegend,'AutoUpdate','off')
                grid on
                
                %*************** timecourse plot
                for xtarg = it_targ
                    if (xcond==2) && (xsame == 2)
                        xrsa_corr = grp_rsa{xph}.react.cond{xcond}.same{xsame}.newsame{xnewsame}.targ{xtarg}.corr_w;
                    else
                        xrsa_corr = grp_rsa{xph}.react.cond{xcond}.same{xsame}.targ{xtarg}.corr_w;
                    end
                    
                    %*************** mean line plots
                    clear fity xmean xse
                    xmean = mean(xrsa_corr);
                    xse   = std(xrsa_corr)/sqrt(n_subjs);
                    
                    fitx  = linspace(1, n_trs, n_trs*10);
                    fity  = interp1(x_tick, xmean, fitx,'spline');
                    plot(fitx, fity, '-','Color', xtarg_color{xtarg}, 'LineWidth', 2); hold on;
                    
                    if xfill
                        %*************** std error-bar filling
                        clear xerr fit_err
                        xerr(1,:) = xmean - xse;
                        xerr(2,:) = xmean + xse;
                        
                        for i = 1:2, fit_err(:,i) = interp1(x_tick, xerr(i,:), fitx,'spline'); end
                        
                        in_between = [fit_err(:,1)', fliplr(fit_err(:,2)')];
                        h = fill([fitx fliplr(fitx)], in_between, xtarg_color{xtarg});
                        set(h,'facealpha', .4)
                        
                        plot(fitx, fit_err(:,1), '-', 'Color', xtarg_color{xtarg})
                        plot(fitx, fit_err(:,2), '-', 'Color', xtarg_color{xtarg})
                    end
                end
                
                %*************** real onset lines
                plot([dur_stim dur_stim]+1, y_lim, '-', 'Color', xonset_color)
                plot([(dur_stim + dur_manip) (dur_stim + dur_manip)]+1, y_lim, '-', 'Color', xonset_color)
                
                h1 = text(dur_stim + 1.5, y_lim(1)+0.05, '(2) operation', 'Color', xonset_color, 'FontSize', f_size);
                h2 = text(dur_stim + dur_manip + 1.5, y_lim(1)+0.05, '(3) fixation', 'Color', xonset_color, 'FontSize', f_size);
                set(h1,'Rotation', 90); set(h2,'Rotation', 90);
                
                %*************** set gca
%                 set(gca,'xlim', x_lim, 'ylim', y_lim);
                set(gca,'XTick', x_tick,'YTick', y_tick)
                set(gca,'XTickLabel', x_ticklable);
                set(gca,'FontSize', f_size/1.5)
                
                xlabel('Volume (tr)', 'FontSize', f_size);
                ylabel('Similarity (z-R)', 'FontSize', f_size);
                
                %*************** title
                if (xcond==2) && (xsame == 2)
                    title(sprintf('iRSA N in N+1: %s\nN vs. N+1: %s, N-new vs. N+1: %s', ...
                        conds_names{xcond}, same_names{xsame}, same_names{xnewsame}),...
                        'FontSize',f_size,'FontWeight','bold');
                else
                    title(sprintf('iRSA N in N+1: %s\nN vs. N+1: %s', ...
                        conds_names{xcond}, same_names{xsame}),...
                        'FontSize',f_size,'FontWeight','bold');
                end
                
            end
        end
    end
    
    %*************** save fig
    fig_fname = fullfile(output_dir, sprintf('plot_tc_react_n+1_targs_%s_n%s_%dconds', basename, num2str(n_subjs), n_condition));
    
    savefig(xfig, sprintf('%s.fig', fig_fname));
    set(gcf,'PaperUnits','inches','PaperPosition', fig_rect/100)
    print('-dpng', '-r300', sprintf('%s.png',fig_fname))
    saveas(xfig, sprintf('%s.png',fig_fname), 'png')
    
%     close(xfig);
end

%% *************** figure: bar
xwindow        = args.bin_trs{2};%peack of targ (for both N & N+1)

y_lim          = [-0.05 0.6];
y_tick         = y_lim(1):0.05:y_lim(end);

for xf = 1:2
    
    n_condition = length(it_conds{xf});
    fig_rect    = [0 0 500*3 500*n_condition];
    xfig        = figure;
    set(xfig, 'Position', fig_rect)
    
    for jj = 1:length(it_conds{xf})
        xcond = it_conds{xf}(jj);
        for xsame = 1:2
            
            if (xcond==2) && (xsame == 2)
                it_newsame = 1:2;
            else
                it_newsame = 1;
            end
            
            for xnewsame = it_newsame
                
                if (xcond==2) && (xsame == 2)
                    xsplt   = 4 + xnewsame;
                    it_targ = it_targs{xcond}{xsame}{xnewsame};
                else
                    xsplt   = xsame + 3 * (jj-1);
                    it_targ = it_targs{xcond}{xsame};
                end
                
                x_tick = 1:length(it_targ);
                x_lim  = [x_tick(1) x_tick(end)] + [-1 1];
                
                subplot(n_condition, 3, xsplt)
                
                clear xlegend xrsa_corr xmean xse xbin_alpha
                
                for it = 1:length(it_targ)
                    xtarg = it_targ(it);
                    
                    if (xcond==2) && (xsame == 2)
                        xrsa_corr = mean(grp_rsa{xph}.react.cond{xcond}.same{xsame}.newsame{xnewsame}.targ{xtarg}.corr_w(:,xwindow),2);
                        n_trial   = mean(grp_rsa{xph}.react.cond{xcond}.same{xsame}.newsame{xnewsame}.targ{xtarg}.n_trials);
                    else
                        xrsa_corr = mean(grp_rsa{xph}.react.cond{xcond}.same{xsame}.targ{xtarg}.corr_w(:,xwindow),2);
                        n_trial   = mean(grp_rsa{xph}.react.cond{xcond}.same{xsame}.targ{xtarg}.n_trials);
                    end
                    
                    %*************** bar plot
                    xmean(it) = mean(xrsa_corr);
                    xse(it)   = std(xrsa_corr)/sqrt(n_subjs);
                    
                    b = bar(it, xmean(it)); hold on
                    set(b, 'facecolor', xtarg_color{xtarg});
                    
                    %*************** legend
                    xlegend{it} = targ_names{xtarg};
                end
                
                %*************** alpha level
                n_pw       = size(combntns(it_targ, 2),1);
                xbin_alpha = 0.05/n_pw;
                text (0.1, y_lim(2)-0.05, ['\alpha', sprintf('(%s)=%1.4f', num2str(n_pw), xbin_alpha)]);
                
                %*************** n trials
                text(0.1, y_lim(2)-0.1, sprintf('n trials=%1.2f', n_trial));
                
                n = 1;
                for xcol = 1:(length(it_targ)-1)
                    for xrow = (xcol+1):length(it_targ)
                        clear xp
                        xcol_targ = it_targ(xcol);
                        xrow_targ = it_targ(xrow);
                        
                        if (xcond==2) && (xsame == 2)
                            xcol_corr = mean(grp_rsa{xph}.react.cond{xcond}.same{xsame}.newsame{xnewsame}.targ{xcol_targ}.corr_w(:,xwindow),2);
                            xrow_corr = mean(grp_rsa{xph}.react.cond{xcond}.same{xsame}.newsame{xnewsame}.targ{xrow_targ}.corr_w(:,xwindow),2);
                        else
                            xcol_corr = mean(grp_rsa{xph}.react.cond{xcond}.same{xsame}.targ{xcol_targ}.corr_w(:,xwindow),2);
                            xrow_corr = mean(grp_rsa{xph}.react.cond{xcond}.same{xsame}.targ{xrow_targ}.corr_w(:,xwindow),2);
                        end
                        
                        [~, xp] = ttest(xcol_corr, xrow_corr, 'Alpha', xbin_alpha);
                        
                        %*************** ttest
                        if xp <= xbin_alpha%xbin_alpha
                            n = n+1;
                            yy_sig = y_lim(2)-0.1 - (0.03 * n);
                            
                            if (xp <= xbin_alpha) && (xp > (xbin_alpha/5))
                                xsig = '*';
                            elseif (xp <= (xbin_alpha/5)) && (xp > (xbin_alpha/50))
                                xsig = '**';
                            elseif xp <= (xbin_alpha/50)%0.001
                                xsig = '***';
                            end
                            
                            plot([xcol xrow],[yy_sig yy_sig]-0.005,'k-');
                            text(xcol+(xrow-xcol)/2 - (length(xsig)*0.05), yy_sig, xsig, 'FontSize', f_size+10);
%                             text(xcol+(xrow-xcol)/2 + 0.3, yy_sig-0.015, sprintf('p=%4.4f', xp), 'FontSize', 10);
                        end
                    end
                end
                
                %*************** legend
                lg          = legend(xlegend);
                lg.Location = 'NorthEast';%'SouthEast';
                lg.FontSize = 8;
                legend(xlegend,'AutoUpdate','off')
                
                errorbar(1:length(it_targ), xmean, xse,'k.')
                
                %*************** set gca
                set(gca,'xlim', x_lim, 'ylim', y_lim);
                set(gca,'XTick', x_tick,'YTick', y_tick)
                set(gca,'FontSize', f_size/1.5)
                
                xlabel('targets', 'FontSize', f_size);
                ylabel('Similarity (z-R)', 'FontSize', f_size);
                
                %*************** title
                if (xcond==2) && (xsame == 2)
                    title(sprintf('iRSA N in N+1: %s: %s-%sTR\nN vs. N+1: %s, N-new vs. N+1: %s', ...
                        conds_names{xcond}, num2str(xwindow(1)), num2str(xwindow(end)), ...
                        same_names{xsame}, same_names{xnewsame}),...
                        'FontSize',f_size,'FontWeight','bold');
                else
                    title(sprintf('iRSA N in N+1: %s: %s-%sTR\nN vs. N+1: %s', ...
                        conds_names{xcond}, num2str(xwindow(1)), num2str(xwindow(end)), ...
                        same_names{xsame}),...
                        'FontSize',f_size,'FontWeight','bold');
                end
                
            end
        end
    end
    
    %*************** save fig
    fig_fname = fullfile(output_dir, sprintf('plot_bar_react_n+1_%s_n%s_%dconds', basename, num2str(n_subjs), n_condition));
    
    savefig(xfig, sprintf('%s.fig', fig_fname));
    set(gcf,'PaperUnits','inches','PaperPosition', fig_rect/100)
    print('-dpng', '-r300', sprintf('%s.png',fig_fname))
    saveas(xfig, sprintf('%s.png',fig_fname), 'png')
    
%     close(xfig);
end

%% *************** figure: bar: N-targ vs. N+1-targ
% 1_N-target, 2_N-newtarg, 3_N+1-target,
% 4_N-relative, 5_N+1-relative, 6_others
clear it_targ
% it_targs{xcond}{xsame}{xnewsame}
xwindow = args.bin_trs{2};%peack of targ (for both N & N+1)

if strcmp(args.mask_name, 'hippocampus')
    y_lim   = [-0.05 0.2];
else
    y_lim   = [-0.05 0.3];
end

y_tick  = y_lim(1):0.05:y_lim(end);

for xf =  2
    
    n_condition = length(it_conds{xf});
    fig_rect    = [0 0 500*3 500*n_condition];
    xfig        = figure;
    set(xfig, 'Position', fig_rect)
    
    for jj = 1:length(it_conds{xf})
        xcond = it_conds{xf}(jj);
        for xsame = 1:2
            
            if (xcond==2) && (xsame == 2)
                it_newsame = 1:2;
            else
                it_newsame = 1;
            end
            
            for xnewsame = it_newsame
                
                if (xcond==2) && (xsame == 2) 
                    xsplt   = 4 + xnewsame;
                    if (xnewsame==2)
                        it_targ = [1 3];
                    else
                        it_targ = [1 3 6];
                    end
                else
                    xsplt   = xsame + 3 * (jj-1);
                    it_targ = [1 3 6];
                end
                
                x_tick = 1:length(it_targ);
                x_lim  = [x_tick(1) x_tick(end)] + [-1 1];
                
                subplot(n_condition, 3, xsplt)
                
                clear xlegend xrsa_corr xmean xse xbin_alpha
                
                for it = 1:length(it_targ)
                    xtarg = it_targ(it);
                    
                    if (xcond==2) && (xsame == 2)
                        xrsa_corr = mean(grp_rsa{xph}.react.cond{xcond}.same{xsame}.newsame{xnewsame}.targ{xtarg}.corr_w(:,xwindow),2);
                        n_trial   = mean(grp_rsa{xph}.react.cond{xcond}.same{xsame}.newsame{xnewsame}.targ{xtarg}.n_trials);
                    else
                        xrsa_corr = mean(grp_rsa{xph}.react.cond{xcond}.same{xsame}.targ{xtarg}.corr_w(:,xwindow),2);
                        n_trial   = mean(grp_rsa{xph}.react.cond{xcond}.same{xsame}.targ{xtarg}.n_trials);
                    end
                    
                    %*************** bar plot
                    xmean(it) = mean(xrsa_corr);
                    xse(it)   = std(xrsa_corr)/sqrt(n_subjs);
                    
                    b = bar(it, xmean(it)); hold on
                    set(b, 'facecolor', xtarg_color{xtarg});
                    
                    %*************** legend
                    xlegend{it} = targ_names{xtarg};
                end
                
                %*************** alpha level
                n_pw       = size(combntns(it_targ, 2),1);
                xbin_alpha = 0.05/n_pw;
%                 text (0.1, y_lim(2)-0.05, ['\alpha', sprintf('(%s)=%1.4f', num2str(n_pw), xbin_alpha)]);
                
                %*************** n trials
                text(0.1, y_lim(2)-0.1, sprintf('n trials=%1.2f', n_trial));
                
                n = 1;
                for xcol = 1:(length(it_targ)-1)
                    for xrow = (xcol+1):length(it_targ)
                        clear xp
                        xcol_targ = it_targ(xcol);
                        xrow_targ = it_targ(xrow);
                        
                        if (xcond==2) && (xsame == 2)
                            xcol_corr = mean(grp_rsa{xph}.react.cond{xcond}.same{xsame}.newsame{xnewsame}.targ{xcol_targ}.corr_w(:,xwindow),2);
                            xrow_corr = mean(grp_rsa{xph}.react.cond{xcond}.same{xsame}.newsame{xnewsame}.targ{xrow_targ}.corr_w(:,xwindow),2);
                        else
                            xcol_corr = mean(grp_rsa{xph}.react.cond{xcond}.same{xsame}.targ{xcol_targ}.corr_w(:,xwindow),2);
                            xrow_corr = mean(grp_rsa{xph}.react.cond{xcond}.same{xsame}.targ{xrow_targ}.corr_w(:,xwindow),2);
                        end
                        
                        [~, xp] = ttest(xcol_corr, xrow_corr, 'Alpha', xbin_alpha);
                        
                        %*************** ttest
                        if xp <= xbin_alpha%xbin_alpha
                            n = n+1;
                            yy_sig = y_lim(2) - (0.03 * n);
                            
                            if (xp <= xbin_alpha) && (xp > (xbin_alpha/5))
                                xsig = '*';
                            elseif (xp <= (xbin_alpha/5)) && (xp > (xbin_alpha/50))
                                xsig = '**';
                            elseif xp <= (xbin_alpha/50)%0.001
                                xsig = '***';
                            end
                            
                            plot([xcol xrow],[yy_sig yy_sig]-0.005,'k-');
                            text(xcol+(xrow-xcol)/2 - (length(xsig)*0.05), yy_sig, xsig, 'FontSize', f_size+10);
                        end
                    end
                end
                
                %*************** legend
                if (xcond==1) && (xsame==1)
                    lg          = legend(xlegend);
                    lg.Location = 'NorthEast';%'SouthEast';
                    lg.FontSize = 8;
                    legend(xlegend,'AutoUpdate','off')
                end
                
                errorbar(1:length(it_targ), xmean, xse,'k.')
                
                %*************** set gca
                set(gca,'xlim', x_lim, 'ylim', y_lim);
                set(gca,'XTick', x_tick,'YTick', y_tick)
                set(gca,'FontSize', f_size/1.5)
                
                xlabel('targets', 'FontSize', f_size);
                ylabel('Similarity (z-R)', 'FontSize', f_size);
                
                %*************** title
                if (xcond==2) && (xsame == 2)
                    title(sprintf('iRSA N in N+1: %s: %s-%sTR\nN vs. N+1: %s, N-new vs. N+1: %s', ...
                        conds_names{xcond}, num2str(xwindow(1)), num2str(xwindow(end)), ...
                        same_names{xsame}, same_names{xnewsame}),...
                        'FontSize',f_size,'FontWeight','bold');
                else
                    title(sprintf('iRSA N in N+1: %s: %s-%sTR\nN vs. N+1: %s', ...
                        conds_names{xcond}, num2str(xwindow(1)), num2str(xwindow(end)), ...
                        same_names{xsame}),...
                        'FontSize',f_size,'FontWeight','bold');
                end
                
            end
        end
    end
    
    %*************** save fig
    fig_fname = fullfile(output_dir, sprintf('plot_bar_react_n+1_targs_%s_n%s_%dconds', basename, num2str(n_subjs), n_condition));
    
    savefig(xfig, sprintf('%s.fig', fig_fname));
    set(gcf,'PaperUnits','inches','PaperPosition', fig_rect/100)
    print('-dpng', '-r300', sprintf('%s.png',fig_fname))
    saveas(xfig, sprintf('%s.png',fig_fname), 'png')
    
%     close(xfig);
end

%% ============= SAVE GROUP RSA
g_fname = fullfile(dirs.rsa.group.pattern, sprintf('group_%s_blk%d.mat', basename, n_tr_blks));

fprintf('\n....saving grp_rsa patterns\n')
save(g_fname, 'grp_rsa','-v7.3')

g_fsize = dir(g_fname);
fprintf('file size: %s GB\n', num2str(g_fsize.bytes/(10^9)));

end
