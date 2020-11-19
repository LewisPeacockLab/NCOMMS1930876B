function[] = clearmem_rsa_decode_04_bootstrap(args, dirs)
%*************** bootstrapping
%% ============= SETUP DIRECTORY
xph             = 2;
xsubj_grp       = args.filtered_subs;
n_subjs         = length(xsubj_grp);
output_dir      = sprintf('%s/bootstrap', dirs.rsa.group.parse);
if ~isdir(output_dir), mkdir(output_dir); end

%% ============= SETUP PARAMETERS
conds_names     = args.index{xph}.param.conds_names;
same_names      = {'same','different','same-diff'};
spm_p_thresh    = args.spm_p_thresh;
spm_v_extent    = args.spm_v_extent;
spm_thresh_type = args.spm_thresh_type;
rsa_mask        = args.rsa_mask;

it_conds = [1 2 4 5];
n_conds  = length(it_conds);

xtarg       = 1;% 1_target, 2_nontarget, 3_baseline
xlevel      = 1;% 1_cate, 2_subcate
xbin        = 2;% 1_n trial (args.bin_trs{1}), 2_n+1 trial (args.bin_trs{2}), 3_n+1 baseline
it_sames    = 1:3;% 1_same, 2_differ, 3_same-differ
n_sample    = args.n_sample;
n_iteration = args.n_iteration;
xalpha      = 0.05;

xcond_color = args.cond_color;
xcolor{1}   = [232 14 138]/255;%same
xcolor{2}   = [101 47 142]/255;%differ

basename = sprintf('rsa_%s_spmT_%s_thresh%s_ext%s_%s_mask', ...
    args.level, spm_thresh_type, num2str(spm_p_thresh), ...
    num2str(spm_v_extent), rsa_mask);

g_fname = fullfile(dirs.rsa.group.pattern, sprintf('group_%s.mat', basename));

%% ============= LOAD GROUP STRUCTURE
fprintf('\n....loading grp_rsa patterns\n')
xx = load(g_fname);%grp_rsa.subj{xsub}

g_fsize = dir(g_fname);
fprintf('file size: %s GB\n', num2str(g_fsize.bytes/(10^9)));

xsubj = xx.grp_rsa{xph}.subj;
clear grp_rsa

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%*************** STUDY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ============= ORIGINAL MEAN VALUE: MU
clear grp_mean_corr grp_mean_corr_rep

for xsame = it_sames
    for it = 1:n_conds
        xcond = it_conds(it);
        
        for it_sub = 1:n_subjs
            xsub = xsubj_grp(it_sub);
            clear xrsa xrep_rsa
            
            if xsame~=3
                if (xcond == 2) && (xsame == 2)
                    xrsa     = xsubj{xsub}.timecourse.sameseq{xlevel}{xsame}.sync.next.cond{xcond}.newitem{2}.targ{xtarg};
                    xrep_rsa = xsubj{xsub}.timecourse.sameseq{xlevel}{xsame}.sync.next.cond{xcond}.newitem{1}.targ{xtarg};
                else
                    xrsa     = xsubj{xsub}.timecourse.sameseq{xlevel}{xsame}.sync.next.cond{xcond}.targ{xtarg};
                end
                
                n_trials = length(xrsa.tr{1}.corr_w);
                
                for xtrial = 1:n_trials
                    xbin_corr = [];
                    for it_tr = 1:length(args.bin_trs{xbin})
                        xtr = args.bin_trs{xbin}(it_tr);
                        
                        if length(xrsa.tr{xtr}.corr_w) >= xtrial
                            xbin_corr(it_tr) = xrsa.tr{xtr}.corr_w(xtrial); %#ok<*AGROW>
                        end
                    end
                    
                    if ~isempty(xbin_corr)
                        xtrial_corr(xtrial) = mean(xbin_corr);
                    else
                        fprintf('\n ... no data: cond %d, same %d, sub %s, trial %s\n',...
                            xcond, xsame, num2str(xsub), num2str(xtrial));
                    end
                end
                
                grp_mean_corr{xsame}{it}(it_sub) = mean(xtrial_corr);
            else
                grp_mean_corr{xsame}{it}(it_sub) = ...
                    grp_mean_corr{1}{it}(it_sub) - grp_mean_corr{2}{it}(it_sub);
            end
            
            %************** XREPLACE: SAMPLING
            if (xcond == 2) && (xsame == 2)
                clear n_trials xsample
                n_trials = length(xrep_rsa.tr{1}.corr_w);
                
                for xtrial = 1:n_trials
                    for it_tr = 1:length(args.bin_trs{xbin})
                        xtr = args.bin_trs{xbin}(it_tr);
                        
                        if length(xrep_rsa.tr{xtr}.corr_w) >= xtrial
                            xbin_corr_rep(it_tr) = xrep_rsa.tr{xtr}.corr_w(xtrial);
                        end
                    end
                    
                    if ~isempty(xbin_corr)
                        xtrial_corr_rep(xtrial) = mean(xbin_corr_rep);
                    else
                        fprintf('\n ... no data: cond %d, same %d, sub %s, trial %s\n',...
                            xcond, xsame, num2str(xsub), num2str(xtrial));
                    end
                end
                
                grp_mean_corr_rep(it_sub) = mean(xtrial_corr_rep);
            end
        end
        
        xmu_rsa.mean(it, xsame) = mean(grp_mean_corr{xsame}{it});
        xmu_rsa.se(it, xsame)   = std(grp_mean_corr{xsame}{it})/sqrt(n_subjs);
        
        if (xcond == 2) && (xsame == 2)
            xmu_rsa_rep.mean = mean(grp_mean_corr_rep);
            xmu_rsa_rep.se   = std(grp_mean_corr_rep)/sqrt(n_subjs);
        end
    end
    
    %*************** ttest
    xmu_rsa.pvalue{xsame} = nan(n_conds, n_conds);
    for xcol = 1:n_conds
        for xrow = 1:n_conds
            
            col_rsa = grp_mean_corr{xsame}{xcol};
            row_rsa = grp_mean_corr{xsame}{xrow};
            
            [~, xmu_rsa.pvalue{xsame}(xcol, xrow), ~, xstats] = ...
                ttest(col_rsa, row_rsa, 'Alpha', xalpha);
            xmu_rsa.tvalue{xsame}(xcol, xrow) = xstats.tstat;
        end
    end
end

%*************** ttest replace
col_rsa = grp_mean_corr{2}{2};
row_rsa = grp_mean_corr_rep;
[~, xmu_rsa_rep.pvalue, ~, xstats] = ttest(col_rsa, row_rsa, 'Alpha', xalpha);
xmu_rsa_rep.tvalue = xstats.tstat;

%% ============= BOOTSTRAPPING
%*************** N & N+1: SAME/DIFF
%*************** Bin mean: 11:16 from N+1 onset
%*************** sample (replace) 15 trials
clear xboot_rsa xboot_rsa_rep

for xsame = it_sames
    xboot_rsa.pvalue{xsame} = cell(n_conds, n_conds);
    xboot_rsa.tvalue{xsame} = cell(n_conds, n_conds);
end

for xboot = 1:n_iteration
    
    if ~(mod(xboot-1, 20))
        fprintf('%s ', num2str(xboot));
    end
    
    clear grp_mean_corr grp_mean_corr_rep
    
    for xsame = it_sames
        for it = 1:n_conds
            xcond = it_conds(it);
            clear xbin_corr xtrial_corr xbin_corr_rep xtrial_corr_rep
            
            for it_sub = 1:n_subjs
                xsub = xsubj_grp(it_sub);
                clear xrsa xrep_rsa
                
                if xsame~=3
                    if (xcond == 2) && (xsame == 2)
                        xrsa     = xsubj{xsub}.timecourse.sameseq{xlevel}{xsame}.sync.next.cond{xcond}.newitem{2}.targ{xtarg};
                        xrep_rsa = xsubj{xsub}.timecourse.sameseq{xlevel}{xsame}.sync.next.cond{xcond}.newitem{1}.targ{xtarg};
                    else
                        xrsa     = xsubj{xsub}.timecourse.sameseq{xlevel}{xsame}.sync.next.cond{xcond}.targ{xtarg};
                    end
                    
                    %************** XRSA: SAMPLING
                    clear n_trials xsample
                    n_trials = length(xrsa.tr{1}.corr_w);
                    xsample  = sort(datasample(1:n_trials, n_sample, 'Replace', true));
                    
                    for it_trial = 1:n_sample
                        xtrial = xsample(it_trial);
                        
                        xbin_corr = [];
                        for it_tr = 1:length(args.bin_trs{xbin})
                            xtr = args.bin_trs{xbin}(it_tr);
                            
                            if length(xrsa.tr{xtr}.corr_w) >= xtrial
                                xbin_corr(it_tr) = xrsa.tr{xtr}.corr_w(xtrial); %#ok<*AGROW>
                            end
                        end
                        
                        if ~isempty(xbin_corr)
                            xtrial_corr(it_trial) = mean(xbin_corr);
                        else
                            fprintf('\n ... no data: cond %d, same %d, sub %s, trial %s\n',...
                                xcond, xsame, num2str(xsub), num2str(xtrial));
                        end
                    end
                    
                    grp_mean_corr{xsame}{it}(it_sub) = mean(xtrial_corr);
                else
                    grp_mean_corr{xsame}{it}(it_sub) = ...
                        grp_mean_corr{1}{it}(it_sub) - grp_mean_corr{2}{it}(it_sub);
                end
                
                %************** XREPLACE: SAMPLING
                if (xcond == 2) && (xsame == 2)
                    clear n_trials xsample
                    n_trials = length(xrep_rsa.tr{1}.corr_w);
                    xsample  = sort(datasample(1:n_trials, n_sample, 'Replace', true));
                    
                    for it_trial = 1:n_sample
                        xtrial = xsample(it_trial);
                        
                        for it_tr = 1:length(args.bin_trs{xbin})
                            xtr = args.bin_trs{xbin}(it_tr);
                            
                            if length(xrep_rsa.tr{xtr}.corr_w) >= xtrial
                                xbin_corr_rep(it_tr) = xrep_rsa.tr{xtr}.corr_w(xtrial);
                            end
                        end
                        
                        if ~isempty(xbin_corr)
                            xtrial_corr_rep(it_trial) = mean(xbin_corr_rep);
                        else
                            fprintf('\n ... no data: cond %d, same %d, sub %s, trial %s\n',...
                                xcond, xsame, num2str(xsub), num2str(xtrial));
                        end
                    end
                    
                    grp_mean_corr_rep(it_sub) = mean(xtrial_corr_rep);
                end
            end
            
            xboot_rsa.mean{xsame}{it}(xboot) = mean(grp_mean_corr{xsame}{it});
            
            if (xcond == 2) && (xsame == 2)
                xboot_rsa_rep.mean(xboot) = mean(grp_mean_corr_rep);
            end
        end
    end
end

fprintf('\n');

xfname = fullfile(output_dir, sprintf('bootstrap_%s_n%s.mat', basename, num2str(n_subjs)));
save(xfname, 'xboot_rsa','xboot_rsa_rep','xmu_rsa','xmu_rsa_rep','-v7.3')

g_fsize = dir(xfname);
fprintf('file size: %s GB\n', num2str(g_fsize.bytes/(10^9)));

%% ============= BOOTSTRAP RESULTS NORMALITY
% Anderson-Darling test: h = 0 indicates that it's from normal distribution
% null hypothesis: the data is from a population with a normal distribution
% h=0: adtest fails to reject the null hypothesis at the 5% significance level.
normality = {'normal','non-normal'};
fig_rect  = [0 0 900 1200];
xfig      = figure;
set(xfig, 'Position', fig_rect)

for it = 1:n_conds
    xcond = it_conds(it);
    for xsame = it_sames
        
        xsp = xsame + (4 * (it-1));
        subplot(n_conds, 4, xsp)
        
        qqplot(xboot_rsa.mean{xsame}{it});
        
        %*************** Anderson-Darling test
        [xh, xp, xadstat, ~] = adtest(xboot_rsa.mean{xsame}{it});
        
        title(sprintf('%s N & N+1 %s\n %s: p=%1.3f, adstats=%1.3f', ...
            conds_names{xcond}, same_names{xsame},...
            normality{xh+1}, xp, xadstat));
    end
end

%*************** replace new-item == next-item
subplot(n_conds, 4, 8)
qqplot(xboot_rsa_rep.mean);

%*************** Anderson-Darling test
[xh, xp, xadstat, ~] = adtest(xboot_rsa_rep.mean);

title(sprintf('%s N & N+1 %s: new == next\n %s: p=%1.3f, adstats=%1.3f', ...
    conds_names{2}, same_names{2},...
    normality{xh+1}, xp, xadstat));

%*************** save fig
fig_fname = fullfile(output_dir, sprintf('qqplot_normality_%s_n%s', ...
    basename, num2str(n_subjs)));

set(gcf,'PaperUnits','inches','PaperPosition', fig_rect/100)
print('-dpng', '-r300', sprintf('%s.png',fig_fname))
saveas(xfig, sprintf('%s.png',fig_fname), 'png')

close(xfig);

%% ============= BOOTSTRAP STATS
%*************** stats: ci
for xsame = it_sames
    for it = 1:n_conds
        xcond = it_conds(it);
        clear lower_ci upper_ci
        lower_ci = prctile(sort(xboot_rsa.mean{xsame}{it}), (xalpha * 100));
        upper_ci = prctile(sort(xboot_rsa.mean{xsame}{it}), ((1-xalpha) * 100));
        xboot_rsa.ci{xsame}{it} = [lower_ci upper_ci];
        
        if (xcond == 2) && (xsame == 2)
            clear lower_ci upper_ci
            lower_ci = prctile(sort(xboot_rsa_rep.mean), (xalpha * 100));
            upper_ci = prctile(sort(xboot_rsa_rep.mean), ((1-xalpha) * 100));    
            xboot_rsa_rep.ci = [lower_ci upper_ci];
        end
    end
    
    %*************** ttest
    for xcol = 1:(n_conds-1)
        for xrow = (xcol+1):n_conds
            
            xcol_rsa = xboot_rsa.mean{xsame}{xcol};
            xrow_rsa = xboot_rsa.mean{xsame}{xrow};
            
            if mean(xboot_rsa.mean{xsame}{xcol}) > mean(xboot_rsa.mean{xsame}{xrow})
                xboot_rsa.pw_percent{xsame}(xrow, xcol) = ...
                    1 - (numel(find(xcol_rsa > xrow_rsa))/n_iteration);
            else
                xboot_rsa.pw_percent{xsame}(xrow, xcol) = ...
                    1 - (numel(find(xcol_rsa < xrow_rsa))/n_iteration);
            end
        end
    end
end

%*************** ttest
xboot_rsa_rep.pw_percent = 1-(numel(find(xboot_rsa_rep.mean < xboot_rsa.mean{2}{2}))/n_iteration);

xfname = fullfile(output_dir, sprintf('bootstrap_%s_n%s', basename, num2str(n_subjs)));
save(xfname,'xmu_rsa','xmu_rsa_rep','-v7.3')

g_fsize = dir(xfname);
fprintf('file size: %s GB\n', num2str(g_fsize.bytes/(10^9)));

%% *************** figure
f_size   = 10;
fig_rect = [0 0 1200 500];

xfig = figure;
set(xfig, 'Position', fig_rect)

for xsame = it_sames

    x_lim = [1 n_conds] + [-0.5 0.5];
    if xsame ~= 3
        y_lim  = [0 0.4];
        y_off  = 0.05;
        y_tick = y_lim(1):y_off:y_lim(end);
    else
        y_lim  = [-0.06 0.09];
        y_off  = 0.01;
        y_tick = y_lim(1):y_off:y_lim(end);
    end
    
    subplot(1, 4, xsame)
    
    clear xmean xerr    
    for it = 1:n_conds
        xcond = it_conds(it);
        
        xmean(it)  = mean(xboot_rsa.mean{xsame}{it});
        xerr(:,it) = (xboot_rsa.ci{xsame}{it} - xmean(it))';
        
        b = bar(it, xmean(it)); hold on
        set(b, 'facecolor', xcond_color{xcond});
        xlegend{it} = conds_names{xcond};
    end
    
    %*************** legend
    if xsame == 2
        lg             = legend(xlegend);
        lg.Location    = 'SouthEast';%'SouthEast';
        legend(xlegend,'AutoUpdate','off')
    end
    
    errorbar(1:n_conds, xmean, xerr(1,:), xerr(2,:),'k.')
    
    %*************** ttest
    n = 1;
    for xcol = 1:(n_conds-1)
        for xrow = (xcol+1):n_conds
            xp = xboot_rsa.pw_percent{xsame}(xrow, xcol);
            
            xx_sig = xcol + (xrow-xcol)/2;
            yy_sig = y_lim(2) - (n * (y_off/2));
            n = n+1;
            
            plot([xcol xrow],[yy_sig yy_sig],'k-');
            if (xp <= xalpha)
                text(xx_sig, yy_sig, '*', 'FontSize', f_size);
            end
            text(xx_sig, yy_sig - (y_off/5), sprintf('p=%4.3f', xp), 'FontSize', 6);
        end
    end

    %*************** set
    set(gca,'xlim', x_lim,'ylim', y_lim);
    set(gca,'YTick', y_tick);
    
    title(sprintf('RSA: %s', same_names{xsame}));
    ylabel('RSA');
end

%*************** replace
clear xmean xerr
x_lim  = [1 2] + [-0.5 0.5];
y_lim  = [0 0.3];
y_off  = 0.05;
y_tick = y_lim(1):y_off:y_lim(end);

subplot(1, 4, 4)

xmean(1)  = mean(xboot_rsa_rep.mean);
xerr(:,1) = (xboot_rsa_rep.ci - xmean(1))';
xmean(2)  = mean(xboot_rsa.mean{2}{2});
xerr(:,2) = (xboot_rsa.ci{2}{2} - xmean(2))';

for xsame = 1:2
    b = bar(xsame, xmean(xsame)); hold on
    set(b, 'facecolor', xcolor{xsame});
end

%*************** legend
xlegend        = {'same','differ'};
lg             = legend(xlegend);
lg.Location    = 'SouthEast';%'SouthEast';
legend(xlegend,'AutoUpdate','off')

errorbar(1:2, xmean, xerr(1,:), xerr(2,:),'k.')

%*************** ttest
xp = xboot_rsa_rep.pw_percent;
yy_sig = y_lim(2) - (y_off/2);

plot([1 2],[yy_sig yy_sig],'k-');
if (xp <= xalpha)
    text(1.5, yy_sig, '*', 'FontSize', f_size);
end
text(1.5, yy_sig - (y_off/5), sprintf('p=%4.3f', xp), 'FontSize', 7);

%*************** set
set(gca,'xlim', x_lim,'ylim', y_lim);
set(gca,'YTick', y_tick);

title(sprintf('RSA: replace: new == next'));
ylabel('RSA');

%*************** save fig
fig_fname = fullfile(output_dir, sprintf('bs_plot_next_sameseq_%s_n%s', ...
    basename, num2str(n_subjs)));

savefig(xfig, sprintf('%s.fig', fig_fname));
set(gcf,'PaperUnits','inches','PaperPosition', fig_rect/100)
print('-dpng', '-r300', sprintf('%s.png',fig_fname))
saveas(xfig, sprintf('%s.png',fig_fname), 'png')

close(xfig);

