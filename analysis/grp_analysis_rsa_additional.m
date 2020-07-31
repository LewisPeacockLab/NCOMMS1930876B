function[] = grp_analysis_rsa_additional(grp_rsa_add, args_rsa, dirs)

xph           = 1;
basename      = args_rsa.basename;
cate_color    = args_rsa.cate_color;
xparam        = args_rsa.index{xph}.param;
n_category    = xparam.n_category;
n_subcategory = xparam.n_subcategory;
category_name = xparam.category_name;
tsubj_grp     = args_rsa.filtered_subs;
rm_subs       = [1 5 21:23];
xsubj_grp     = tsubj_grp(~ismember(tsubj_grp, rm_subs));

n_subjs       = length(xsubj_grp);
n_runs        = xparam.n_runs;
n_targs       = 3;

%% ============= DECODING STUDY
clear xrsa_decode
xph = 2;

for xcate = 1:n_category
    for xsubcate = 1:n_subcategory
        for xtarg = 1:n_targs
            for xit = 1:n_runs
                for i = 1:2%1_rep, 2_run
                    xrsa_decode.cate{xcate}.subcate{xsubcate}.targ{xtarg}.run{i}{xit} = [];
                end
            end
        end
    end
end

for it = 1:length(xsubj_grp)
    xsub = xsubj_grp(it);
    
    for xcate = 1:n_category
        for xsubcate = 1:n_subcategory
            for xtarg = 1:n_targs
                
                xx = grp_rsa_add{xph}.subj{xsub}.decoding.cate{xcate}.subcate{xsubcate}.targ{xtarg};
                
                for xit = 1:length(xx.rep)
                    if ~isempty(xx.rep{xit}.corr_w)
                        tmp = isnan(xx.rep{xit}.corr_w);
                        if sum(tmp)
                            fprintf('%s, cate:%d, subcate:%d, targ:%d\n', num2str(xsub), xcate, xsubcate, xtarg);
                        else
                            xrsa_decode.cate{xcate}.subcate{xsubcate}.targ{xtarg}.run{1}{xit} = ...
                                horzcat(xrsa_decode.cate{xcate}.subcate{xsubcate}.targ{xtarg}.run{1}{xit}, ...
                                mean(xx.rep{xit}.corr_w));
                        end
                    end
                end
                
                
                for xit = 1:length(xx.run)
                    if ~isempty(xx.run{xit}.corr_w)
                        tmp = isnan(xx.run{xit}.corr_w);
                        if sum(tmp)
                            fprintf('%s, cate:%d, subcate:%d, targ:%d\n', num2str(xsub), xcate, xsubcate, xtarg);
                        else
                            xrsa_decode.cate{xcate}.subcate{xsubcate}.targ{xtarg}.run{2}{xit} = ...
                                horzcat(xrsa_decode.cate{xcate}.subcate{xsubcate}.targ{xtarg}.run{2}{xit}, ...
                                mean(xx.run{xit}.corr_w));
                        end
                    end
                end
            end
        end
    end
end

%%
for xcate = 1:n_category
    for xtarg = 1:n_targs
        for xit = 1:n_runs
            for i = 1:2%1_rep, 2_run
                t_corr = [];
                for xsubcate = 1:n_subcategory
                    t_corr = vertcat(t_corr, xrsa_decode.cate{xcate}.subcate{xsubcate}.targ{xtarg}.run{i}{xit});
                end
                
                xrsa_decode.cate{xcate}.targ{xtarg}.run{i}{xit} = mean(t_corr);
            end
        end
    end
end

%%

clear xcorr
for xcate = 1:n_category
    for xtarg = 1:n_targs
        for xit = 1:n_runs
            for i = 1:2%1_rep, 2_run
                xcorr{i}{xcate}.mean(xtarg, xit) = mean(xrsa_decode.cate{xcate}.targ{xtarg}.run{i}{xit});
                xcorr{i}{xcate}.se(xtarg, xit)   = std(xrsa_decode.cate{xcate}.targ{xtarg}.run{i}{xit})/sqrt(n_subjs); 
            end
        end
    end
end

%% ============= figure
xtarg = 1;
fig_rect = [0 0 1200 500];
y_lim    = [0 0.5];
xfig     = figure;
set(xfig, 'Position', fig_rect)

for i = 1:2%1_rep, 2_run
    clear xmean xse b
    
    subplot(1, 2, i);
    
    for xcate = 1:n_category
        xx = (1:5) + (5 * (xcate-1));
        b{xcate} = bar(xx, xcorr{i}{xcate}.mean(xtarg, :)); hold on %#ok<*AGROW>
        set(b{xcate}, 'facecolor', cate_color{xcate});
    end
    
    %*************** legend
    lg          = legend(category_name);
    lg.Location = 'BestOutside';%'SouthEast';
    lg.FontSize = 8;
    legend(category_name,'AutoUpdate','off')
    
    for xcate = 1:n_category
        xx = (1:5) + (5 * (xcate-1));
        errorbar(xx, xcorr{i}{xcate}.mean(xtarg, :), xcorr{i}{xcate}.se(xtarg, :), 'k.')
    end
    
    %*************** set
    set(gca,'xlim', [0 16],'ylim', y_lim);
    set(gca,'XTick', 1:(n_category * n_runs), 'XTickLabel', repmat(1:5, [1,3]), 'FontSize', 10);
  
    ylabel('RSA');
    
    if i==1    
        title(sprintf('RSA decoding: mean repetitions (N=%s)', num2str(n_subjs)));
    else
        title(sprintf('RSA decoding: run (N=%s)', num2str(n_subjs)));
    end
    
    %*************** stats
    for xcate = 1:n_category
        clear xp
        n = 0;
        for xit = 1:(n_runs-1)
            clear xx xx_sig
            xx = xrsa_decode.cate{xcate}.targ{xtarg}.run{i}{xit};
            
            for yit = xit+1:n_runs
                n = n+1;
                
                yy = xrsa_decode.cate{xcate}.targ{xtarg}.run{i}{yit};
                
                [~, xp(n)] = ttest(xx, yy);
            end     
        end
        
        [~, ~, ~, xadj_p] = fdr_bh(xp, 0.05, 'pdep', 'no');

        n = 0; sig_n = 0;
        for xit = 1:(n_runs-1)
            for yit = xit+1:n_runs
                n = n+1;
                
                if xadj_p(n) <= .05
                    sig_n = sig_n+1;
                    xx_sig(1) = xit + (n_runs * (xcate-1));
                    xx_sig(2) = yit + (n_runs * (xcate-1));
                    yy_sig    = y_lim(2) - (0.02*sig_n);
                    plot(xx_sig, [yy_sig yy_sig],'k-');
                    text(xx_sig(1) + (xx_sig(2)-xx_sig(1))/2, yy_sig, '*', 'FontSize', 10);
                end
            end
        end
    end
end

%*************** save fig
fig_fname = fullfile(dirs.rsa.group.parse, sprintf('RSA_decoding_n_repeat_%s', basename));

savefig(xfig, sprintf('%s.fig', fig_fname));
set(gcf,'PaperUnits','inches','PaperPosition', fig_rect/100)
print('-dpng', '-r300', sprintf('%s.png',fig_fname))
saveas(xfig, sprintf('%s.png',fig_fname), 'png')

close(xfig);

end
