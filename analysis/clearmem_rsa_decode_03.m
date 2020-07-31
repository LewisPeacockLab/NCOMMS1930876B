function[] = clearmem_rsa_decode_03(args, dirs)
%*************** RSA timecourse: test on study

%% ============= SETUP DIRECTORY
output_dir      = dirs.rsa.group.parse;
xsubj_grp       = args.filtered_subs;

%% ============= SETUP PARAMETERS
xparam          = args.index{1}.param;
n_category      = xparam.n_category;
category_names  = {'face','fruit','scene'};
conds_names     = {'maintain','replace category','replace subcategory','target suppress','global clear'};
n_condition     = length(conds_names);
n_tc_trs        = args.tc_tr_disp;
cate_members    = 1:n_category;

spm_p_thresh    = args.spm_p_thresh;
spm_v_extent    = args.spm_v_extent;
spm_thresh_type = args.spm_thresh_type;

rsa_mask        = args.rsa_mask;

basename = sprintf('patterns_%s_spmT_%s_thresh%s_ext%s_%s_mask', ...
    args.level, spm_thresh_type, num2str(spm_p_thresh), num2str(spm_v_extent), rsa_mask);

%% *************** plot parameter
xparam = args.index{2}.param;

for xcond = 1:n_condition
    
    if xcond~=2
        xtarg_color{xcond}{1} = [238, 20, 91]/255;% targ
        xtarg_color{xcond}{2} = [144, 144, 144]/255;% non_target
    else
        xtarg_color{xcond}{1} = [238, 20, 91]/255;% targ
        xtarg_color{xcond}{2} = [0, 188, 182]/255;% newtarg
        xtarg_color{xcond}{3} = [144, 144, 144]/255;% non_target
    end
end

xcond_color  = args.cond_color;
xcate_color  = args.cate_color;
xonset_color = args.onset_color;

dur_stim     = xparam.dur_stim;
dur_manip    = xparam.dur_manipulation;
dur_sync     = args.dur_sync_disp;

x_ticklable  = 1:n_tc_trs;

n_tr_blks    = args.n_tr_blks;%trs in one stats block: 5=2.3s/3=1.38s
n_stat_blks  = n_tc_trs/n_tr_blks;
xalpha       = 0.05/n_stat_blks;

%% ============= LOGGING
fprintf('(+) using scratch dir: %s\n', output_dir);
%*************** turn on diary to capture analysis output
diary off;
diary(sprintf('%s/diary.txt', output_dir));
fprintf('running code: %s at %s\n', mfilename, datestr(now, 0))
fprintf('#####################################################################\n\n');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%============= 1ST LEVEL RSA > GROUP STRUCTURE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
g_fname = fullfile(dirs.rsa.group.pattern, sprintf('group_%s.mat', basename));

if args.grp_pattern
    for xsub = args.g_sub
        
        %% ============= LOAD PATTERNS
        
        fprintf('\n(+) reseting grp_rsa patterns: %s\n', args.subject_list(xsub).name)
        
        fname = fullfile(dirs.rsa.group.pattern, sprintf('%s_%s.mat', basename, args.subject_list(xsub).name));
        load(fname);% 'patterns'
        
        %*************** LOCALIZER
        xph = 1;
        for xcate = 1:n_category
            
            fprintf('... %s: mask %s\n', args.phase_name{xph}, category_names{xcate});
            
            grp_rsa.subj{xsub}.patterns{xph}.maskcate{xcate}.beta = ...
                patterns{xph}.maskcate{xcate}.beta; %#ok<*USENS>
            
            for ycate = 1:3
                
                grp_rsa.subj{xsub}.patterns{xph}.maskcate{xcate}.targcate{ycate}.meanPat = ...
                    patterns{xph}.maskcate{xcate}.targcate{ycate}.meanPat;
                
                grp_rsa.subj{xsub}.patterns{xph}.maskcate{xcate}.targcate{ycate}.weighted.meanPat = ...
                    patterns{xph}.maskcate{xcate}.targcate{ycate}.weighted.meanPat;
                
                grp_rsa.subj{xsub}.patterns{xph}.maskcate{xcate}.targcate{ycate}.template = ...
                    patterns{xph}.maskcate{xcate}.targcate{ycate}.template;
                
                grp_rsa.subj{xsub}.patterns{xph}.maskcate{xcate}.targcate{ycate}.weighted.template = ...
                    patterns{xph}.maskcate{xcate}.targcate{ycate}.weighted.template;
            end
        end
        
        grp_rsa.subj{xsub}.patterns{xph}.info = patterns{xph}.info;
        
        %% *************** STUDY
        clear xbeta
        xph   = 2;
        
        %*************** beta
        for xcate = 1:n_category, xbeta{xcate} = patterns{1}.maskcate{xcate}.beta; end
        
        for xcate = 1:n_category
            
            fprintf('... %s: mask %s\n', args.phase_name{xph}, category_names{xcate});
            
            for ycate = 1:3
                
                grp_rsa.subj{xsub}.patterns{xph}.maskcate{xcate}.targcate{ycate}.meanPat = ...
                    patterns{xph}.maskcate{xcate}.targcate{ycate}.meanPat;
                
                grp_rsa.subj{xsub}.patterns{xph}.maskcate{xcate}.targcate{ycate}.template = ...
                    patterns{xph}.maskcate{xcate}.targcate{ycate}.template;
                
                n_runs = size(patterns{xph}.maskcate{xcate}.targcate{ycate}.meanPat,2);
                
                clear xpat_w
                for xrun = 1:n_runs
                    
                    %*************** weighted with beta
                    xpat           = patterns{xph}.maskcate{xcate}.targcate{ycate}.meanPat(:,xrun)';
                    ttpat          = xpat.*xbeta{xcate}'; %#ok<*AGROW>
                    xpat_w(:,xrun) = ttpat';
                end
                
                grp_rsa.subj{xsub}.patterns{xph}.maskcate{xcate}.targcate{ycate}.weighted.meanPat  = xpat_w;
                grp_rsa.subj{xsub}.patterns{xph}.maskcate{xcate}.targcate{ycate}.weighted.template = mean(xpat_w,2);
            end
        end
        
        %% *************** condition: timecourse
        for xcond = 1:n_condition
            
            fprintf('... %s timecourse: mask %s\n', args.phase_name{xph}, conds_names{xcond});
            
            if xcond~=2
                %*************** 1,3,4 condition
                for xcate = 1:n_category
                    for xtarg = 1:3
                        clear ycate xpat xpat_w
                        
                        noncate = cate_members(~ismember(cate_members, xcate));
                        
                        if xtarg==1,     ycate = xcate;
                        elseif xtarg==2, ycate = noncate(1);
                        elseif xtarg==3, ycate = noncate(2);
                        end
                        
                        it_n_trs = length(patterns{xph}.timecourse.cond{xcond}.targcate{xcate}.targ{xtarg}.mean_trPats);
                        
                        %*************** extract/weighting patterns
                        for xtr = 1:it_n_trs
                            tpat  = patterns{xph}.timecourse.cond{xcond}.targcate{xcate}.targ{xtarg}.mean_trPats{xtr};
                            ttpat = tpat';
                            
                            xpat(:, xtr)   = tpat;
                            xpat_w(:, xtr) = ttpat.*xbeta{ycate}';
                        end
                        
                        grp_rsa.subj{xsub}.patterns{xph}.timecourse.cond{xcond}.targcate{xcate}.targ{xtarg}.mean_trPats          = xpat;
                        grp_rsa.subj{xsub}.patterns{xph}.timecourse.cond{xcond}.targcate{xcate}.targ{xtarg}.weighted.mean_trPats = xpat_w;
                        
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        %*************** sync TRs
                        it_n_trs = length(patterns{xph}.timecourse.sync.cond{xcond}.targcate{xcate}.targ{xtarg}.mean_trPats);
                        
                        %*************** extract/weighting patterns
                        for xtr = 1:it_n_trs
                            tpat  = patterns{xph}.timecourse.sync.cond{xcond}.targcate{xcate}.targ{xtarg}.mean_trPats{xtr};
                            ttpat = tpat';
                            
                            xpat(:, xtr)   = tpat;
                            xpat_w(:, xtr) = ttpat.*xbeta{ycate}';
                        end
                        
                        grp_rsa.subj{xsub}.patterns{xph}.timecourse.sync.cond{xcond}.targcate{xcate}.targ{xtarg}.mean_trPats          = xpat;
                        grp_rsa.subj{xsub}.patterns{xph}.timecourse.sync.cond{xcond}.targcate{xcate}.targ{xtarg}.weighted.mean_trPats = xpat_w;
                    end
                end
            else
                %*************** replace category
                for xcate = 1:n_category
                    
                    it_newcate = cate_members(~ismember(cate_members, xcate));
                    
                    for xnewcate = it_newcate
                        for xtarg = 1:3
                            clear ycate xpat xpat_w
                            
                            noncate = cate_members(~ismember(cate_members, [xcate xnewcate]));
                            
                            if xtarg==1,     ycate = xcate;
                            elseif xtarg==2, ycate = xnewcate;
                            elseif xtarg==3, ycate = noncate;
                            end
                            
                            n_tc_trs = length(patterns{xph}.timecourse.cond{xcond}.targcate{xcate}.newtarg{xnewcate}.targ{xtarg}.mean_trPats);
                            
                            for xtr = 1:n_tc_trs
                                tpat  = patterns{xph}.timecourse.cond{xcond}.targcate{xcate}.newtarg{xnewcate}.targ{xtarg}.mean_trPats{xtr};
                                ttpat = tpat';
                                
                                xpat(:, xtr)   = tpat;
                                xpat_w(:, xtr) = ttpat.*xbeta{ycate}';
                            end
                            
                            grp_rsa.subj{xsub}.patterns{xph}.timecourse.cond{xcond}.targcate{xcate}.newtarg{xnewcate}.targ{xtarg}.mean_trPats          = xpat;
                            grp_rsa.subj{xsub}.patterns{xph}.timecourse.cond{xcond}.targcate{xcate}.newtarg{xnewcate}.targ{xtarg}.weighted.mean_trPats = xpat_w;
                            
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            %*************** sync TRs
                            n_tc_trs = length(patterns{xph}.timecourse.sync.cond{xcond}.targcate{xcate}.newtarg{xnewcate}.targ{xtarg}.mean_trPats);
                            
                            for xtr = 1:n_tc_trs
                                tpat  = patterns{xph}.timecourse.sync.cond{xcond}.targcate{xcate}.newtarg{xnewcate}.targ{xtarg}.mean_trPats{xtr};
                                ttpat = tpat';
                                
                                xpat(:, xtr)   = tpat;
                                xpat_w(:, xtr) = ttpat.*xbeta{ycate}';
                            end
                            
                            grp_rsa.subj{xsub}.patterns{xph}.timecourse.sync.cond{xcond}.targcate{xcate}.newtarg{xnewcate}.targ{xtarg}.mean_trPats          = xpat;
                            grp_rsa.subj{xsub}.patterns{xph}.timecourse.sync.cond{xcond}.targcate{xcate}.newtarg{xnewcate}.targ{xtarg}.weighted.mean_trPats = xpat_w;
                            
                        end
                    end
                end
            end
        end
        
        %% *************** info
        grp_rsa.subj{xsub}.patterns{xph}.info = patterns{xph}.info;
        
    end
    
    fprintf('\n....saving grp_rsa patterns\n')
    save(g_fname, 'grp_rsa','-v7.3')
    
    g_fsize = dir(g_fname);
    fprintf('file size: %s GB\n', num2str(g_fsize.bytes/(10^9)));
    
else
    
    fprintf('\n....loading existing grp_rsa patterns\n')
    load(g_fname);%grp_rsa.subj{xsub}
    
    g_fsize = dir(g_fname);
    fprintf('file size: %s GB\n', num2str(g_fsize.bytes/(10^9)));
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(args.cluster, 'local')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%============= EXTRACT PARAMS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    clear vox_nums
    
    xph      = 1;
    vox_nums = [];
    
    for xsub = xsubj_grp
        vox_nums = vertcat(vox_nums, grp_rsa.subj{xsub}.patterns{xph}.info.nVox);
    end
    
    n_selected_voxels.mean = mean(vox_nums);
    n_selected_voxels.std  = std(vox_nums);
    
    %*************** plot
    fig_rect = [0 0 500 400];
    
    xfig     = figure;
    set(xfig, 'Position', fig_rect)
    
    for xcate = 1:n_category
        h{xcate} = bar(xcate, n_selected_voxels.mean(xcate)); hold on
        set(h{xcate}, 'facecolor', xcate_color{xcate});
    end
    
    %*************** legend
    xlegend        = category_names;
    lg             = legend(xlegend);
    lg.Location    = 'BestOutside';
    legend(xlegend,'AutoUpdate','off')
    grid on
    
    for xcate = 1:n_category
        errorbar(xcate, n_selected_voxels.mean(xcate), n_selected_voxels.std(xcate),'k.')
        text(xcate, n_selected_voxels.mean(xcate) + 100, sprintf('%4.2f',n_selected_voxels.mean(xcate)));
    end
    
    set(gca,'XTick', 1:3)
    set(gca,'XTickLabel', category_names)
    
    title('# of voxels for ROI: uncorrected p<.05, extent vox = 10');
    ylabel('number of voxels');
    
    %%*************** save fig
    fig_fname = fullfile(output_dir, sprintf('plot_ROI_seleced_voxels_%s_n%s', ...
        basename, num2str(length(xsubj_grp))));
    
    savefig(xfig, sprintf('%s.fig', fig_fname));
    set(gcf,'PaperUnits','inches','PaperPosition', fig_rect/100)
    print('-djpeg', sprintf('%s.jpg',fig_fname), '-r100')
    saveas(xfig, sprintf('%s.jpg',fig_fname), 'jpg')
    
    close(xfig);
    
    %%*************** write table
    xtable = array2table(vox_nums, 'VariableNames', category_names);
    
    %*************** write tables to csv files
    csv_name = sprintf('%s/selected_voxs_%s_n%s.csv', output_dir, basename, num2str(length(xsubj_grp)));
    writetable(xtable, csv_name, 'WriteRowNames', true)
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%============= LOCALIZER SETUP STRUCTURES
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    clear n_runs
    
    xph         = 1;
    n_runs{xph} = 5;%args.index{xph}.param.n_runs;
    
    %*************** Pearson correlation
    % run based
    for xmaskcate = 1:n_category
        for xcate = 1:n_category
            for xrun = 1:n_runs{xph}
                for ycate = 1:n_category
                    for yrun = 1:n_runs{xph}
                        grp_rsa.loc.mask{xmaskcate}.run.corr.raw{xcate}{xrun}{ycate}{yrun}       = [];
                        grp_rsa.loc.mask{xmaskcate}.run.corr.weighted{xcate}{xrun}{ycate}{yrun}  = [];
                    end
                end
            end
        end
    end
    
    % template (mean of runs)
    for xmaskcate = 1:n_category
        for xcate = 1:n_category
            for ycate = 1:n_category
                grp_rsa.loc.mask{xmaskcate}.template.corr.raw{xcate}{ycate}      = [];
                grp_rsa.loc.mask{xmaskcate}.template.corr.weighted{xcate}{ycate} = [];
            end
        end
    end
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%============= LOCALIZER
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    fprintf('#####################################################################\n');
    fprintf('=========== RSA for within Localizer: %s \n', args.level);
    fprintf('#####################################################################\n');
    
    for xsub = xsubj_grp
        %*************** LOAD PATTERNS
        
        fprintf('%s ', num2str(xsub));
        
        xpat = grp_rsa.subj{xsub}.patterns;
        
        select_runs   = unique(xpat{xph}.info.matrix(findCol(xpat{xph}.info.header, {'it_run'}),:));
        
        %% ============= CONFUSION MATRIX WITHIN LOCALIZER: PER MASK: PER RUN
        %*************** verification of RSA
        for xmaskcate = 1:n_category
            for xcate = 1:n_category
                for xrun = 1:size(xpat{xph}.maskcate{xmaskcate}.targcate{xcate}.meanPat,2)
                    
                    xit_run = select_runs(xrun);
                    
                    %*************** reference pattern
                    xref_pat   = xpat{xph}.maskcate{xmaskcate}.targcate{xcate}.meanPat(:,xrun); %#ok<*NODEF>
                    xref_pat_w = xpat{xph}.maskcate{xmaskcate}.targcate{xcate}.weighted.meanPat(:,xrun);
                    
                    for ycate = 1:n_category
                        for yrun = 1:size(xpat{xph}.maskcate{xmaskcate}.targcate{ycate}.meanPat, 2)
                            
                            yit_run = select_runs(yrun);
                            
                            %*************** test pattern
                            y_pat   = xpat{xph}.maskcate{xmaskcate}.targcate{ycate}.meanPat(:,yrun);
                            y_pat_w = xpat{xph}.maskcate{xmaskcate}.targcate{ycate}.weighted.meanPat(:,yrun);
                            
                            %*************** similarity: z-transformed pearson correlation
                            xr = corr2(xref_pat, y_pat); xwr = corr2(xref_pat_w, y_pat_w);
                            
                            grp_rsa.loc.mask{xmaskcate}.run.corr.raw{xcate}{xit_run}{ycate}{yit_run} = ...
                                horzcat(grp_rsa.loc.mask{xmaskcate}.run.corr.raw{xcate}{xit_run}{ycate}{yit_run}, ...
                                .5*(log(1+xr) - log(1-xr)));
                            
                            grp_rsa.loc.mask{xmaskcate}.run.corr.weighted{xcate}{xit_run}{ycate}{yit_run} = ...
                                horzcat(grp_rsa.loc.mask{xmaskcate}.run.corr.weighted{xcate}{xit_run}{ycate}{yit_run}, ...
                                .5*(log(1+xwr) - log(1-xwr)));
                        end
                    end
                end
            end
        end
        
        %% ============= CONFUSION MATRIX WITHIN LOCALIZER: MEAN
        for xmaskcate = 1:n_category
            for xcate = 1:n_category
                %*************** reference pattern
                xref_pat   = xpat{xph}.maskcate{xmaskcate}.targcate{xcate}.template;
                xref_pat_w = xpat{xph}.maskcate{xmaskcate}.targcate{xcate}.weighted.template;
                
                for ycate = 1:n_category
                    %*************** test pattern
                    y_pat   = xpat{xph}.maskcate{xmaskcate}.targcate{ycate}.template;
                    y_pat_w = xpat{xph}.maskcate{xmaskcate}.targcate{ycate}.weighted.template;
                    
                    %*************** similarity: z-transformed pearson correlation
                    xr = corr2(xref_pat, y_pat); xwr = corr2(xref_pat_w, y_pat_w);
                    
                    grp_rsa.loc.mask{xmaskcate}.template.corr.raw{xcate}{ycate} = ...
                        horzcat(grp_rsa.loc.mask{xmaskcate}.template.corr.raw{xcate}{ycate}, ...
                        .5*(log(1+xr) - log(1-xr)));
                    
                    grp_rsa.loc.mask{xmaskcate}.template.corr.weighted{xcate}{ycate} = ...
                        horzcat(grp_rsa.loc.mask{xmaskcate}.template.corr.weighted{xcate}{ycate}, ...
                        .5*(log(1+xwr) - log(1-xwr)));
                end
            end
        end
    end
    
    fprintf('\n');
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%============= LOCALIZER PLOT
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%============= CONFUSION MATRIX WITHIN LOCALIZER: PER MASK: PER RUN
    xph = 1;
    
    for xmaskcate = 1:n_category
        
        grp_rsa.loc.mask{xmaskcate}.run.corr.heatmap.raw      = zeros(n_category*n_runs{xph},n_category*n_runs{xph});
        grp_rsa.loc.mask{xmaskcate}.run.corr.heatmap.weighted = zeros(n_category*n_runs{xph},n_category*n_runs{xph});
        
        for xcate = 1:n_category
            for xrun = 1:n_runs{xph}
                
                xunit = xrun + (n_runs{1} * (xcate-1));% targ
                
                for ycate = 1:n_category
                    for yrun = 1:n_runs{xph}
                        
                        yunit = yrun + (n_runs{1} * (ycate-1));% decoded
                        
                        %*************** similarity: pearson correlation
                        grp_rsa.loc.mask{xmaskcate}.run.corr.heatmap.raw(xunit,yunit) =...
                            mean(grp_rsa.loc.mask{xmaskcate}.run.corr.raw{xcate}{xrun}{ycate}{yrun});
                        
                        grp_rsa.loc.mask{xmaskcate}.run.corr.heatmap.weighted(xunit,yunit) =...
                            mean(grp_rsa.loc.mask{xmaskcate}.run.corr.weighted{xcate}{xrun}{ycate}{yrun});
                    end
                end
            end
        end
    end
    
    %% --------- figure
    fig_rect    = [0 0 2000 1200];
    xticks      = round(n_runs{xph}/2):n_runs{xph}:(n_runs{xph} * n_category);
    bar_lim     = [-0.9 1.5];%0.9
    bar_lim_w   = [-0.9 1.5];%0.7
    
    heatmap_fig = figure;
    set(heatmap_fig, 'Position', fig_rect)
    
    for xmaskcate = 1:n_category
        
        %*************** similarity: pearson correlation
        subplot(2, 3, xmaskcate)
        imagesc(grp_rsa.loc.mask{xmaskcate}.run.corr.heatmap.raw);
        colormap('jet'); colorbar; caxis(bar_lim); hold on
        
        title(sprintf('Localizer RSA: ROI for %s (N=%s)', ...
            category_names{xmaskcate}, num2str(length(xsubj_grp))),...
            'FontSize',15,'FontWeight','bold');
        xlabel('decoded category', 'FontSize', 15);
        ylabel('target category', 'FontSize', 15);
        
        set(gca,'YTick', xticks, 'YTickLabel', category_names, 'FontSize', 12);
        set(gca,'XTick', xticks, 'XTickLabel', category_names, 'FontSize', 12);
        set(gca,'xaxisLocation','top')
        
        %*************** similarity: pearson correlation: weighted
        subplot(2, 3, xmaskcate + 3)
        imagesc(grp_rsa.loc.mask{xmaskcate}.run.corr.heatmap.weighted);
        colormap('jet'); colorbar; caxis(bar_lim_w); hold on
        
        title(sprintf('Localizer weighted RSA: ROI for %s (N=%s)', ...
            category_names{xmaskcate}, num2str(length(xsubj_grp))),...
            'FontSize',15,'FontWeight','bold');
        xlabel('target category', 'FontSize', 15);
        ylabel('decode category', 'FontSize', 15);
        
        set(gca,'YTick', xticks, 'YTickLabel', category_names, 'FontSize', 12);
        set(gca,'XTick', xticks, 'XTickLabel', category_names, 'FontSize', 12);
        set(gca,'xaxisLocation','top')
        
    end
    
    %*************** save fig
    fig_fname = fullfile(output_dir, sprintf('plot_similarity_matrix_%s_run_%s_n%s', ...
        args.phase_name{xph}, basename, num2str(length(xsubj_grp))));
    
    savefig(heatmap_fig, sprintf('%s.fig', fig_fname));
    set(gcf,'PaperUnits','inches','PaperPosition', fig_rect/100)
    print('-djpeg', sprintf('%s.jpg',fig_fname), '-r100')
    saveas(heatmap_fig, sprintf('%s.jpg',fig_fname), 'jpg')
    
    close(heatmap_fig);
    
    %% ============= CONFUSION MATRIX WITHIN LOCALIZER: MEAN
    for xmaskcate = 1:n_category
        
        grp_rsa.loc.mask{xmaskcate}.template.corr.heatmap.raw      = zeros(n_category, n_category);
        grp_rsa.loc.mask{xmaskcate}.template.corr.heatmap.weighted = zeros(n_category, n_category);
        
        for xcate = 1:n_category
            for ycate = 1:n_category
                
                %*************** similarity: pearson correlation
                grp_rsa.loc.mask{xmaskcate}.template.corr.heatmap.raw(xcate,ycate) =...
                    mean(grp_rsa.loc.mask{xmaskcate}.template.corr.raw{xcate}{ycate});
                
                grp_rsa.loc.mask{xmaskcate}.template.corr.heatmap.weighted(xcate,ycate) =...
                    mean(grp_rsa.loc.mask{xmaskcate}.template.corr.weighted{xcate}{ycate});
            end
        end
    end
    
    %% --------- figure
    fig_rect    = [0 0 2000 1200];
    xticks      = 1:n_category;
    bar_lim     = [-0.9 3];
    
    heatmap_fig = figure;
    set(heatmap_fig, 'Position', fig_rect)
    
    for xmaskcate = 1:n_category
        
        %*************** similarity: pearson correlation
        subplot(2, 3, xmaskcate)
        imagesc(grp_rsa.loc.mask{xmaskcate}.template.corr.heatmap.raw);
        colormap('jet'); colorbar; caxis(bar_lim); hold on
        
        %*************** z-transformed r
        for xrow = 1:3
            for xcol = 1:3
                text(xrow, xcol,...
                    sprintf('%4.4f', grp_rsa.loc.mask{xmaskcate}.template.corr.heatmap.raw(xrow, xcol)),...
                    'FontSize', 15)
            end
        end
        
        title(sprintf('Localizer RSA: ROI for %s (N=%s)', ...
            category_names{xmaskcate}, num2str(length(xsubj_grp))),...
            'FontSize',15,'FontWeight','bold');
        xlabel('target category', 'FontSize', 15);
        ylabel('decode category', 'FontSize', 15);
        
        set(gca,'YTick', xticks, 'YTickLabel', category_names, 'FontSize', 12);
        set(gca,'XTick', xticks, 'XTickLabel', category_names, 'FontSize', 12);
        set(gca,'xaxisLocation','top')
        
        %*************** similarity: pearson correlation: weighted
        subplot(2, 3, xmaskcate + 3)
        imagesc(grp_rsa.loc.mask{xmaskcate}.template.corr.heatmap.weighted);
        colormap('jet'); colorbar; caxis(bar_lim); hold on
        
        %*************** z-transformed r
        for xrow = 1:3
            for xcol = 1:3
                text(xrow, xcol,...
                    sprintf('%4.4f', grp_rsa.loc.mask{xmaskcate}.template.corr.heatmap.weighted(xrow, xcol)),...
                    'FontSize', 15)
            end
        end
        
        title(sprintf('Localizer weighted RSA: ROI for %s (N=%s)', ...
            category_names{xmaskcate}, num2str(length(xsubj_grp))),...
            'FontSize',15,'FontWeight','bold');
        xlabel('target category', 'FontSize', 15);
        ylabel('decode category', 'FontSize', 15);
        
        set(gca,'YTick', xticks, 'YTickLabel', category_names, 'FontSize', 12);
        set(gca,'XTick', xticks, 'XTickLabel', category_names, 'FontSize', 12);
        set(gca,'xaxisLocation','top')
        
    end
    
    %*************** save fig
    fig_fname = fullfile(output_dir, sprintf('plot_similarity_matrix_%s_template_%s_n%s', ...
        args.phase_name{xph}, basename, num2str(length(xsubj_grp))));
    
    savefig(heatmap_fig, sprintf('%s.fig', fig_fname));
    set(gcf,'PaperUnits','inches','PaperPosition', fig_rect/100)
    print('-djpeg', sprintf('%s.jpg',fig_fname), '-r100')
    saveas(heatmap_fig, sprintf('%s.jpg',fig_fname), 'jpg')
    
    close(heatmap_fig);
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%============= STUDY SETUP STRUCTURES
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    xph = 2;
    
    %*************** Pearson correlation
    % xtarg: 1_target, 2_newtarg/nontarg, 3_nontarg
    % timecourse
    for xcond = 1:n_condition
        if xcond~=2
            for xcate = 1:n_category
                for xtarg = 1:3
                    for xtr = 1:n_tc_trs
                        grp_rsa.study.timecourse.cond{xcond}.targcate{xcate}.targ{xtarg}.corr.raw{xtr}      = [];
                        grp_rsa.study.timecourse.cond{xcond}.targcate{xcate}.targ{xtarg}.corr.weighted{xtr} = [];
                        grp_rsa.study.timecourse.cond{xcond}.targcate{xcate}.targ{xtarg}.corr.beta{xtr}     = [];
                    end
                    
                    for xtr = 1:dur_sync
                        grp_rsa.study.timecourse.sync.cond{xcond}.targcate{xcate}.targ{xtarg}.corr.raw{xtr}      = [];
                        grp_rsa.study.timecourse.sync.cond{xcond}.targcate{xcate}.targ{xtarg}.corr.weighted{xtr} = [];
                        grp_rsa.study.timecourse.sync.cond{xcond}.targcate{xcate}.targ{xtarg}.corr.beta{xtr}     = [];
                    end
                end
            end
        else
            for xcate = 1:n_category
                for xnewcate = cate_members(~ismember(cate_members, xcate))
                    for xtarg = 1:3
                        for xtr = 1:n_tc_trs
                            grp_rsa.study.timecourse.cond{xcond}.targcate{xcate}.newtarg{xnewcate}.targ{xtarg}.corr.raw{xtr}      = [];
                            grp_rsa.study.timecourse.cond{xcond}.targcate{xcate}.newtarg{xnewcate}.targ{xtarg}.corr.weighted{xtr} = [];
                            grp_rsa.study.timecourse.cond{xcond}.targcate{xcate}.newtarg{xnewcate}.targ{xtarg}.corr.beta{xtr}     = [];
                        end
                        
                        for xtr = 1:dur_sync
                            grp_rsa.study.timecourse.sync.cond{xcond}.targcate{xcate}.newtarg{xnewcate}.targ{xtarg}.corr.raw{xtr}      = [];
                            grp_rsa.study.timecourse.sync.cond{xcond}.targcate{xcate}.newtarg{xnewcate}.targ{xtarg}.corr.weighted{xtr} = [];
                            grp_rsa.study.timecourse.sync.cond{xcond}.targcate{xcate}.newtarg{xnewcate}.targ{xtarg}.corr.beta{xtr}     = [];
                        end
                    end
                end
            end
        end
    end
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%============= STUDY
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%============= RSA
    
    fprintf('#####################################################################\n');
    fprintf('=========== TIMECOURSE RSA for Study: %s \n', args.level);
    fprintf('#####################################################################\n');
    
    for xsub = xsubj_grp
        %*************** LOAD PATTERNS
        
        fprintf('%s ', num2str(xsub));
        
        xpat = grp_rsa.subj{xsub}.patterns;
        
        %% ============= TIMECOURSE PER CATEGORY
        % patterns{xph}.timecourse.cond{xcond}.target{xtarg}.cate{xcate}.trPats{xtr}
        % xtarg: 1_target, 2_newtarg/nontarg, 3_nontarg
        
        for xcond = 1:n_condition
            for xcate = 1:n_category
                clear xnoncate xnewcate xref_pat xref_pat_w
                %*************** targ: reference pattern
                xref_pat{1}   = xpat{1}.maskcate{xcate}.targcate{xcate}.template;
                xref_pat_w{1} = xpat{1}.maskcate{xcate}.targcate{xcate}.weighted.template;
                xref_pat_b{1} = xpat{1}.maskcate{xcate}.beta;
                
                if xcond ~=2
                    %*************** non target category
                    xnoncate = cate_members(~ismember(cate_members, xcate));
                    
                    %*************** nontarg: reference pattern
                    for xx = 1:length(xnoncate)
                        xref_pat{xx+1}   = xpat{1}.maskcate{xnoncate(xx)}.targcate{xnoncate(xx)}.template;
                        xref_pat_w{xx+1} = xpat{1}.maskcate{xnoncate(xx)}.targcate{xnoncate(xx)}.weighted.template;
                        xref_pat_b{xx+1} = xpat{1}.maskcate{xnoncate(xx)}.beta;
                    end
                    
                    for xtarg = 1:3
                        for xtr = 1:n_tc_trs
                            clear xr xwr
                            %*************** test pattern
                            y_pat   = xpat{xph}.timecourse.cond{xcond}.targcate{xcate}.targ{xtarg}.mean_trPats(:, xtr);
                            y_pat_w = xpat{xph}.timecourse.cond{xcond}.targcate{xcate}.targ{xtarg}.weighted.mean_trPats(:, xtr);
                            
                            %*************** similarity: z-transformed pearson correlation
                            xr  = corr2(xref_pat{xtarg}, y_pat);
                            xwr = corr2(xref_pat_w{xtarg}, y_pat_w);
                            xbr = corr2(xref_pat_b{xtarg}, y_pat_w);
                            
                            grp_rsa.study.timecourse.cond{xcond}.targcate{xcate}.targ{xtarg}.corr.raw{xtr} = ...
                                horzcat(grp_rsa.study.timecourse.cond{xcond}.targcate{xcate}.targ{xtarg}.corr.raw{xtr}, ...
                                .5*(log(1+xr) - log(1-xr)));
                            
                            grp_rsa.study.timecourse.cond{xcond}.targcate{xcate}.targ{xtarg}.corr.weighted{xtr} = ...
                                horzcat(grp_rsa.study.timecourse.cond{xcond}.targcate{xcate}.targ{xtarg}.corr.weighted{xtr}, ...
                                .5*(log(1+xwr) - log(1-xwr)));
                            
                            grp_rsa.study.timecourse.cond{xcond}.targcate{xcate}.targ{xtarg}.corr.beta{xtr} = ...
                                horzcat(grp_rsa.study.timecourse.cond{xcond}.targcate{xcate}.targ{xtarg}.corr.beta{xtr}, ...
                                .5*(log(1+xbr) - log(1-xbr)));
                        end
                        
                        %*************** sync trs
                        for xtr = 1:dur_sync
                            clear xr xwr
                            %*************** test pattern
                            y_pat   = xpat{xph}.timecourse.sync.cond{xcond}.targcate{xcate}.targ{xtarg}.mean_trPats(:, xtr);
                            y_pat_w = xpat{xph}.timecourse.sync.cond{xcond}.targcate{xcate}.targ{xtarg}.weighted.mean_trPats(:, xtr);
                            
                            %*************** similarity: z-transformed pearson correlation
                            xr  = corr2(xref_pat{xtarg}, y_pat);
                            xwr = corr2(xref_pat_w{xtarg}, y_pat_w);
                            xbr = corr2(xref_pat_b{xtarg}, y_pat_w);
                            
                            grp_rsa.study.timecourse.sync.cond{xcond}.targcate{xcate}.targ{xtarg}.corr.raw{xtr} = ...
                                horzcat(grp_rsa.study.timecourse.sync.cond{xcond}.targcate{xcate}.targ{xtarg}.corr.raw{xtr}, ...
                                .5*(log(1+xr) - log(1-xr)));
                            
                            grp_rsa.study.timecourse.sync.cond{xcond}.targcate{xcate}.targ{xtarg}.corr.weighted{xtr} = ...
                                horzcat(grp_rsa.study.timecourse.sync.cond{xcond}.targcate{xcate}.targ{xtarg}.corr.weighted{xtr}, ...
                                .5*(log(1+xwr) - log(1-xwr)));
                            
                            grp_rsa.study.timecourse.sync.cond{xcond}.targcate{xcate}.targ{xtarg}.corr.beta{xtr} = ...
                                horzcat(grp_rsa.study.timecourse.sync.cond{xcond}.targcate{xcate}.targ{xtarg}.corr.beta{xtr}, ...
                                .5*(log(1+xbr) - log(1-xbr)));
                            
                        end
                    end
                else
                    %*************** new category
                    xnewcate = cate_members(~ismember(cate_members, xcate));
                    
                    for xx = 1:length(xnewcate)
                        %*************** non target category
                        xnoncate = cate_members(~ismember(cate_members, [xcate xnewcate(xx)]));
                        
                        %*************** newtarg: reference pattern
                        xref_pat{2}   = xpat{1}.maskcate{xnewcate(xx)}.targcate{xnewcate(xx)}.template;
                        xref_pat_w{2} = xpat{1}.maskcate{xnewcate(xx)}.targcate{xnewcate(xx)}.weighted.template;
                        xref_pat_b{2} = xpat{1}.maskcate{xnewcate(xx)}.beta;
                        
                        %*************** nontarg: reference pattern
                        xref_pat{3}   = xpat{1}.maskcate{xnoncate}.targcate{xnoncate}.template;
                        xref_pat_w{3} = xpat{1}.maskcate{xnoncate}.targcate{xnoncate}.weighted.template;
                        xref_pat_b{3} = xpat{1}.maskcate{xnoncate}.beta;
                        
                        for xtarg = 1:3
                            for xtr = 1:n_tc_trs
                                clear xr xwr
                                %*************** test pattern
                                y_pat   = xpat{xph}.timecourse.cond{xcond}.targcate{xcate}.newtarg{xnewcate(xx)}.targ{xtarg}.mean_trPats(:, xtr);
                                y_pat_w = xpat{xph}.timecourse.cond{xcond}.targcate{xcate}.newtarg{xnewcate(xx)}.targ{xtarg}.weighted.mean_trPats(:, xtr);
                                
                                %*************** similarity: z-transformed pearson correlation
                                xr  = corr2(xref_pat{xtarg}, y_pat);
                                xwr = corr2(xref_pat_w{xtarg}, y_pat_w);
                                xbr = corr2(xref_pat_b{xtarg}, y_pat_w);
                                
                                grp_rsa.study.timecourse.cond{xcond}.targcate{xcate}.newtarg{xnewcate(xx)}.targ{xtarg}.corr.raw{xtr} = ...
                                    horzcat(grp_rsa.study.timecourse.cond{xcond}.targcate{xcate}.newtarg{xnewcate(xx)}.targ{xtarg}.corr.raw{xtr}, ...
                                    .5*(log(1+xr) - log(1-xr)));
                                
                                grp_rsa.study.timecourse.cond{xcond}.targcate{xcate}.newtarg{xnewcate(xx)}.targ{xtarg}.corr.weighted{xtr} = ...
                                    horzcat(grp_rsa.study.timecourse.cond{xcond}.targcate{xcate}.newtarg{xnewcate(xx)}.targ{xtarg}.corr.weighted{xtr}, ...
                                    .5*(log(1+xwr) - log(1-xwr)));
                                
                                grp_rsa.study.timecourse.cond{xcond}.targcate{xcate}.newtarg{xnewcate(xx)}.targ{xtarg}.corr.beta{xtr} = ...
                                    horzcat(grp_rsa.study.timecourse.cond{xcond}.targcate{xcate}.newtarg{xnewcate(xx)}.targ{xtarg}.corr.beta{xtr}, ...
                                    .5*(log(1+xbr) - log(1-xbr)));
                            end
                            
                            %*************** sync trs
                            for xtr = 1:dur_sync
                                clear xr xwr
                                %*************** test pattern
                                y_pat   = xpat{xph}.timecourse.sync.cond{xcond}.targcate{xcate}.newtarg{xnewcate(xx)}.targ{xtarg}.mean_trPats(:, xtr);
                                y_pat_w = xpat{xph}.timecourse.sync.cond{xcond}.targcate{xcate}.newtarg{xnewcate(xx)}.targ{xtarg}.weighted.mean_trPats(:, xtr);
                                
                                %*************** similarity: z-transformed pearson correlation
                                xr  = corr2(xref_pat{xtarg}, y_pat);
                                xwr = corr2(xref_pat_w{xtarg}, y_pat_w);
                                xbr = corr2(xref_pat_b{xtarg}, y_pat_w);
                                
                                grp_rsa.study.timecourse.sync.cond{xcond}.targcate{xcate}.newtarg{xnewcate(xx)}.targ{xtarg}.corr.raw{xtr} = ...
                                    horzcat(grp_rsa.study.timecourse.sync.cond{xcond}.targcate{xcate}.newtarg{xnewcate(xx)}.targ{xtarg}.corr.raw{xtr}, ...
                                    .5*(log(1+xr) - log(1-xr)));
                                
                                grp_rsa.study.timecourse.sync.cond{xcond}.targcate{xcate}.newtarg{xnewcate(xx)}.targ{xtarg}.corr.weighted{xtr} = ...
                                    horzcat(grp_rsa.study.timecourse.sync.cond{xcond}.targcate{xcate}.newtarg{xnewcate(xx)}.targ{xtarg}.corr.weighted{xtr}, ...
                                    .5*(log(1+xwr) - log(1-xwr)));
                                
                                grp_rsa.study.timecourse.sync.cond{xcond}.targcate{xcate}.newtarg{xnewcate(xx)}.targ{xtarg}.corr.beta{xtr} = ...
                                    horzcat(grp_rsa.study.timecourse.sync.cond{xcond}.targcate{xcate}.newtarg{xnewcate(xx)}.targ{xtarg}.corr.beta{xtr}, ...
                                    .5*(log(1+xbr) - log(1-xbr)));
                            end
                        end
                    end
                end
            end
        end
    end
    
    fprintf('\n')
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%============= TIMECOURSE: PER TARGET
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %*************** target/ non-target
    % grp_rsa.study.timecourse.cond{xcond}.corr.mean.raw(xtarg, xtr)
    % cond~=2: targ:1, non-targ:2
    % cond==2: targ:1, new_targ:2, non-targ:3
    for xcond = 1:n_condition
        
        if xcond == 2, n_targ = 3; else, n_targ = 2; end
        
        grp_rsa.study.targcorr.cond{xcond}.mean   = zeros(n_targ, n_tc_trs);
        grp_rsa.study.targcorr.cond{xcond}.se     = zeros(n_targ, n_tc_trs);
        grp_rsa.study.targcorr_w.cond{xcond}.mean = zeros(n_targ, n_tc_trs);
        grp_rsa.study.targcorr_w.cond{xcond}.se   = zeros(n_targ, n_tc_trs);
        grp_rsa.study.targcorr_b.cond{xcond}.mean = zeros(n_targ, n_tc_trs);
        grp_rsa.study.targcorr_b.cond{xcond}.se   = zeros(n_targ, n_tc_trs);
        
        for xtr = 1:n_tc_trs
            if xcond~=2
                %*************** 1 target: corr b/w target template & target trial (study)
                % e.g.: face_temp & face trial under face mask
                xtarg = 1;
                xtarg_corr = []; xtarg_corr_w = []; xtarg_corr_b = [];
                
                for xcate = 1:n_category
                    xtarg_corr = vertcat(xtarg_corr, ...
                        grp_rsa.study.timecourse.cond{xcond}.targcate{xcate}.targ{xtarg}.corr.raw{xtr});
                    xtarg_corr_w = vertcat(xtarg_corr_w, ...
                        grp_rsa.study.timecourse.cond{xcond}.targcate{xcate}.targ{xtarg}.corr.weighted{xtr});
                    xtarg_corr_b = vertcat(xtarg_corr_b, ...
                        grp_rsa.study.timecourse.cond{xcond}.targcate{xcate}.targ{xtarg}.corr.beta{xtr});
                end
                
                grp_rsa.study.targcorr.cond{xcond}.targ{xtarg}.tr{xtr}   = mean(xtarg_corr);
                grp_rsa.study.targcorr_w.cond{xcond}.targ{xtarg}.tr{xtr} = mean(xtarg_corr_w);
                grp_rsa.study.targcorr_b.cond{xcond}.targ{xtarg}.tr{xtr} = mean(xtarg_corr_b);
                
                grp_rsa.study.targcorr.cond{xcond}.mean(xtarg, xtr)      = mean(mean(xtarg_corr));
                grp_rsa.study.targcorr.cond{xcond}.se(xtarg, xtr)        = std(mean(xtarg_corr))/sqrt(length(xsubj_grp));
                
                grp_rsa.study.targcorr_w.cond{xcond}.mean(xtarg, xtr)    = mean(mean(xtarg_corr_w));
                grp_rsa.study.targcorr_w.cond{xcond}.se(xtarg, xtr)      = std(mean(xtarg_corr_w))/sqrt(length(xsubj_grp));
                
                grp_rsa.study.targcorr_b.cond{xcond}.mean(xtarg, xtr)    = mean(mean(xtarg_corr_b));
                grp_rsa.study.targcorr_b.cond{xcond}.se(xtarg, xtr)      = std(mean(xtarg_corr_b))/sqrt(length(xsubj_grp));
                
                %*************** 2 nontarget
                xtarg_corr = []; xtarg_corr_w = []; xtarg_corr_b = [];
                
                for xtarg = 2:3
                    for xcate = 1:n_category
                        xtarg_corr = vertcat(xtarg_corr, ...
                            grp_rsa.study.timecourse.cond{xcond}.targcate{xcate}.targ{xtarg}.corr.raw{xtr});
                        xtarg_corr_w = vertcat(xtarg_corr_w, ...
                            grp_rsa.study.timecourse.cond{xcond}.targcate{xcate}.targ{xtarg}.corr.weighted{xtr});
                        xtarg_corr_b = vertcat(xtarg_corr_b, ...
                            grp_rsa.study.timecourse.cond{xcond}.targcate{xcate}.targ{xtarg}.corr.beta{xtr});
                    end
                end
                
                xtarg = 2;
                grp_rsa.study.targcorr.cond{xcond}.targ{xtarg}.tr{xtr}   = mean(xtarg_corr);
                grp_rsa.study.targcorr_w.cond{xcond}.targ{xtarg}.tr{xtr} = mean(xtarg_corr_w);
                grp_rsa.study.targcorr_b.cond{xcond}.targ{xtarg}.tr{xtr} = mean(xtarg_corr_b);
                
                grp_rsa.study.targcorr.cond{xcond}.mean(xtarg, xtr)      = mean(mean(xtarg_corr));
                grp_rsa.study.targcorr.cond{xcond}.se(xtarg, xtr)        = std(mean(xtarg_corr))/sqrt(length(xsubj_grp));
                
                grp_rsa.study.targcorr_w.cond{xcond}.mean(xtarg, xtr)    = mean(mean(xtarg_corr_w));
                grp_rsa.study.targcorr_w.cond{xcond}.se(xtarg, xtr)      = std(mean(xtarg_corr_w))/sqrt(length(xsubj_grp));
                
                grp_rsa.study.targcorr_b.cond{xcond}.mean(xtarg, xtr)    = mean(mean(xtarg_corr_b));
                grp_rsa.study.targcorr_b.cond{xcond}.se(xtarg, xtr)      = std(mean(xtarg_corr_b))/sqrt(length(xsubj_grp));
                
            else
                for xtarg = 1:3
                    %*************** 2 newtarg, 3 nontarget
                    xtarg_corr = []; xtarg_corr_w = []; xtarg_corr_b = [];
                    
                    for xcate = 1:n_category
                        for xnewcate = cate_members(~ismember(cate_members, xcate))
                            xtarg_corr = vertcat(xtarg_corr, ...
                                grp_rsa.study.timecourse.cond{xcond}.targcate{xcate}.newtarg{xnewcate}.targ{xtarg}.corr.raw{xtr});
                            xtarg_corr_w = vertcat(xtarg_corr_w, ...
                                grp_rsa.study.timecourse.cond{xcond}.targcate{xcate}.newtarg{xnewcate}.targ{xtarg}.corr.weighted{xtr});
                            xtarg_corr_b = vertcat(xtarg_corr_b, ...
                                grp_rsa.study.timecourse.cond{xcond}.targcate{xcate}.newtarg{xnewcate}.targ{xtarg}.corr.beta{xtr});
                        end
                    end
                    
                    grp_rsa.study.targcorr.cond{xcond}.targ{xtarg}.tr{xtr}   = mean(xtarg_corr);
                    grp_rsa.study.targcorr_w.cond{xcond}.targ{xtarg}.tr{xtr} = mean(xtarg_corr_w);
                    grp_rsa.study.targcorr_b.cond{xcond}.targ{xtarg}.tr{xtr} = mean(xtarg_corr_b);
                    
                    grp_rsa.study.targcorr.cond{xcond}.mean(xtarg, xtr)      = mean(mean(xtarg_corr));
                    grp_rsa.study.targcorr.cond{xcond}.se(xtarg, xtr)        = std(mean(xtarg_corr))/sqrt(length(xsubj_grp));
                    
                    grp_rsa.study.targcorr_w.cond{xcond}.mean(xtarg, xtr)    = mean(mean(xtarg_corr_w));
                    grp_rsa.study.targcorr_w.cond{xcond}.se(xtarg, xtr)      = std(mean(xtarg_corr_w))/sqrt(length(xsubj_grp));
                    
                    grp_rsa.study.targcorr_b.cond{xcond}.mean(xtarg, xtr)    = mean(mean(xtarg_corr_b));
                    grp_rsa.study.targcorr_b.cond{xcond}.se(xtarg, xtr)      = std(mean(xtarg_corr_b))/sqrt(length(xsubj_grp));
                    
                end
            end
        end
    end
    
    %% ============= STATS: RANDOM
    for xcond = 1:n_condition
        for xblk = 1:n_stat_blks
            
            if xcond == 2, n_targ = 3; else, n_targ = 2; end
            
            grp_rsa.study.targcorr.ttest.cond{xcond}.blk{xblk}.pvalue    = nan(n_targ, n_targ);
            grp_rsa.study.targcorr_w.ttest.cond{xcond}.xblk{xblk}.pvalue = nan(n_targ, n_targ);
            grp_rsa.study.targcorr_b.ttest.cond{xcond}.xblk{xblk}.pvalue = nan(n_targ, n_targ);
            
            for xcol = 1:(n_targ-1)
                for xrow = (xcol+1):n_targ
                    
                    it_trs = (1:n_tr_blks) + (n_tr_blks * (xblk-1));
                    
                    for xx = 1:3
                        clear xpat
                        if xx == 1,     xpat = grp_rsa.study.targcorr.cond{xcond};
                        elseif xx == 2, xpat = grp_rsa.study.targcorr_w.cond{xcond};
                        elseif xx == 3, xpat = grp_rsa.study.targcorr_b.cond{xcond};
                        end
                        
                        %#########################################################
                        %*************** similarity: pearson correlation
                        xtarg_col = []; xtarg_row = [];
                        
                        for xtr = it_trs
                            xtarg_col = vertcat(xtarg_col, xpat.targ{xcol}.tr{xtr});
                            xtarg_row = vertcat(xtarg_row, xpat.targ{xrow}.tr{xtr});
                        end
                        
                        %*************** ttest
                        [~, xpvalue, ~, xstats] = ttest(mean(xtarg_col), mean(xtarg_row), 'Alpha', xalpha);
                        
                        if xx == 1
                            grp_rsa.study.targcorr.ttest.cond{xcond}.xblk{xblk}.pvalue(xrow, xcol)   = xpvalue;
                            grp_rsa.study.targcorr.ttest.cond{xcond}.xblk{xblk}.stats{xrow, xcol}    = xstats;
                        elseif xx == 2
                            grp_rsa.study.targcorr_w.ttest.cond{xcond}.xblk{xblk}.pvalue(xrow, xcol) = xpvalue;
                            grp_rsa.study.targcorr_w.ttest.cond{xcond}.xblk{xblk}.stats{xrow, xcol}  = xstats;
                        elseif xx == 3
                            grp_rsa.study.targcorr_b.ttest.cond{xcond}.xblk{xblk}.pvalue(xrow, xcol) = xpvalue;
                            grp_rsa.study.targcorr_b.ttest.cond{xcond}.xblk{xblk}.stats{xrow, xcol}  = xstats;
                        end
                    end
                end
            end
        end
    end
    
    %% ============= STATS: B/W CONDS WITH ONLY-TARGET
    for xblk = 1:n_stat_blks
        
        grp_rsa.study.onlytargcorr.ttest.blk{xblk}.pvalue    = nan(n_condition, n_condition);
        grp_rsa.study.onlytargcorr_w.ttest.xblk{xblk}.pvalue = nan(n_condition, n_condition);
        grp_rsa.study.onlytargcorr_b.ttest.xblk{xblk}.pvalue = nan(n_condition, n_condition);
        
        for xcol = 1:(n_condition-1)
            for xrow = (xcol+1):n_condition
                it_trs = (1:n_tr_blks) + (n_tr_blks * (xblk-1));
                
                for xx = 1:3
                    clear xpat
                    if xx == 1,     xpat = grp_rsa.study.targcorr;
                    elseif xx == 2, xpat = grp_rsa.study.targcorr_w;
                    elseif xx == 3, xpat = grp_rsa.study.targcorr_b;
                    end
                    
                    %*************** similarity: pearson correlation\
                    xtarg_col = []; xtarg_row = [];
                    
                    for xtr = it_trs
                        xtarg_col = vertcat(xtarg_col, xpat.cond{xcol}.targ{1}.tr{xtr});
                        xtarg_row = vertcat(xtarg_row, xpat.cond{xrow}.targ{1}.tr{xtr});
                    end
                    
                    %*************** ttest
                    [~, xpvalue, ~, xstats] = ttest(mean(xtarg_col), mean(xtarg_row), 'Alpha', xalpha);
                    
                    if xx == 1
                        grp_rsa.study.onlytargcorr.ttest.blk{xblk}.pvalue(xrow, xcol)   = xpvalue;
                        grp_rsa.study.onlytargcorr.ttest.blk{xblk}.stats{xrow, xcol}    = xstats;
                    elseif xx == 2
                        grp_rsa.study.onlytargcorr_w.ttest.blk{xblk}.pvalue(xrow, xcol) = xpvalue;
                        grp_rsa.study.onlytargcorr_w.ttest.blk{xblk}.stats{xrow, xcol}  = xstats;
                    elseif xx == 3
                        grp_rsa.study.onlytargcorr_b.ttest.blk{xblk}.pvalue(xrow, xcol) = xpvalue;
                        grp_rsa.study.onlytargcorr_b.ttest.blk{xblk}.stats{xrow, xcol}  = xstats;
                    end
                end
            end
        end
    end
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%============= STUDY PLOT
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% ============= TIMECOURSE: PER TARGET
    pat_names = {'raw','weighted','beta'};
    rect_w    = 1200;
    fig_rect  = [0 0 rect_w rect_w];
    x_lim     = [1 n_tc_trs];
    x_tick    = 1:n_tc_trs;
    y_lim     = [-0.7 1.3];
    
    for xcond = 1:n_condition
        
        xfig = figure;
        set(xfig, 'Position', fig_rect)
        
        %*************** mean line plots
        fitx = linspace(1, n_tc_trs, n_tc_trs*10);
        
        if xcond == 2, n_targ = 3; else, n_targ = 2; end
        
        for xsubplt = 1:2
            clear xmean xstd fity
            
            %         subplot(1,3, xsubplt)
            subplot(2,1, xsubplt)
            
            for xtarg = 1:n_targ
                if xsubplt == 1
                    xmean{xtarg} = grp_rsa.study.targcorr.cond{xcond}.mean(xtarg, x_tick);
                    xstd{xtarg}  = grp_rsa.study.targcorr.cond{xcond}.se(xtarg, x_tick);
                elseif xsubplt == 2
                    xmean{xtarg} = grp_rsa.study.targcorr_w.cond{xcond}.mean(xtarg, x_tick);
                    xstd{xtarg}  = grp_rsa.study.targcorr_w.cond{xcond}.se(xtarg, x_tick);
                elseif xsubplt == 3
                    xmean{xtarg} = grp_rsa.study.targcorr_b.cond{xcond}.mean(xtarg, x_tick);
                    xstd{xtarg}  = grp_rsa.study.targcorr_b.cond{xcond}.se(xtarg, x_tick);
                end
                
                fity{xtarg}  = interp1(x_tick, xmean{xtarg}, fitx,'spline');
                
                plot(fitx, fity{xtarg}, '-','Color', xtarg_color{xcond}{xtarg}, 'LineWidth', 2); hold on;
            end
            
            %*************** legend
            if xcond~=2
                xlegend = {'target','nontarget'};
            else
                xlegend = {'target','newtarget','nontarget'};
            end
            lg             = legend(xlegend);
            lg.Location    = 'bestoutside';
            legend(xlegend,'AutoUpdate','off')
            grid on
            
            %*************** std error-bar filling
            for xtarg = 1:n_targ
                clear xerr
                xerr(1,:) = xmean{xtarg} - xstd{xtarg};
                xerr(2,:) = xmean{xtarg} + xstd{xtarg};
                
                for i = 1:2, fit_err(:,i) = interp1(x_tick, xerr(i,:), fitx,'spline'); end
                
                in_between = [fit_err(:,1)', fliplr(fit_err(:,2)')];
                h = fill([fitx fliplr(fitx)], in_between, xtarg_color{xcond}{xtarg});
                set(h,'facealpha', .4)
                
                plot(fitx, fit_err(:,1), '-', 'Color', xtarg_color{xcond}{xtarg})
                plot(fitx, fit_err(:,2), '-', 'Color', xtarg_color{xcond}{xtarg})
            end
            
            %*************** stats stars
            n = 0;
            for xcol = 1:(n_targ-1)
                for xrow = (xcol+1):n_targ
                    yy_sig = (y_lim(2)-0.1) - (n * 0.1);
                    text(-3, yy_sig, sprintf('%s vs. %s', xlegend{xcol}, ...
                        xlegend{xrow}), 'FontSize', 12);
                    
                    for xblk = 1:n_stat_blks
                        if xsubplt == 1
                            xpvalue = grp_rsa.study.targcorr.ttest.cond{xcond}.xblk{xblk}.pvalue(xrow, xcol);
                        elseif xsubplt == 2
                            xpvalue = grp_rsa.study.targcorr_w.ttest.cond{xcond}.xblk{xblk}.pvalue(xrow, xcol);
                        elseif xsubplt == 3
                            xpvalue = grp_rsa.study.targcorr_b.ttest.cond{xcond}.xblk{xblk}.pvalue(xrow, xcol);
                        end
                        
                        xx_sig = (n_tr_blks/2) + (n_tr_blks * (xblk-1)) + 0.5 - 0.08;
                        
                        if xpvalue <= xalpha
                            text(xx_sig, yy_sig, '*', 'FontSize', 20);
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
            
            h1 = text(1.5, y_lim(1)+0.05, '(1) stim onset', 'Color', xonset_color, 'FontSize', 12);
            h2 = text(dur_stim + 1.5, y_lim(1)+0.05, '(2) operation onset', 'Color', xonset_color, 'FontSize', 12);
            h3 = text(dur_stim + dur_manip + 1.5, y_lim(1)+0.05, '(3) fixation onset', 'Color', xonset_color, 'FontSize', 12);
            set(h1,'Rotation', 90); set(h2,'Rotation', 90); set(h3,'Rotation', 90);
            
            %*************** set gca
            title(sprintf('RSA timecourse %s: %s (N=%s, p<%1.4f)', ...
                pat_names{xsubplt}, conds_names{xcond}, num2str(length(xsubj_grp)), xalpha),...
                'FontSize',15,'FontWeight','bold');
            xlabel('Volume (tr)', 'FontSize', 12);
            ylabel('Similarity (Z Pearson R)', 'FontSize', 12);
            
            set(gca,'xlim', x_lim, 'ylim', y_lim);
            set(gca,'XTick', x_tick)
            set(gca,'YTick', y_lim(1):0.1:y_lim(end))
            set(gca,'XTickLabel', x_ticklable)
        end
        
        %% *************** save fig
        fig_fname = fullfile(output_dir, sprintf('plot_timecourse_blk%s_%s_cond%d_%s_n%s', ...
            num2str(n_tr_blks), args.phase_name{xph}, xcond, basename, num2str(length(xsubj_grp))));
        
        %     savefig(xfig, sprintf('%s.fig', fig_fname));
        set(gcf,'PaperUnits','inches','PaperPosition', fig_rect/100)
        print('-dpng', '-r300', sprintf('%s.png',fig_fname))
        saveas(xfig, sprintf('%s.png',fig_fname), 'png')
        
        close(xfig);
        
    end
    
    %% ============= TIMECOURSE: TARGET-Only
    fig_rect  = [0 0 rect_w rect_w];
    y_lim     = [-0.8 1.3];
    x_tick    = 1:n_tc_trs;
    x_lim     = [1 n_tc_trs];
    
    sel_conds{1}  = 1:n_condition;
    sel_conds{2}  = [1 2 3];
    sel_conds{3}  = [1 4];
    sel_conds{4}  = [1 2];
    sel_conds{5}  = [2 3];
    sel_conds{6}  = [2 4];
    sel_conds{7}  = [2 5];
    sel_conds{8}  = [4 5];
    
    for xfill = 1%0:1
        for xsel = 1:length(sel_conds)
            
            xsel_conds = sel_conds{xsel};
            
            xfig = figure;
            set(xfig, 'Position', fig_rect)
            
            %*************** mean line plots
            fitx = linspace(1, n_tc_trs, n_tc_trs*10);
            
            for xsubplt = 1:2
                clear xmean xstd fity
                
                %             subplot(1,3, xsubplt)
                subplot(2,1, xsubplt)
                
                for xcond = xsel_conds
                    if xsubplt == 1
                        xmean{xcond} = grp_rsa.study.targcorr.cond{xcond}.mean(1, x_tick);
                        xstd{xcond}  = grp_rsa.study.targcorr.cond{xcond}.se(1, x_tick);
                    elseif xsubplt == 2
                        xmean{xcond} = grp_rsa.study.targcorr_w.cond{xcond}.mean(1, x_tick);
                        xstd{xcond}  = grp_rsa.study.targcorr_w.cond{xcond}.se(1, x_tick);
                    elseif xsubplt == 3
                        xmean{xcond} = grp_rsa.study.targcorr_b.cond{xcond}.mean(1, x_tick);
                        xstd{xcond}  = grp_rsa.study.targcorr_b.cond{xcond}.se(1, x_tick);
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
                lg.Location    = 'bestoutside';
                legend(xlegend,'AutoUpdate','off')
                grid on
                
                if xfill
                    %*************** std error-bar filling
                    for xcond = xsel_conds %#ok<*UNRCH>
                        clear xerr fit_err
                        xerr(1,:) = xmean{xcond} - xstd{xcond};
                        xerr(2,:) = xmean{xcond} + xstd{xcond};
                        
                        for i = 1:2, fit_err(:,i) = interp1(x_tick, xerr(i,:), fitx,'spline'); end
                        
                        in_between = [fit_err(:,1)', fliplr(fit_err(:,2)')];
                        h = fill([fitx fliplr(fitx)], in_between, xcond_color{xcond});
                        h.FaceAlpha = .4;
                        h.EdgeAlpha = .4;
                        h.EdgeColor = xcond_color{xcond};
                        
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
                        
                        yy_sig = (y_lim(2)-0.1) - (n * 0.1);
                        
                        if xsel==1
                            text(-3, yy_sig, sprintf('%d vs. %d', xcol, xrow), 'FontSize', 12);
                        end
                        
                        for xblk = 1:n_stat_blks
                            
                            if xsubplt == 1
                                xpvalue = grp_rsa.study.onlytargcorr.ttest.blk{xblk}.pvalue(xrow, xcol);
                            elseif xsubplt == 2
                                xpvalue = grp_rsa.study.onlytargcorr_w.ttest.blk{xblk}.pvalue(xrow, xcol);
                            elseif xsubplt == 3
                                xpvalue = grp_rsa.study.onlytargcorr_b.ttest.blk{xblk}.pvalue(xrow, xcol);
                            end
                            
                            xx_sig = (n_tr_blks/2) + (n_tr_blks * (xblk-1)) + 0.5 - 0.08;
                            if xpvalue <= xalpha
                                text(xx_sig, yy_sig, '*', 'FontSize', 20);
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
                
                h1 = text(1.5, y_lim(1)+0.05, '(1) stim onset', 'Color', xonset_color, 'FontSize', 12);
                h2 = text(dur_stim + 1.5, y_lim(1)+0.05, '(2) operation onset', 'Color', xonset_color, 'FontSize', 12);
                h3 = text(dur_stim + dur_manip + 1.5, y_lim(1)+0.05, '(3) fixation onset', 'Color', xonset_color, 'FontSize', 12);
                set(h1,'Rotation', 90); set(h2,'Rotation', 90); set(h3,'Rotation', 90);
                
                %*************** set gca
                title(sprintf('RSA target-timecourse %s (N=%s, p<%1.4f)', ...
                    pat_names{xsubplt}, num2str(length(xsubj_grp)), xalpha),...
                    'FontSize',15,'FontWeight','bold');
                xlabel('Volume (tr)', 'FontSize', 12);
                ylabel('Similarity (Pearson R)', 'FontSize', 12);
                
                set(gca,'xlim', x_lim, 'ylim', y_lim);
                set(gca,'XTick', x_tick)
                set(gca,'YTick', y_lim(1):0.1:y_lim(end))
                set(gca,'XTickLabel', x_ticklable)
            end
            
            %% *************** save fig
            fig_fname = fullfile(output_dir, sprintf('plot_target_timecourse_blk%s_%s_sel%d_fill%d_%s_n%s', ...
                num2str(n_tr_blks), args.phase_name{xph}, xsel, xfill, basename, num2str(length(xsubj_grp))));
            
            savefig(xfig, sprintf('%s.fig', fig_fname));
            set(gcf,'PaperUnits','inches','PaperPosition', fig_rect/100)
            print('-dpng', '-r300', sprintf('%s.png',fig_fname))
            saveas(xfig, sprintf('%s.png',fig_fname), 'png')
            
            close(xfig);
            
        end
    end
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%============= TIMECOURSE: SYNC PER TARGET
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %*************** target/ non-target
    % grp_rsa.study.timecourse.cond{xcond}.corr.mean.raw(xtarg, xtr)
    % cond~=2: targ:1, non-targ:2
    % cond==2: targ:1, new_targ:2, non-targ:3
    for xcond = 1:n_condition
        
        if xcond == 2, n_targ = 3; else, n_targ = 2; end
        
        %*************** sync trs
        grp_rsa.study.targcorr.sync.cond{xcond}.mean   = zeros(n_targ, dur_sync);
        grp_rsa.study.targcorr.sync.cond{xcond}.se     = zeros(n_targ, dur_sync);
        grp_rsa.study.targcorr_w.sync.cond{xcond}.mean = zeros(n_targ, dur_sync);
        grp_rsa.study.targcorr_w.sync.cond{xcond}.se   = zeros(n_targ, dur_sync);
        grp_rsa.study.targcorr_b.sync.cond{xcond}.mean = zeros(n_targ, dur_sync);
        grp_rsa.study.targcorr_b.sync.cond{xcond}.se   = zeros(n_targ, dur_sync);
        
        for xtr = 1:dur_sync
            if xcond~=2
                %*************** 1 target: corr b/w target template & target trial (study)
                % e.g.: face_temp & face trial under face mask
                xtarg = 1;
                xtarg_corr = []; xtarg_corr_w = []; xtarg_corr_b = [];
                
                for xcate = 1:n_category
                    xtarg_corr = vertcat(xtarg_corr, ...
                        grp_rsa.study.timecourse.sync.cond{xcond}.targcate{xcate}.targ{xtarg}.corr.raw{xtr});
                    xtarg_corr_w = vertcat(xtarg_corr_w, ...
                        grp_rsa.study.timecourse.sync.cond{xcond}.targcate{xcate}.targ{xtarg}.corr.weighted{xtr});
                    xtarg_corr_b = vertcat(xtarg_corr_b, ...
                        grp_rsa.study.timecourse.sync.cond{xcond}.targcate{xcate}.targ{xtarg}.corr.beta{xtr});
                end
                
                grp_rsa.study.targcorr.sync.cond{xcond}.targ{xtarg}.tr{xtr}   = mean(xtarg_corr);
                grp_rsa.study.targcorr_w.sync.cond{xcond}.targ{xtarg}.tr{xtr} = mean(xtarg_corr_w);
                grp_rsa.study.targcorr_b.sync.cond{xcond}.targ{xtarg}.tr{xtr} = mean(xtarg_corr_b);
                
                grp_rsa.study.targcorr.sync.cond{xcond}.mean(xtarg, xtr)      = mean(mean(xtarg_corr));
                grp_rsa.study.targcorr.sync.cond{xcond}.se(xtarg, xtr)        = std(mean(xtarg_corr))/sqrt(length(xsubj_grp));
                
                grp_rsa.study.targcorr_w.sync.cond{xcond}.mean(xtarg, xtr)    = mean(mean(xtarg_corr_w));
                grp_rsa.study.targcorr_w.sync.cond{xcond}.se(xtarg, xtr)      = std(mean(xtarg_corr_w))/sqrt(length(xsubj_grp));
                
                grp_rsa.study.targcorr_b.sync.cond{xcond}.mean(xtarg, xtr)    = mean(mean(xtarg_corr_b));
                grp_rsa.study.targcorr_b.sync.cond{xcond}.se(xtarg, xtr)      = std(mean(xtarg_corr_b))/sqrt(length(xsubj_grp));
                
                %*************** 2 nontarget
                xtarg_corr = []; xtarg_corr_w = []; xtarg_corr_b = [];
                
                for xtarg = 2:3
                    for xcate = 1:n_category
                        xtarg_corr = vertcat(xtarg_corr, ...
                            grp_rsa.study.timecourse.sync.cond{xcond}.targcate{xcate}.targ{xtarg}.corr.raw{xtr});
                        xtarg_corr_w = vertcat(xtarg_corr_w, ...
                            grp_rsa.study.timecourse.sync.cond{xcond}.targcate{xcate}.targ{xtarg}.corr.weighted{xtr});
                        xtarg_corr_b = vertcat(xtarg_corr_b, ...
                            grp_rsa.study.timecourse.sync.cond{xcond}.targcate{xcate}.targ{xtarg}.corr.beta{xtr});
                    end
                end
                
                xtarg = 2;
                grp_rsa.study.targcorr.sync.cond{xcond}.targ{xtarg}.tr{xtr}   = mean(xtarg_corr);
                grp_rsa.study.targcorr_w.sync.cond{xcond}.targ{xtarg}.tr{xtr} = mean(xtarg_corr_w);
                grp_rsa.study.targcorr_b.sync.cond{xcond}.targ{xtarg}.tr{xtr} = mean(xtarg_corr_b);
                
                grp_rsa.study.targcorr.sync.cond{xcond}.mean(xtarg, xtr)      = mean(mean(xtarg_corr));
                grp_rsa.study.targcorr.sync.cond{xcond}.se(xtarg, xtr)        = std(mean(xtarg_corr))/sqrt(length(xsubj_grp));
                
                grp_rsa.study.targcorr_w.sync.cond{xcond}.mean(xtarg, xtr)    = mean(mean(xtarg_corr_w));
                grp_rsa.study.targcorr_w.sync.cond{xcond}.se(xtarg, xtr)      = std(mean(xtarg_corr_w))/sqrt(length(xsubj_grp));
                
                grp_rsa.study.targcorr_b.sync.cond{xcond}.mean(xtarg, xtr)    = mean(mean(xtarg_corr_b));
                grp_rsa.study.targcorr_b.sync.cond{xcond}.se(xtarg, xtr)      = std(mean(xtarg_corr_b))/sqrt(length(xsubj_grp));
                
            else
                for xtarg = 1:3
                    %*************** 2 newtarg, 3 nontarget
                    xtarg_corr = []; xtarg_corr_w = []; xtarg_corr_b = [];
                    
                    for xcate = 1:n_category
                        for xnewcate = cate_members(~ismember(cate_members, xcate))
                            xtarg_corr = vertcat(xtarg_corr, ...
                                grp_rsa.study.timecourse.sync.cond{xcond}.targcate{xcate}.newtarg{xnewcate}.targ{xtarg}.corr.raw{xtr});
                            xtarg_corr_w = vertcat(xtarg_corr_w, ...
                                grp_rsa.study.timecourse.sync.cond{xcond}.targcate{xcate}.newtarg{xnewcate}.targ{xtarg}.corr.weighted{xtr});
                            xtarg_corr_b = vertcat(xtarg_corr_b, ...
                                grp_rsa.study.timecourse.sync.cond{xcond}.targcate{xcate}.newtarg{xnewcate}.targ{xtarg}.corr.beta{xtr});
                        end
                    end
                    
                    grp_rsa.study.targcorr.sync.cond{xcond}.targ{xtarg}.tr{xtr}   = mean(xtarg_corr);
                    grp_rsa.study.targcorr_w.sync.cond{xcond}.targ{xtarg}.tr{xtr} = mean(xtarg_corr_w);
                    grp_rsa.study.targcorr_b.sync.cond{xcond}.targ{xtarg}.tr{xtr} = mean(xtarg_corr_b);
                    
                    grp_rsa.study.targcorr.sync.cond{xcond}.mean(xtarg, xtr)      = mean(mean(xtarg_corr));
                    grp_rsa.study.targcorr.sync.cond{xcond}.se(xtarg, xtr)        = std(mean(xtarg_corr))/sqrt(length(xsubj_grp));
                    
                    grp_rsa.study.targcorr_w.sync.cond{xcond}.mean(xtarg, xtr)    = mean(mean(xtarg_corr_w));
                    grp_rsa.study.targcorr_w.sync.cond{xcond}.se(xtarg, xtr)      = std(mean(xtarg_corr_w))/sqrt(length(xsubj_grp));
                    
                    grp_rsa.study.targcorr_b.sync.cond{xcond}.mean(xtarg, xtr)    = mean(mean(xtarg_corr_b));
                    grp_rsa.study.targcorr_b.sync.cond{xcond}.se(xtarg, xtr)      = std(mean(xtarg_corr_b))/sqrt(length(xsubj_grp));
                    
                end
            end
        end
    end
    
    %% ============= STATS: B/W CONDS WITH ONLY-TARGET
    for xblk = 1:(dur_sync/n_tr_blks)
        
        grp_rsa.study.onlytargcorr.ttest.sync.blk{xblk}.pvalue    = nan(n_condition, n_condition);
        grp_rsa.study.onlytargcorr_w.ttest.sync.xblk{xblk}.pvalue = nan(n_condition, n_condition);
        grp_rsa.study.onlytargcorr_b.ttest.sync.xblk{xblk}.pvalue = nan(n_condition, n_condition);
        
        for xcol = 1:(n_condition-1)
            for xrow = (xcol+1):n_condition
                it_trs = (1:n_tr_blks) + (n_tr_blks * (xblk-1));
                
                for xx = 1:3
                    clear xpat
                    if xx == 1,     xpat = grp_rsa.study.targcorr.sync;
                    elseif xx == 2, xpat = grp_rsa.study.targcorr_w.sync;
                    elseif xx == 3, xpat = grp_rsa.study.targcorr_b.sync;
                    end
                    
                    %*************** similarity: pearson correlation\
                    xtarg_col = []; xtarg_row = [];
                    
                    for xtr = it_trs
                        xtarg_col = vertcat(xtarg_col, xpat.cond{xcol}.targ{1}.tr{xtr});
                        xtarg_row = vertcat(xtarg_row, xpat.cond{xrow}.targ{1}.tr{xtr});
                    end
                    
                    %*************** ttest
                    [~, xpvalue, ~, xstats] = ttest(mean(xtarg_col), mean(xtarg_row), 'Alpha', xalpha);
                    
                    if xx == 1
                        grp_rsa.study.onlytargcorr.ttest.sync.blk{xblk}.pvalue(xrow, xcol)   = xpvalue;
                        grp_rsa.study.onlytargcorr.ttest.sync.blk{xblk}.stats{xrow, xcol}    = xstats;
                    elseif xx == 2
                        grp_rsa.study.onlytargcorr_w.ttest.sync.blk{xblk}.pvalue(xrow, xcol) = xpvalue;
                        grp_rsa.study.onlytargcorr_w.ttest.sync.blk{xblk}.stats{xrow, xcol}  = xstats;
                    elseif xx == 3
                        grp_rsa.study.onlytargcorr_b.ttest.sync.blk{xblk}.pvalue(xrow, xcol) = xpvalue;
                        grp_rsa.study.onlytargcorr_b.ttest.sync.blk{xblk}.stats{xrow, xcol}  = xstats;
                    end
                end
            end
        end
    end
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%============= SYNC STUDY PLOT
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% ============= TIMECOURSE: TARGET-Only
    fig_rect  = [0 0 (rect_w/2.3) rect_w];
    y_lim     = [-0.8 1.3];
    x_tick    = 1:dur_sync;
    x_lim     = [1 dur_sync];
    
    sel_conds{1} = 1:n_condition;
    sel_conds{2} = [1 4];
    sel_conds{3} = [2 3];
    sel_conds{4} = [1 2];
    sel_conds{5} = [4 5];
    
    for xfill = 0:1
        for xsel = 1:length(sel_conds)
            
            xsel_conds = sel_conds{xsel};
            
            xfig = figure;
            set(xfig, 'Position', fig_rect)
            
            %*************** mean line plots
            fitx = linspace(1, dur_sync, dur_sync*10);
            
            for xsubplt = 1:2
                clear xmean xstd fity
                
                subplot(2,1, xsubplt)
                
                for xcond = xsel_conds
                    if xsubplt == 1
                        xmean{xcond} = grp_rsa.study.targcorr.sync.cond{xcond}.mean(1, x_tick);
                        xstd{xcond}  = grp_rsa.study.targcorr.sync.cond{xcond}.se(1, x_tick);
                    elseif xsubplt == 2
                        xmean{xcond} = grp_rsa.study.targcorr_w.sync.cond{xcond}.mean(1, x_tick);
                        xstd{xcond}  = grp_rsa.study.targcorr_w.sync.cond{xcond}.se(1, x_tick);
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
                lg.Location    = 'bestoutside';
                legend(xlegend,'AutoUpdate','off')
                grid on
                
                if xfill
                    %*************** std error-bar filling
                    for xcond = xsel_conds %#ok<*UNRCH>
                        clear xerr fit_err
                        xerr(1,:) = xmean{xcond} - xstd{xcond};
                        xerr(2,:) = xmean{xcond} + xstd{xcond};
                        
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
                        
                        yy_sig = (y_lim(2)-0.1) - (n * 0.1);
                        
                        for xblk = 1:(dur_sync/n_tr_blks)
                            
                            if xsubplt == 1
                                xpvalue = grp_rsa.study.onlytargcorr.ttest.sync.blk{xblk}.pvalue(xrow, xcol);
                            elseif xsubplt == 2
                                xpvalue = grp_rsa.study.onlytargcorr_w.ttest.sync.blk{xblk}.pvalue(xrow, xcol);
                            end
                            
                            xx_sig = (n_tr_blks/2) + (n_tr_blks * (xblk-1)) + 0.5 - 0.08;
                            if xpvalue <= xalpha
                                text(xx_sig, yy_sig, '*', 'FontSize', 20);
                            end
                            
                            xx = n_tr_blks + (n_tr_blks * (xblk-1)) + 0.5;
                            plot([xx xx], y_lim, '--', 'Color', 'b')
                            
                        end
                        
                        n = n + 1;
                    end
                end
                
                %*************** real onset lines
                plot([dur_stim dur_stim]+1, y_lim, '-', 'Color', xonset_color)
                
                h1 = text(1.5, y_lim(1)+0.05, '(1) stim onset', 'Color', xonset_color, 'FontSize', 12);
                h2 = text(dur_stim + 1.5, y_lim(1)+0.05, '(2) operation onset', 'Color', xonset_color, 'FontSize', 12);
                set(h1,'Rotation', 90); set(h2,'Rotation', 90);
                
                %*************** set gca
                title(sprintf('RSA target-timecourse %s (N=%s, p<%1.4f)', ...
                    pat_names{xsubplt}, num2str(length(xsubj_grp)), xalpha),...
                    'FontSize',15,'FontWeight','bold');
                xlabel('Volume (tr)', 'FontSize', 12);
                ylabel('Similarity (Pearson R)', 'FontSize', 12);
                
                set(gca,'xlim', x_lim, 'ylim', y_lim);
                set(gca,'XTick', x_tick)
                set(gca,'YTick', y_lim(1):0.1:y_lim(end))
                set(gca,'XTickLabel', x_ticklable)
            end
            
            %% *************** save fig
            fig_fname = fullfile(output_dir, sprintf('plot_target_sync_timecourse_blk%s_%s_sel%d_fill%d_%s_n%s', ...
                num2str(n_tr_blks), args.phase_name{xph}, xsel, xfill, basename, num2str(length(xsubj_grp))));
            
            savefig(xfig, sprintf('%s.fig', fig_fname));
            set(gcf,'PaperUnits','inches','PaperPosition', fig_rect/100)
            print('-dpng', '-r300', sprintf('%s.png',fig_fname))
            saveas(xfig, sprintf('%s.png',fig_fname), 'png')
            
            close(xfig);
            
        end
    end
    
    %% ============= SAVE GROUP RSA
    fprintf('\n....saving grp_rsa patterns\n')
    save(g_fname, 'grp_rsa','-v7.3')
    
    g_fsize = dir(g_fname);
    fprintf('file size: %s GB\n', num2str(g_fsize.bytes/(10^9)));
end

