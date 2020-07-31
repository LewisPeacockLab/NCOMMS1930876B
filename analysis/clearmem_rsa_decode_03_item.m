function[] = clearmem_rsa_decode_03_item(args, dirs)
%*************** RSA timecourse: test on study

%% ============= SETUP DIRECTORY
output_dir      = dirs.rsa.group.parse;
xsubj_grp       = args.filtered_subs;

%% ============= SETUP PARAMETERS
n_runs{1}       = 5;
n_runs{2}       = 6;
        
xparam          = args.index{1}.param;
n_category      = xparam.n_category;
n_subcategory   = xparam.n_subcategory;
n_targ          = 3;%1_targ, 2_nontarg, 3_baseline
n_item          = xparam.n_item;

category_names  = {'face','fruit','scene'};
subcate_names   = xparam.subcategory_name;
conds_names     = {'maintain','replace (category)','replace (subcategory)','suppress','clear'};
n_condition     = length(conds_names);
n_trs           = args.tc_tr_disp;
cate_members    = 1:n_category;

xmatrix = args.index{1}.matrix;
xheader = args.index{1}.header;

for xcate = 1:n_category
    for xsubcate = 1:n_subcategory
        for xitem = 1:n_item
            it_runs{xcate}{xsubcate}{xitem} = ...
                unique(getDATA(xmatrix', xheader, ...
                {'category','subcategory','item','stimulus'}, ...
                {xcate, xsubcate, xitem, 1}, findCol(xheader,{'run'}))); %#ok<*AGROW>
            
            nn_runs{xcate}{xsubcate}(xitem) = length(it_runs{xcate}{xsubcate}{xitem});
            
            %*************** subcategory name
            item_names{xcate}{xsubcate}{xitem} = ...
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
xparam    = args.index{2}.param;
dur_sync  = args.dur_sync;
n_tr_blks = args.n_tr_blks;%trs in one stats block: 5=2.3s/3=1.38s

%% ============= LOGGING
fprintf('(+) using scratch dir: %s\n', output_dir);
%*************** turn on diary to capture analysis output
diary off;
diary(sprintf('%s/diary.txt', output_dir));
fprintf('running code: %s at %s\n\n', mfilename, datestr(now, 0))
fprintf('#####################################################################\n\n');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%============= 1ST LEVEL RSA > GROUP STRUCTURE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
g_fname = fullfile(dirs.rsa.group.pattern, sprintf('group_%s.mat', basename));

if args.grp_pattern
    %% ============= SETUP STRUCTURE  
    %*************** LOCALIZER
    xph = 1;
    
    for xcate = 1:n_category
        for xsubcate = 1:n_subcategory
            
            n_units = sum(nn_runs{xcate}{xsubcate});
            
            grp_rsa{xph}.item.cate{xcate}.subcate{xsubcate}.corr   = cell(n_units, n_units);
            grp_rsa{xph}.item.cate{xcate}.subcate{xsubcate}.corr_w = cell(n_units, n_units);
        end
    end
    
    for xcate = 1:n_category
        n_units = n_category * n_subcategory * n_item;
        grp_rsa{xph}.template.maskcate{xcate}.corr   = cell(n_units, n_units);
        grp_rsa{xph}.template.maskcate{xcate}.corr_w = cell(n_units, n_units);
    end
    
    for xcate = 1:n_category
        for xsubcate = 1:n_subcategory
            for xtarg = 1:3
                grp_rsa{xph}.item.cate{xcate}.subcate{xsubcate}.targ{xtarg} = [];
            end
        end
    end
    
    %*************** STUDY
    xph = 2;
    for xcond = 1:n_condition
        for xtarg = 1:n_targ
            for xtr = 1:n_trs
                %*************** targ vs. nontarg
                grp_rsa{xph}.timecourse.cond{xcond}.targ{xtarg}.tr{xtr}.corr_w = [];
                
                for xcate = 1:n_category
                    for xsubcate = 1:n_subcategory
                        grp_rsa{xph}.timecourse.cond{xcond}.cate{xcate}.subcate{xsubcate}.targ{xtarg}.tr{xtr}.corr_w = [];
                    end
                end
                
                %*************** same vs. diff
                for xlevel = 1:2%1_cate, 2_subcate
                    
                    if xlevel==1, it_sames = 1:2; else, it_sames = 1:3; end
                    
                    for xsame = it_sames % 1_same, 2_differ, 3_related
                        grp_rsa{xph}.timecourse.sameseq{xlevel}{xsame}.cond{xcond}.targ{xtarg}.tr{xtr}.corr_w = [];
                        
                        %*************** n trials
                        if (xtarg==1) && (xtr==1)
                            grp_rsa{xph}.timecourse.sameseq{xlevel}{xsame}.cond{xcond}.n_trial = [];
                        end
                    end
                end
            end
            
            for xtr = 1:dur_sync
                grp_rsa{xph}.timecourse.sync.cond{xcond}.targ{xtarg}.tr{xtr}.corr_w = [];
                
                for xlevel = 1:2%1_cate, 2_subcate
                    
                    if xlevel==1, it_sames = 1:2; else, it_sames = 1:3; end
                    
                    for xsame = it_sames % 1_same, 2_differ, 3_related
                        grp_rsa{xph}.timecourse.sameseq{xlevel}{xsame}.sync.cond{xcond}.targ{xtarg}.tr{xtr}.corr_w = [];
                        grp_rsa{xph}.timecourse.sameseq{xlevel}{xsame}.sync.next.cond{xcond}.targ{xtarg}.tr{xtr}.corr_w = [];
                    end
                end
            end
        end
    end
    
    %% ============= LOAD SUBJECT DATA
    for it = 1:length(xsubj_grp)
        xsub = xsubj_grp(it);
        clear sub_args subj_rsa sub_dirs rsa_fname
        
        %*************** setup subject & directories
        sub_args             = args;
        sub_args.subject_num = xsub;
        sub_args.subject_id  = args.subject_list(xsub).name;
        sub_dirs             = setup_directory(dirs, sub_args);
        
        %*************** load subj_rsa
        
        fprintf('(+) reseting grp_rsa patterns: %s\n', args.subject_list(xsub).name)
        
        if args.item_within
            rsa_fname = fullfile(sub_dirs.rsa.pattern, sprintf('%s_within.mat', basename));
        else
            rsa_fname = fullfile(sub_dirs.rsa.pattern, sprintf('%s.mat', basename));
        end
        load(rsa_fname);% 'subj_rsa'
        
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%*************** BETA WEIGHT
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        for xcate = 1:n_category
            grp_rsa{1}.subj{xsub}.item_beta{xcate} = subj_rsa{1}.subj{xcate};
        end
        
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%*************** LOCALIZER
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        xph = 1;
        
        for xcate = 1:n_category
            for xsubcate = 1:n_subcategory
                for xitem = 1:n_item
                    for xrun = 1:nn_runs{xcate}{xsubcate}(xitem)
                        
                        it_xrun = it_runs{xcate}{xsubcate}{xitem}(xrun);
                        
                        xcell    = xrun + sum(nn_runs{xcate}{xsubcate}(1:(xitem-1)));
                        it_xcell = it_xrun + (n_runs{xph} * (xitem-1));
                        
                        for yitem = 1:n_item
                            for yrun = 1:nn_runs{xcate}{xsubcate}(yitem)
                            
                                it_yrun = it_runs{xcate}{xsubcate}{yitem}(yrun);
                                
                                ycell    = yrun + sum(nn_runs{xcate}{xsubcate}(1:(yitem-1)));
                                it_ycell = it_yrun + (n_runs{xph} * (yitem-1));
                                
                                grp_rsa{xph}.item.cate{xcate}.subcate{xsubcate}.corr{xcell, ycell} = ...
                                    horzcat(grp_rsa{xph}.item.cate{xcate}.subcate{xsubcate}.corr{xcell, ycell}, ...
                                    subj_rsa{xph}.rsa.item.cate{xcate}.subcate{xsubcate}.corr{it_xcell, it_ycell}); %#ok<*NODEF,*USENS>
                                
                                grp_rsa{xph}.item.cate{xcate}.subcate{xsubcate}.corr_w{xcell, ycell} = ...
                                    horzcat(grp_rsa{xph}.item.cate{xcate}.subcate{xsubcate}.corr_w{xcell, ycell}, ...
                                    subj_rsa{xph}.rsa.item.cate{xcate}.subcate{xsubcate}.corr_w{it_xcell, it_ycell});
                            end
                        end
                    end
                end
            end
        end
        
        %% ************** within/between items
        % grp_rsa{xph}.item.cate{xcate}.subcate{xsubcate}.ttest.mean(xsub, xcol)
        % xcol: 1_within, 2_between
        clear xwithin xbetween
        for xcate = 1:n_category
            for xsubcate = 1:n_subcategory
                xwithin{xcate}{xsubcate} = []; xbetween{xcate}{xsubcate} = [];
            end
        end
        
        for xcate = 1:n_category
            for xsubcate = 1:n_subcategory
                for xitem = 1:n_item
                    for xrun = 1:nn_runs{xcate}{xsubcate}(xitem)
                        
                        it_xrun = it_runs{xcate}{xsubcate}{xitem}(xrun);
                        
                        xcell    = xrun + sum(nn_runs{xcate}{xsubcate}(1:(xitem-1)));
                        it_xcell = it_xrun + (n_runs{xph} * (xitem-1));
                        
                        for yitem = xitem:n_item
                            
                            if xitem == yitem, yy_runs = (xrun+1):nn_runs{xcate}{xsubcate}(yitem);
                            else,              yy_runs = 1:nn_runs{xcate}{xsubcate}(yitem); end
                            
                            for yrun = yy_runs
                                
                                it_yrun = it_runs{xcate}{xsubcate}{yitem}(yrun);
                                
                                ycell    = yrun + sum(nn_runs{xcate}{xsubcate}(1:(yitem-1)));
                                it_ycell = it_yrun + (n_runs{xph} * (yitem-1));
                                
                                %************** within/between items
                                if xitem == yitem % within items
                                    xwithin{xcate}{xsubcate} = horzcat(xwithin{xcate}{xsubcate}, ...
                                        subj_rsa{xph}.rsa.item.cate{xcate}.subcate{xsubcate}.corr_w{it_xcell, it_ycell});
                                else % between items
                                    xbetween{xcate}{xsubcate} = horzcat(xbetween{xcate}{xsubcate}, ...
                                        subj_rsa{xph}.rsa.item.cate{xcate}.subcate{xsubcate}.corr_w{it_xcell, it_ycell});
                                end 
                            end
                        end
                    end
                end
            end
        end
        
        for xcate = 1:n_category
            for xsubcate = 1:n_subcategory
                grp_rsa{xph}.item.ttest.cate{xcate}.subcate{xsubcate}.mean(it, 1) = mean(xwithin{xcate}{xsubcate});
                grp_rsa{xph}.item.ttest.cate{xcate}.subcate{xsubcate}.mean(it, 2) = mean(xbetween{xcate}{xsubcate});
            end
        end
        
        %% *************** template
        for xmaskcate = 1:n_category
            for xcate = 1:n_category
                for xsubcate = 1:n_subcategory
                    for xitem = 1:n_item
                        
                        xcell = xitem + (n_item * (xsubcate-1)) + (n_item * n_subcategory * (xcate-1));
                        
                        for ycate = 1:n_category
                            for ysubcate = 1:n_subcategory
                                for yitem = 1:n_item
                                    
                                    ycell = yitem + (n_item * (ysubcate-1)) + (n_item * n_subcategory * (ycate-1));
                                    
                                    grp_rsa{xph}.template.maskcate{xmaskcate}.corr{xcell, ycell} = ...
                                        horzcat(grp_rsa{xph}.template.maskcate{xmaskcate}.corr{xcell, ycell},...
                                        subj_rsa{xph}.rsa.template.maskcate{xmaskcate}.corr{xcell, ycell});
                                    
                                    grp_rsa{xph}.template.maskcate{xmaskcate}.corr_w{xcell, ycell} = ...
                                        horzcat(grp_rsa{xph}.template.maskcate{xmaskcate}.corr_w{xcell, ycell},...
                                        subj_rsa{xph}.rsa.template.maskcate{xmaskcate}.corr_w{xcell, ycell});
                                end
                                
                            end
                        end
                    end
                end
            end
        end
        
        %% *************** ITEM ACROSS CATEGORY: targ vs. non-targ
        for xcate = 1:n_category
            for xsubcate = 1:n_subcategory
                for xtarg = 1:3
                    grp_rsa{xph}.item.cate{xcate}.subcate{xsubcate}.targ{xtarg} = ...
                        horzcat(grp_rsa{xph}.item.cate{xcate}.subcate{xsubcate}.targ{xtarg},...
                        mean(subj_rsa{xph}.rsa.item.cate{xcate}.subcate{xsubcate}.targ{xtarg}));
                end
            end
        end
                
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%*************** STUDY
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        xph = 2;
        
        grp_rsa{xph}.subj{xsub}.timecourse = subj_rsa{xph}.timecourse;
        grp_rsa{xph}.subj{xsub}.react      = subj_rsa{xph}.react;

        for xcond = 1:n_condition
            for xtarg = 1:n_targ
                for xtr = 1:n_trs
                    %*************** targ vs. nontarg
                    grp_rsa{xph}.timecourse.cond{xcond}.targ{xtarg}.tr{xtr}.corr_w = ...
                        horzcat(grp_rsa{xph}.timecourse.cond{xcond}.targ{xtarg}.tr{xtr}.corr_w, ...
                        subj_rsa{xph}.rsa.timecourse.cond{xcond}.targ{xtarg}.tr{xtr}.corr_w);
                    %*************** targ vs. nontarg: per cate/subcate
                    for xcate = 1:n_category
                        for xsubcate = 1:n_subcategory
                            grp_rsa{xph}.timecourse.cond{xcond}.cate{xcate}.subcate{xsubcate}.targ{xtarg}.tr{xtr}.corr_w = ...
                                horzcat(grp_rsa{xph}.timecourse.cond{xcond}.cate{xcate}.subcate{xsubcate}.targ{xtarg}.tr{xtr}.corr_w, ...
                                mean(subj_rsa{xph}.timecourse.cond{xcond}.cate{xcate}.subcate{xsubcate}.targ{xtarg}.tr{xtr}.corr_w));
                        end
                    end
                    
                    for xlevel = 1:2%1_cate, 2_subcate
                        
                        if xlevel==1, it_sames = 1:2; else, it_sames = 1:3; end
                
                        for xsame = it_sames % 1_same, 2_differ, 3_related
                            grp_rsa{xph}.timecourse.sameseq{xlevel}{xsame}.cond{xcond}.targ{xtarg}.tr{xtr}.corr_w = ...
                                horzcat(grp_rsa{xph}.timecourse.sameseq{xlevel}{xsame}.cond{xcond}.targ{xtarg}.tr{xtr}.corr_w, ...
                                subj_rsa{xph}.rsa.timecourse.sameseq{xlevel}{xsame}.cond{xcond}.targ{xtarg}.tr{xtr}.corr_w);
                            
                            %*************** n trials
                            if (xtarg==1) && (xtr==1)
                                grp_rsa{xph}.timecourse.sameseq{xlevel}{xsame}.cond{xcond}.n_trial = ...
                                    horzcat(grp_rsa{xph}.timecourse.sameseq{xlevel}{xsame}.cond{xcond}.n_trial, ...
                                    subj_rsa{xph}.rsa.timecourse.sameseq{xlevel}{xsame}.cond{xcond}.n_trial);
                            end
                        end
                    end
                end
                
                %*************** sync TRs
                for xtr = 1:dur_sync
                    
                    grp_rsa{xph}.timecourse.sync.cond{xcond}.targ{xtarg}.tr{xtr}.corr_w = ...
                        horzcat(grp_rsa{xph}.timecourse.sync.cond{xcond}.targ{xtarg}.tr{xtr}.corr_w, ...
                        subj_rsa{xph}.rsa.timecourse.sync.cond{xcond}.targ{xtarg}.tr{xtr}.corr_w);
                    
                    for xlevel = 1:2%1_cate, 2_subcate
                        
                        if xlevel==1, it_sames = 1:2; else, it_sames = 1:3; end
                
                        for xsame = it_sames % 1_same, 2_differ, 3_related
                            grp_rsa{xph}.timecourse.sameseq{xlevel}{xsame}.sync.cond{xcond}.targ{xtarg}.tr{xtr}.corr_w = ...
                                horzcat(grp_rsa{xph}.timecourse.sameseq{xlevel}{xsame}.sync.cond{xcond}.targ{xtarg}.tr{xtr}.corr_w, ...
                                subj_rsa{xph}.rsa.timecourse.sameseq{xlevel}{xsame}.sync.cond{xcond}.targ{xtarg}.tr{xtr}.corr_w);
                            
                            %*************** N+1
                            grp_rsa{xph}.timecourse.sameseq{xlevel}{xsame}.sync.next.cond{xcond}.targ{xtarg}.tr{xtr}.corr_w = ...
                                horzcat(grp_rsa{xph}.timecourse.sameseq{xlevel}{xsame}.sync.next.cond{xcond}.targ{xtarg}.tr{xtr}.corr_w, ...
                                subj_rsa{xph}.rsa.timecourse.sameseq{xlevel}{xsame}.sync.next.cond{xcond}.targ{xtarg}.tr{xtr}.corr_w);
                            
                        end
                    end
                end
            end
        end 
    end
    
    %% ============= SAVE
    
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

%% ============= ADDITIONAL
g_fname_add = fullfile(dirs.rsa.group.pattern, sprintf('group_add_%s.mat', basename));

if args.additional_rsa
    
    for it = 1:length(xsubj_grp)
        xsub = xsubj_grp(it);
        clear sub_args subj_rsa sub_dirs rsa_fname
        
        %*************** setup subject & directories
        sub_args             = args;
        sub_args.subject_num = xsub;
        sub_args.subject_id  = args.subject_list(xsub).name;
        sub_dirs             = setup_directory(dirs, sub_args);
        
        %*************** load subj_rsa
        
        fprintf('(+) reseting grp_rsa_add patterns: %s\n', args.subject_list(xsub).name)
        
        if args.item_within
            rsa_fname = fullfile(sub_dirs.rsa.pattern, sprintf('%s_within.mat', basename));
        else
            rsa_fname = fullfile(sub_dirs.rsa.pattern, sprintf('%s.mat', basename));
        end
        load(rsa_fname);% 'subj_rsa'
        
        xph = 2;
        grp_rsa_add{xph}.subj{xsub}.decoding = subj_rsa{xph}.decoding;
    end
    
    %% ============= SAVE
    
    fprintf('\n....saving grp_rsa_add patterns\n')
    save(g_fname_add, 'grp_rsa_add','-v7.3')
    
    g_fsize = dir(g_fname_add);
    fprintf('file size: %s GB\n', num2str(g_fsize.bytes/(10^9)));
    
else
    fprintf('\n....loading existing grp_rsa_add patterns\n')
    load(g_fname_add);%grp_rsa_add.subj{xsub}
    
    g_fsize = dir(g_fname_add);
    fprintf('file size: %s GB\n', num2str(g_fsize.bytes/(10^9)));
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%============= 2ND LEVEL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(args.cluster, 'local') 
    %% ============= WEIGHT MEAN & DECODABILITY
    args_rsa            = args;
    args_rsa.item_names = item_names;
    args_rsa.basename   = basename;
    
    grp_analysis_rsa_additional(grp_rsa_add, args_rsa, dirs);
    
    %% ============= TIMECOURSE
    grp_rsa = grp_analysis_timecourse_rsa(grp_rsa, args, dirs);
    grp_analysis_timecourse_rsa_sorted(grp_rsa, args, dirs);
%     grp_analysis_timecourse_rsa_sorted_selected_cate(grp_rsa, args, dirs);
 
    %% ============= SAVE GROUP RSA
    g_fname = fullfile(dirs.rsa.group.pattern, sprintf('group_%s.mat', basename));
    
    fprintf('\n....saving grp_rsa patterns\n')
    save(g_fname, 'grp_rsa','-v7.3')
    
    g_fsize = dir(g_fname);
    fprintf('file size: %s GB\n', num2str(g_fsize.bytes/(10^9)));
 
end    
end
