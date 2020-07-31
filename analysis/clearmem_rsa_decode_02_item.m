function[] = clearmem_rsa_decode_02_item(args, dirs)
%*************** create template by averaging the extracted pattern
%*************** category level
%*************** SPM / MVPA PRINCETON TOOLBOX
if strcmp(args.xcluster, 'local')
    addpath ~/github/lewpealab/mvpa/
    addpath(genpath('~/github/GroupICATv4.0b/icatb/'))
    
elseif strcmp(args.xcluster, 'blanca')
    %*************** MVPA PRINCETON TOOLBOX
    addpath('/work/ics/data/projects/banichlab/studies/wmem/fmri/mvpa/utaustin/cu_src/mvpa')
    
    %*************** ICA TOOLBOX
    addpath(genpath('/work/ics/data/projects/banichlab/studies/wmem/fmri/mvpa/utaustin/cu_src/GroupICATv4.0b/icatb'))    
end
mvpa_add_paths;

%% ============= SETUP DIRECTORY
spm_mask_dir    = fullfile(dirs.rsa.roi.home, 'spm');
spm_dir         = dirs.rsa.roi.spm;
output_dir      = dirs.rsa.pattern;

%% ============= SETUP PARAMETERS
xparam          = args.index{1}.param;
category_names  = {'face','fruit','scene'};
n_category      = xparam.n_category;
n_subcategory   = xparam.n_subcategory;
n_item          = xparam.n_item;
n_items         = n_category * n_subcategory * xparam.n_item;
item_array      = 1:n_items;
cate_array      = 1:n_category;
conds_names     = args.index{2}.param.conds_names;
n_condition     = length(conds_names);
n_trials        = args.index{2}.param.n_trials;%of study
n_trs           = args.tc_tr;

spm_p_thresh    = args.spm_p_thresh;
spm_v_extent    = args.spm_v_extent;
spm_thresh_type = args.spm_thresh_type;

rsa_mask        = args.rsa_mask;

for xcate = 1:xparam.n_category
    for xsubcate = 1:xparam.n_subcategory
        for xitem = 1:xparam.n_item
            xunit = xitem + (xparam.n_item * (xsubcate - 1)) + ...
                ((xparam.n_item * xparam.n_subcategory) * (xcate - 1));
            item_names{xunit} = xparam.item_name{xcate}{xsubcate}{xitem}(1:end-4);  %#ok<*AGROW>
        end
    end
end

%% ============= LOADING EXISTING FILE
if args.blk
    basename = sprintf('rsa_%s_spmT_%s_thresh%s_ext%s_%s_mask', ...
        args.level, spm_thresh_type, num2str(spm_p_thresh), num2str(spm_v_extent), rsa_mask);
else% for localizer: TRs when stimulus is presented for 3TRs
    basename = sprintf('rsa_3tr_%s_spmT_%s_thresh%s_ext%s_%s_mask', ...
        args.level, spm_thresh_type, num2str(spm_p_thresh), num2str(spm_v_extent), rsa_mask);
end

if args.item_within
    rsa_fname = fullfile(dirs.rsa.pattern, sprintf('%s_within.mat', basename));
else
    rsa_fname = fullfile(dirs.rsa.pattern, sprintf('%s.mat', basename));
end

%% ============= LOGGING
fprintf('(+) using scratch dir: %s\n', output_dir);
%*************** turn on diary to capture analysis output
diary off;
diary(sprintf('%s/diary.txt', output_dir));
fprintf('running code: %s at %s\n', mfilename, datestr(now, 0))
fprintf('#####################################################################\n\n');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%============= EXTRACT PATTERNS FROM TSPM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('#####################################################################\n');
fprintf('=========== Extract pattern for RSA \n');
fprintf('#####################################################################\n');

%% ============= 01: LOCALIZER EPI PATTERNS
%*************** load mask + read in epis
%*************** zsoring epis
%*************** option: read epis in mask | wholebrain
%*************** define: args.wholebrain '0 | 1'
xph = 1; args.xphase = xph;

xmatrix = args.regs{xph}.matrix;
xheader = args.regs{xph}.header;
n_runs  = args.index{xph}.param.n_runs;

g_runs   = 5;
g_items  = 6;
loc_runs = unique(xmatrix(findCol(xheader, {'it_run'}), :));

for xmaskcate = 1:n_category
    %% ============= INITIATE SUBJ. READ MASK FROM TSPM
    xmask_name = sprintf('spmT_000%d_%s_thresh%s_ext%s_%s_mask',...
        xmaskcate + 1, spm_thresh_type, num2str(spm_p_thresh), num2str(spm_v_extent), rsa_mask);
    
    ph1.basename = sprintf('item_%s_%s_%s_zscored_%s', args.phase_name{xph}, category_names{xmaskcate}, xmask_name, args.epi_name);
    if args.item_within
        ph1.basename = sprintf('item_%s_within_%s_%s_zscored_%s', args.phase_name{xph}, category_names{xmaskcate}, xmask_name, args.epi_name);
    end
    fname        = sprintf('%s/ph1_%s.mat', dirs.rsa.scratch, ph1.basename);
    
    if args.load_pattern_rsa{xph}
        %% ============= LOAD PATTERNS UNDER CATEGORY MASK
        clear tt subj
        
        t_ph1.basename = sprintf('%s_%s_%s_zscored_%s', args.phase_name{xph}, category_names{xmaskcate}, xmask_name, args.epi_name);
        t_fname        = sprintf('%s/ph1_%s.mat', dirs.rsa.scratch, t_ph1.basename);
        
        fprintf('  ... load patterns from ph1\n');
        tt   = load(t_fname);
        subj = tt.ph1.subj; nm = tt.ph1.nm;

        %*************** reset subj.
        subj = remove_object(subj,'pattern', ...
            sprintf('%s_%s_patterns', args.phase_name{xph}, category_names{xmaskcate}));
        subj = remove_object(subj,'pattern','spm_beta');
        
        %*************** EPI PATTERN: subj.patterns{1}
        mask_obj = get_object(subj, 'mask', xmask_name); %'mask object' contents
        fprintf('\n(+) load epi data under mask %s with %s voxels\n', ...
            xmask_name, num2str(mask_obj.nvox));
        
        summarize(subj, 'objtype', 'pattern');
        
        %% ============= EXTRACTING WEIGHT(BETA) FROM SPM
        for xitem = 1:n_items
            subj = load_spm_pattern(subj, sprintf('spm_beta_%s', item_names{xitem}), xmask_name, ...
                fullfile(dirs.rsa.roi.spm, sprintf('spmT_%04i.%s', xitem, args.epiext)));
        end
        
        summarize(subj, 'objtype', 'pattern');
        
        %% ============= save phase 1.
        ph1.args = args; ph1.nm = nm; ph1.subj = subj;
        save(fname,'ph1','-v7.3');
        
        fsize = dir(fname);
        fprintf('file size: %s GB\n', num2str(fsize.bytes/(10^9)));
        
    else
        fprintf('  ... load patterns from ph1\n'); load(fname);
        subj = ph1.subj; nm = ph1.nm;
        
        fsize = dir(fname);
        fprintf('file size: %s GB\n', num2str(fsize.bytes/(10^9)));
    end
    
    summarize(subj)
    
    xsubj{xph}{xmaskcate} = subj;
    
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%============= EXTRACT PATTERNS FROM LOCALIZER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xph           = 1;
loc_n_runs    = 5;
loc_run_array = 1:loc_n_runs;

if args.rsa_pattern(xph)
    for xcate = 1:n_category
        for xsubcate = 1:n_subcategory
            for xrun = 1:g_runs
                rsa.item.cate{xcate}.subcate{xsubcate}.corr   = cell(g_runs * g_items, g_runs * g_items);
                rsa.item.cate{xcate}.subcate{xsubcate}.corr_w = cell(g_runs * g_items, g_runs * g_items);
            end
        end
        
        rsa.template.maskcate{xcate}.corr   = cell(n_items, n_items);
        rsa.template.maskcate{xcate}.corr_w = cell(n_items, n_items);
    end
    
    for xmaskcate = 1:n_category
        for xcate = 1:n_category
            for xsubcate = 1:n_subcategory
                for xitem = 1:g_items
                    for xrun = loc_run_array
                        pat{xph}.maskcate{xmaskcate}.cate{xcate}.subcate{xsubcate}.item{xitem}.run{xrun}   = [];
                        pat_w{xph}.maskcate{xmaskcate}.cate{xcate}.subcate{xsubcate}.item{xitem}.run{xrun} = [];
                    end
                end
            end
        end
    end
    
    %% ============= BETA EXTRACTION
    for xmaskcate = 1:n_category
        for xitem = 1:(n_category * n_subcategory * g_items)
            xbeta = get_mat(xsubj{xph}{xmaskcate}, 'pattern', sprintf('spm_beta_%s', item_names{xitem}));
            
            item_beta{xmaskcate}{xitem}.name = item_names{xitem};
            item_beta{xmaskcate}{xitem}.mat  = xbeta;
        end
    end
    
    %% ============= ITEM WITHIN SUBCATEGORY
    %*************** TR shifted
    for xmaskcate = 1:n_category
        
        t_vox(xmaskcate) = size(xsubj{xph}{xmaskcate}.patterns{1}.mat, 1);
        
        for xcate = 1:n_category
            for xsubcate = 1:n_subcategory
                for xrun = loc_runs
                    
                    %*************** localizer zscored raw patterns
                    tunit = getDATA(xmatrix', xheader, ...
                        {'it_run','category','subcategory'}, ...
                        {xrun, xcate, xsubcate});
                    
                    xit_items = unique(xmatrix(findCol(xheader, {'item'}), tunit));
                    
                    for xitem = xit_items
                        if args.blk
                            xunit = find(getDATA(xmatrix', xheader, ...
                                {'it_run','category','subcategory','item','spike'}, ...
                                {xrun, xcate, xsubcate, xitem, 0}));
                        else % only 3 TR
                            xunit = find(getDATA(xmatrix', xheader, ...
                                {'it_run','category','subcategory','item','stimulus','spike'}, ...
                                {xrun, xcate, xsubcate, xitem, 1, 0}));
                        end
                        
                        xunit_sh = xunit + args.shift_TRs;
                        
                        if xunit_sh(end) > size(xsubj{xph}{xmaskcate}.patterns{1}.mat, 2)
                            endTR    = size(xsubj{xph}{xmaskcate}.patterns{1}.mat, 2);
                            endunit  = find(xunit_sh == endTR);
                            xunit_sh = xunit_sh(1:endunit);
                        end
                        
                        if ~isempty(xunit_sh)
                            clear xpat xxpat xpat_w
                            
                            xpat  = xsubj{xph}{xmaskcate}.patterns{1}.mat(:,xunit_sh);
                            xcell = xitem + (g_items * (xsubcate-1)) + (g_items * n_subcategory * (xcate-1));
                            xbeta = get_mat(xsubj{xph}{xmaskcate}, 'pattern', sprintf('spm_beta_%s', item_names{xcell}));
                            
                            %*************** weighted with beta
                            for i = 1:size(xpat,2)
                                xxpat       = xpat(:,i)';
                                xpat_w(:,i) = xxpat .* xbeta';
                            end
                            
                            pat{xph}.maskcate{xmaskcate}.cate{xcate}.subcate{xsubcate}.item{xitem}.run{xrun}   = mean(xpat, 2);
                            pat_w{xph}.maskcate{xmaskcate}.cate{xcate}.subcate{xsubcate}.item{xitem}.run{xrun} = mean(xpat_w, 2);
                        end
                    end
                end
            end
        end
    end
    
    %% *************** RSA
    it_items = unique(xmatrix(findCol(xheader, {'item'}), :));
    it_items = it_items(~ismember(it_items, 0));
    
    for xcate = 1:n_category
        for xsubcate = 1:n_subcategory
            for xitem = it_items
                for xrun = loc_runs
                    
                    xcell    = xrun + (g_runs * (xitem-1));
                    
                    xx_pat   = pat{xph}.maskcate{xcate}.cate{xcate}.subcate{xsubcate}.item{xitem}.run{xrun};
                    xx_pat_w = pat_w{xph}.maskcate{xcate}.cate{xcate}.subcate{xsubcate}.item{xitem}.run{xrun};
                    
                    if ~isempty(xx_pat)
                        for yitem = it_items
                            for yrun = loc_runs
                                
                                ycell    = yrun + (g_runs * (yitem-1));
                                
                                yy_pat   = pat{xph}.maskcate{xcate}.cate{xcate}.subcate{xsubcate}.item{yitem}.run{yrun};
                                yy_pat_w = pat_w{xph}.maskcate{xcate}.cate{xcate}.subcate{xsubcate}.item{yitem}.run{yrun};
                                
                                if ~isempty(yy_pat)
                                    %*************** similarity: z-transformed pearson correlation
                                    xr = corr2(xx_pat, yy_pat); xwr_rep = corr2(xx_pat_w, yy_pat_w);
                                    
                                    rsa.item.cate{xcate}.subcate{xsubcate}.corr{xcell, ycell} = ...
                                        .5*(log(1+xr) - log(1-xr));
                                    rsa.item.cate{xcate}.subcate{xsubcate}.corr_w{xcell, ycell} = ...
                                        .5*(log(1+xwr_rep) - log(1-xwr_rep));
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    
    %% ============= ITEM ACROSS CATEGORY: targ vs. non-targ
    %*************** setup structure
    % targ:1_targ, 2_related (category), 3_nontarg
    xfid = fopen(fullfile(output_dir, 'rsa_items_vs_items.txt'), 'w+');
    fprintf(xfid, 'RSA_items_vs_items\n\n');
    
    for xcate = 1:n_category
        for xsubcate = 1:n_subcategory
            for xtarg = 1:3
                rsa.item.cate{xcate}.subcate{xsubcate}.targ{xtarg} = [];
            end
        end
    end
    
    %*************** find items used
    clear t_items it_items
    all_items = [];
    for xcate = 1:n_category
        for xsubcate = 1:n_subcategory
            t_items  = unique(getDATA(xmatrix', xheader, ...
                {'category','subcategory'}, {xcate, xsubcate}, ...
                findCol(xheader, {'item'})))';
            all_items = horzcat(all_items, ...
                t_items + (g_items * (xsubcate-1)) + ((g_items * n_subcategory) * (xcate-1)));
            
        end
    end
    
    %*************** RSA    
    for xitem = all_items
        clear xx_pat_w yy_pat_w
        xcate    = fix((xitem-1)/(g_items * n_subcategory)) + 1;
        xsubcate = fix(mod(xitem-1, (g_items * n_subcategory))/g_items) + 1;
        it_xitem = mod(xitem-1, 6) + 1;
        
        for xrun = loc_runs
            %*************** xpat
            xx_pat_w = pat_w{xph}.maskcate{xcate}.cate{xcate}.subcate{xsubcate}.item{it_xitem}.run{xrun};
            
            if ~isempty(xx_pat_w)
                
                fprintf(xfid, '##############################\n');
                fprintf(xfid, 'xitem: %s: run: %s, cate %d, subcate %d\n', num2str(xitem), xrun, xcate, xsubcate);
                
                for yitem = all_items
                    ycate    = fix((yitem-1)/(g_items * n_subcategory)) + 1;
                    ysubcate = fix(mod(yitem-1, (g_items * n_subcategory))/g_items)+1;
                    it_yitem = mod(yitem-1, 6) + 1;
                    
                    for yrun = loc_runs
                        %*************** ypat
                        yy_pat_w = pat_w{xph}.maskcate{xcate}.cate{ycate}.subcate{ysubcate}.item{it_yitem}.run{yrun};
                        
                        if ~isempty(yy_pat_w)
                            %*************** target
                            xtarg = 0;
                            if xitem == yitem
                                xtarg = 1;
                            else
                                if xcate==ycate
                                    xtarg = 2;
                                else
                                    xtarg = 3;
                                end
                            end
                            
                            fprintf(xfid, '... yitem: %s: run: %s, cate %d, subcate %d\n', num2str(yitem), yrun, ycate, ysubcate);
                            fprintf(xfid, '....... targ: %d\n', xtarg);
                            
                            %*************** similarity: z-transformed pearson correlation
                            xwr_rep   = corr2(xx_pat_w, yy_pat_w);
                            z_xwr = .5*(log(1+xwr_rep) - log(1-xwr_rep));
                            
                            if ~isinf(z_xwr)
                                rsa.item.cate{xcate}.subcate{xsubcate}.targ{xtarg} = ...
                                    horzcat(rsa.item.cate{xcate}.subcate{xsubcate}.targ{xtarg}, ...
                                    .5*(log(1+xwr_rep) - log(1-xwr_rep)));
                            end
                        else
                            fprintf(xfid, '** empty: yrun %d\n', yrun);
                        end
                    end
                end
            else
                fprintf(xfid, '** empty: xrun %d\n', xrun);
            end
        end
    end
    
    fclose(xfid);
    
    %% ============= ITEM TEMPLATE ACROSS CATEGORY
    it_items      = unique(xmatrix(findCol(xheader, {'item'}), :));
    it_items      = it_items(~ismember(it_items, 0));
    
    for xmaskcate = 1:n_category
        for xcate = 1:n_category
            for xsubcate = 1:n_subcategory
                for xitem = it_items
                    clear xpat xpat_w
                    xpat = []; xpat_w = [];
                    
                    for xrun = loc_runs
                        xpat   = horzcat(xpat, pat{xph}.maskcate{xmaskcate}.cate{xcate}.subcate{xsubcate}.item{xitem}.run{xrun});
                        xpat_w = horzcat(xpat_w, pat_w{xph}.maskcate{xmaskcate}.cate{xcate}.subcate{xsubcate}.item{xitem}.run{xrun});
                    end
                    
                    template.maskcate{xmaskcate}.cate{xcate}.subcate{xsubcate}.item{xitem}.pat   = mean(xpat, 2);
                    template.maskcate{xmaskcate}.cate{xcate}.subcate{xsubcate}.item{xitem}.pat_w = mean(xpat_w, 2);
                    
                    %*************** different repetition
                    for xrep = loc_run_array
                        clear xpat_w
                        xpat_w = [];
                        
                        for xrun = loc_run_array(1:xrep)
                            
                            xx = pat_w{xph}.maskcate{xmaskcate}.cate{xcate}.subcate{xsubcate}.item{xitem}.run{xrun};
                            
                            if ~isempty(xx)
                                xpat_w = horzcat(xpat_w, xx);
                            end
                        end
                        
                        % mean repetition
                        if ~isempty(xpat_w)
                            template.maskcate{xmaskcate}.cate{xcate}.subcate{xsubcate}.item{xitem}.rep{xrep}.pat_w = mean(xpat_w, 2);
                        end
                        
                        % each repetition
                        template.maskcate{xmaskcate}.cate{xcate}.subcate{xsubcate}.item{xitem}.run{xrep}.pat_w = [];
                        
                        xx = pat_w{xph}.maskcate{xmaskcate}.cate{xcate}.subcate{xsubcate}.item{xitem}.run{xrep};
                        if ~isempty(xx)
                            template.maskcate{xmaskcate}.cate{xcate}.subcate{xsubcate}.item{xitem}.run{xrep}.pat_w = xx;
                        end
                    end
                end
            end
        end
    end
    
    for i = 1:n_items
        xcate    = fix((i-1)/(g_items * n_subcategory)) + 1;
        xsubcate = mod(fix((i-1)/g_items), n_subcategory) + 1;
        xitem    = mod((i-1), g_items) + 1;
        
        template.items{i}.pat   = template.maskcate{xcate}.cate{xcate}.subcate{xsubcate}.item{xitem}.pat;
        template.items{i}.pat_w = template.maskcate{xcate}.cate{xcate}.subcate{xsubcate}.item{xitem}.pat_w;
        
        %*************** different repetition/run
        for xrun = loc_run_array
            xx = template.maskcate{xcate}.cate{xcate}.subcate{xsubcate}.item{xitem}.rep{xrun};
            if ~isempty(xx)
                template.items{i}.rep{xrun}.pat_w = xx.pat_w;
            end
            
            xx = template.maskcate{xcate}.cate{xcate}.subcate{xsubcate}.item{xitem}.run{xrun};
            if ~isempty(xx)
                template.items{i}.run{xrun}.pat_w = xx.pat_w;
            end
        end
    end
    
    %% *************** RSA
    for xmaskcate = 1:n_category
        for xcate = 1:n_category
            for xsubcate = 1:n_subcategory
                for xitem = it_items
                    
                    xcell    = xitem + (g_items * (xsubcate-1)) + (g_items * n_subcategory * (xcate-1));
                    
                    xx_pat   = template.maskcate{xmaskcate}.cate{xcate}.subcate{xsubcate}.item{xitem}.pat;
                    xx_pat_w = template.maskcate{xmaskcate}.cate{xcate}.subcate{xsubcate}.item{xitem}.pat_w;
                    
                    if ~isempty(xx_pat)
                        for ycate = 1:n_category
                            for ysubcate = 1:n_subcategory
                                for yitem = it_items
                                    
                                    ycell    = yitem + (g_items * (ysubcate-1)) + (g_items * n_subcategory * (ycate-1));
                                    
                                    yy_pat   = template.maskcate{xmaskcate}.cate{ycate}.subcate{ysubcate}.item{yitem}.pat;
                                    yy_pat_w = template.maskcate{xmaskcate}.cate{ycate}.subcate{ysubcate}.item{yitem}.pat_w;
                                    
                                    if ~isempty(yy_pat)
                                        %*************** similarity: z-transformed pearson correlation
                                        xr = corr2(xx_pat, yy_pat); xwr_rep = corr2(xx_pat_w, yy_pat_w);
                                        
                                        rsa.template.maskcate{xmaskcate}.corr{xcell, ycell} = ...
                                            .5*(log(1+xr) - log(1-xr));
                                        rsa.template.maskcate{xmaskcate}.corr_w{xcell, ycell} = ...
                                            .5*(log(1+xwr_rep) - log(1-xwr_rep));
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    
    %% *************** RSA
    for xmaskcate = 1:n_category
        for xrun = loc_runs
            for xcate = 1:n_category
                for xsubcate = 1:n_subcategory
                    for xitem = it_items
                        clear xx_pat_w
                        xcell    = xitem + (g_items * (xsubcate-1)) + (g_items * n_subcategory * (xcate-1));
                        xx_pat_w = pat_w{xph}.maskcate{xmaskcate}.cate{xcate}.subcate{xsubcate}.item{xitem}.run{xrun};
                        
                        if ~isempty(xx_pat_w)
                            for ycate = 1:n_category
                                for ysubcate = 1:n_subcategory
                                    for yitem = it_items
                                        clear yy_pat_w
                                        ycell    = yitem + (g_items * (ysubcate-1)) + (g_items * n_subcategory * (ycate-1));
                                        yy_pat_w = pat_w{xph}.maskcate{xmaskcate}.cate{ycate}.subcate{ysubcate}.item{yitem}.run{xrun};
                                        
                                        if ~isempty(yy_pat_w)
                                            %*************** similarity: z-transformed pearson correlation
                                            xwr_rep = corr2(xx_pat_w, yy_pat_w);
                                            
                                            rsa.template.maskcate{xmaskcate}.run{xrun}.corr_w{xcell, ycell} = ...
                                                .5*(log(1+xwr_rep) - log(1-xwr_rep));
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    
    %% *************** subj_rsa
    subj_rsa{xph}.rsa         = rsa;
    subj_rsa{xph}.item_beta   = item_beta;
    subj_rsa{xph}.template    = template;
    subj_rsa{xph}.info.nVox   = t_vox;
    subj_rsa{xph}.info.nTRs   = size(subj.patterns{1}.mat, 2);
    subj_rsa{xph}.info.matrix = xmatrix;
    subj_rsa{xph}.info.header = xheader;
    
else
    %*************** load existing file
    xrsa = load(rsa_fname);%,'subj_rsa'
    
    subj_rsa{xph} = xrsa.subj_rsa{xph};
end

subj_rsa{xph}.subj = xsubj{xph};

%% ============= 02: STUDY EPI PATTERNS
%*************** load mask + read in epis
%*************** zsoring epis
%*************** option: read epis in mask | wholebrain
%*************** define: args.wholebrain '0 | 1'
clear rsa
xph = 2; args.xphase = xph;

xmatrix = args.regs{xph}.matrix;
xheader = args.regs{xph}.header;
n_runs  = args.index{xph}.param.n_runs;

xparam   = args.index{2}.param;
dur_stim = xparam.dur_stim;
dur_mani = xparam.dur_manipulation;
dur_iti  = max(xparam.dur_iti);
dur_sync = n_trs - dur_stim - dur_mani - dur_iti;%next trial

% xtarg: 1_target, 2_related_nontarg 2_nontarg
n_targ = 3;

for xmaskcate = 1:n_category
    %% ============= INITIATE SUBJ. READ MASK FROM TSPM
    xmask_name = sprintf('spmT_000%d_%s_thresh%s_ext%s_%s_mask',...
        xmaskcate + 1, spm_thresh_type, num2str(spm_p_thresh), num2str(spm_v_extent), rsa_mask);
    
    ph2.basename = sprintf('item_%s_%s_%s_zscored_%s', args.phase_name{xph}, category_names{xmaskcate}, xmask_name, args.epi_name);
    fname        = sprintf('%s/ph2_%s.mat', dirs.rsa.scratch, ph2.basename);
    
    if args.load_pattern_rsa{xph}
        %% ============= LOAD PATTERNS UNDER CATEGORY MASK
        clear tt subj
        
        t_ph2.basename = sprintf('%s_%s_%s_zscored_%s', args.phase_name{xph}, category_names{xmaskcate}, xmask_name, args.epi_name);
        t_fname        = sprintf('%s/ph2_%s.mat', dirs.rsa.scratch, t_ph2.basename);
        
        fprintf('  ... load patterns from ph2\n');
        tt   = load(t_fname);
        subj = tt.ph2.subj; nm = tt.ph2.nm;
        
        %*************** reset subj.
        subj = remove_object(subj,'pattern', ...
            sprintf('%s_%s_patterns', args.phase_name{xph}, category_names{xmaskcate}));
        
        %*************** EPI PATTERN: subj.patterns{1}
        mask_obj = get_object(subj, 'mask', xmask_name); %'mask object' contents
        fprintf('\n(+) load epi data under mask %s with %s voxels\n', ...
            xmask_name, num2str(mask_obj.nvox));
        
        summarize(subj, 'objtype', 'pattern');
        
        %% ============= save phase 1.
        ph2.args = args; ph2.nm = nm; ph2.subj = subj;
        save(fname,'ph2','-v7.3');
        
        fsize = dir(fname);
        fprintf('file size: %s GB\n', num2str(fsize.bytes/(10^9)));
        
    else
        fprintf('  ... load patterns from ph1\n'); load(fname);
        subj = ph2.subj; nm = ph2.nm;
        
        fsize = dir(fname);
        fprintf('file size: %s GB\n', num2str(fsize.bytes/(10^9)));
    end
    
    summarize(subj)
    
    xsubj{xph}{xmaskcate} = subj;
    
end

subj_rsa{xph}.info.nTRs   = size(subj.patterns{1}.mat, 2);
subj_rsa{xph}.info.matrix = xmatrix;
subj_rsa{xph}.info.header = xheader; %#ok<*NASGU>

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%============= EXTRACT PATTERNS FROM STUDY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% xtarg: 1_target, 2_related_nontarg 2_nontarg
tic
xph = 2;

if args.rsa_pattern(xph)
    for xcond = 1:n_condition
        for xtarg = 1:n_targ
            for xtr = 1:n_trs
                timecourse.cond{xcond}.targ{xtarg}.tr{xtr}.corr_w = [];
                
                for xcate = 1:n_category
                    for xsubcate = 1:n_subcategory
                        timecourse.cond{xcond}.cate{xcate}.subcate{xsubcate}.targ{xtarg}.tr{xtr}.corr_w = [];
                    end
                end
            end
            
            for xtr = 1:dur_sync
                timecourse.sync.cond{xcond}.targ{xtarg}.tr{xtr}.corr_w = [];
            end
        end
    end
    
    %% *************** STUDY: timecourse: condition
    for xrun = 1:n_runs
        
        fprintf('run: %d | trial: ', xrun)
        
        for xtrial = 1:n_trials(xrun)
            
            fprintf('%s.', num2str(xtrial))
            
            %*************** extract params
            clear item_id xitem_pat xitem_pat_w
            
            tunit      = getDATA(xmatrix', xheader, ...
                {'run', 'trial','presentation'}, {xrun, xtrial, 1});
            xcond      = unique(xmatrix(findCol(xheader, {'condition'}), tunit));
            xcate      = unique(xmatrix(findCol(xheader, {'category'}), tunit));
            xsubcate   = unique(xmatrix(findCol(xheader, {'subcategory'}), tunit));
            xitem      = unique(xmatrix(findCol(xheader, {'item'}), tunit));
            related_items = (1:n_item * n_subcategory) + (n_item * n_subcategory * (xcate - 1));
            
            %*************** 1_target, 2_nontarget, 3_baseline
            item_id{1} = xitem + (g_items * (xsubcate-1)) + (g_items * n_subcategory * (xcate-1));
            item_id{2} = related_items(~ismember(related_items, item_id{1}));
            item_id{3} = item_array(~ismember(item_array, related_items));
            
            %*************** template
            for xtarg = 1:3%1_target, 2_nontarget, 3_baseline
                for i = 1:length(item_id{xtarg})
                    xit = item_id{xtarg}(i);
                    xitem_pat_w{xtarg}{i} = subj_rsa{1}.template.items{xit}.pat_w;
                end
            end
            
            %*************** timecourse
            xunit     = find(getDATA(xmatrix', xheader, {'run', 'trial'}, {xrun, xtrial}));
            xunit_tc  = xunit(1):(xunit(1) + n_trs - 1);
            xspike    = ~(xmatrix(findCol(xheader, {'spike'}), xunit_tc));
            
            if sum(xspike)~=0
                for xtr = 1:length(xunit_tc)
                    if xspike(xtr)
                        for xtarg = 1:3%1_target, 2_nontarget, 3_baseline
                            
                            t_rep_corr = [];
                            
                            for i = 1:length(item_id{xtarg})
                                clear ytarg_pat ytarg_pat_w
                                %*************** target
                                xit         = item_id{xtarg}(i);
                                ycate       = fix((xit-1)/(g_items * n_subcategory)) + 1;
                                
                                if (xtarg==1) && (xcate~=ycate)
                                    fprintf('warning: not matching category for target item\n')
                                end
                                
                                xbeta       = xsubj{1}{ycate}.patterns{xit + 1}.mat;
                                ytarg_pat   = xsubj{2}{ycate}.patterns{1}.mat(:, xunit_tc(xtr));
                                ytarg_pat_w = (ytarg_pat' .* xbeta')';
                                
                                %*************** similarity: z-transformed pearson correlation
                                t_rep_corr = horzcat(t_rep_corr, corr2(xitem_pat_w{xtarg}{i}, ytarg_pat_w));
                            end
                            
                            xwr_rep = mean(t_rep_corr);
                            
                            timecourse.cond{xcond}.targ{xtarg}.tr{xtr}.corr_w = ...
                                horzcat(timecourse.cond{xcond}.targ{xtarg}.tr{xtr}.corr_w,...
                                .5*(log(1+xwr_rep) - log(1-xwr_rep)));
                            
                            timecourse.cond{xcond}.cate{xcate}.subcate{xsubcate}.targ{xtarg}.tr{xtr}.corr_w = ...
                                horzcat(timecourse.cond{xcond}.cate{xcate}.subcate{xsubcate}.targ{xtarg}.tr{xtr}.corr_w,...
                                .5*(log(1+xwr_rep) - log(1-xwr_rep)));
                            
                            if xtr > length(xunit)
                                it_tr = xtr - length(xunit);
                                if it_tr <= dur_sync
                                    
                                    timecourse.sync.cond{xcond}.targ{xtarg}.tr{it_tr}.corr_w = ...
                                        horzcat(timecourse.sync.cond{xcond}.targ{xtarg}.tr{it_tr}.corr_w,...
                                        .5*(log(1+xwr_rep) - log(1-xwr_rep)));
                                end
                            end
                        end
                    end
                end
            end
        end
        
        fprintf('\n')
    end
    
    %% *************** condition: timecourse: mean per TR
    for xcond = 1:n_condition
        for xtarg = 1:n_targ
            for xtr = 1:n_trs
                rsa.timecourse.cond{xcond}.targ{xtarg}.tr{xtr}.corr_w = ...
                    mean(timecourse.cond{xcond}.targ{xtarg}.tr{xtr}.corr_w);
            end
            
            %*************** sync
            for xtr = 1:dur_sync
                rsa.timecourse.sync.cond{xcond}.targ{xtarg}.tr{xtr}.corr_w = ...
                    mean(timecourse.sync.cond{xcond}.targ{xtarg}.tr{xtr}.corr_w);
            end
        end
    end
    
    %% *************** SETUP STUDY: N=N+1: timecourse: condition
    % same for N, N+1: timecourse.sameseq{xlevel}.cond{xcond}.targ{xtarg}.tr{xtr}.corr
    % xlevel: 1_cate, 2_subcate
    
    n_targ   = 3;
    
    for xcond = 1:n_condition
        for xtarg = 1:n_targ
            for xtr = 1:n_trs
                for xlevel = 1:2%1_cate, 2_subcate
                    
                    if xlevel==1, it_sames = 1:2; else, it_sames = 1:3; end
                    
                    for xsame = it_sames % 1_same, 2_differ, 3_related
                        timecourse.sameseq{xlevel}{xsame}.cond{xcond}.targ{xtarg}.tr{xtr}.corr_w = [];
                        
                        if xtr <= dur_sync
                            timecourse.sameseq{xlevel}{xsame}.sync.cond{xcond}.targ{xtarg}.tr{xtr}.corr_w = [];
                            timecourse.sameseq{xlevel}{xsame}.sync.next.cond{xcond}.targ{xtarg}.tr{xtr}.corr_w = [];
                        end
                    end
                end
            end
        end
    end
    
    %% *************** STUDY: N=N+1: timecourse: condition
    % same for N, N+1: timecourse.sameseq{xlevel}{xsame}.cond{xcond}.targ{xtarg}.tr{xtr}.corr
    % xlevel: 1_cate, 2_subcate
    
    for xrun = 1:n_runs
        
        fprintf('run: %d | trial: ', xrun)
        
        for xtrial = 1:n_trials(xrun)-1
            
            fprintf('%s.', num2str(xtrial))
            
            %*************** extract params
            clear item_id xitem_pat xitem_pat_w
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%
            %*************** N trial
            tunit       = getDATA(xmatrix', xheader, ...
                {'run','trial','presentation'}, {xrun, xtrial, 1});
            xcate       = unique(xmatrix(findCol(xheader, {'category'}), tunit));
            xsubcate    = unique(xmatrix(findCol(xheader, {'subcategory'}), tunit));
            xitem       = unique(xmatrix(findCol(xheader, {'item'}), tunit));
            related_items = (1:n_item * n_subcategory) + (n_item * n_subcategory * (xcate - 1));

            %*************** from N trial
            xcond       = unique(xmatrix(findCol(xheader, {'condition'}), tunit));
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%
            %*************** N+1 trial
            ttunit      = getDATA(xmatrix', xheader, ...
                {'run','trial','presentation'}, {xrun, xtrial+1, 1});
            xxcate      = unique(xmatrix(findCol(xheader, {'category'}), ttunit));
            xxsubcate   = unique(xmatrix(findCol(xheader, {'subcategory'}), ttunit));
            xxitem      = unique(xmatrix(findCol(xheader, {'item'}), ttunit));
            ncate_items = (1:n_item * n_subcategory) + (n_item * n_subcategory * (xxcate - 1));
            
            %%%%%%%%%%%%%%%%%%%%%%%%%% N trial
            %*************** 1_target, 2_nontarget, 3_baseline
            item_id{1}  = xitem + (g_items * (xsubcate-1)) + (g_items * n_subcategory * (xcate-1));
            item_id{2}  = related_items(~ismember(related_items, item_id{1}));
            item_id{3}  = item_array(~ismember(item_array, related_items));
            
            %*************** template
            for xtarg = 1:3%1_target, 2_nontarget, 3_baseline
                for i = 1:length(item_id{xtarg})
                    xit = item_id{xtarg}(i);
                    xitem_pat_w{xtarg}{i} = subj_rsa{1}.template.items{xit}.pat_w;
                end
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%% N+1 trial
            %*************** 1_target, 2_nontarget, 3_baseline
            nitem_id{1} = xxitem + (g_items * (xxsubcate-1)) + (g_items * n_subcategory * (xxcate-1));
            nitem_id{2} = ncate_items(~ismember(ncate_items, nitem_id{1}));
            nitem_id{3} = item_array(~ismember(item_array, ncate_items));
            
            %*************** template
            for xtarg = 1:3%1_target, 2_nontarget, 3_baseline
                for i = 1:length(nitem_id{xtarg})
                    xxit = nitem_id{xtarg}(i);
                    xxitem_pat_w{xtarg}{i} = subj_rsa{1}.template.items{xxit}.pat_w;
                end
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%
            %*************** timecourse
            xunit    = find(getDATA(xmatrix', xheader, {'run', 'trial'}, {xrun, xtrial}));
            xunit_tc = xunit(1):(xunit(1) + n_trs - 1);
            xspike   = ~(xmatrix(findCol(xheader, {'spike'}), xunit_tc));
            
            if sum(xspike)~=0
                for xtr = 1:length(xunit_tc)
                    if xspike(xtr)
                        for xtarg = 1:3%1_target, 2_nontarget, 3_baseline
                            t_rep_corr = []; t_xxwr = [];
                            
                            for i = 1:length(item_id{xtarg})
                                %%%%%%%%%%%%%%%%%%%%%%%%%% N trial
                                clear ytarg_pat ytarg_pat_w
                                %*************** target
                                xit         = item_id{xtarg}(i);
                                ycate       = fix((xit-1)/(g_items * n_subcategory)) + 1;
                                
                                xbeta       = xsubj{1}{ycate}.patterns{xit + 1}.mat;
                                ytarg_pat   = xsubj{2}{ycate}.patterns{1}.mat(:, xunit_tc(xtr));
                                ytarg_pat_w = (ytarg_pat' .* xbeta')';
                                
                                %*************** similarity: z-transformed pearson correlation
                                t_rep_corr = horzcat(t_rep_corr, corr2(xitem_pat_w{xtarg}{i}, ytarg_pat_w));
                                
                                %%%%%%%%%%%%%%%%%%%%%%%%%% N+1 trial
                                clear ytarg_pat ytarg_pat_w
                                %*************** target
                                xxit         = nitem_id{xtarg}(i);
                                yycate       = fix((xxit-1)/(g_items * n_subcategory)) + 1;
                                
                                xxbeta       = xsubj{1}{yycate}.patterns{xxit + 1}.mat;
                                yytarg_pat   = xsubj{2}{yycate}.patterns{1}.mat(:, xunit_tc(xtr));
                                yytarg_pat_w = (yytarg_pat' .* xxbeta')';
                                
                                %*************** similarity: z-transformed pearson correlation
                                t_xxwr = horzcat(t_xxwr, corr2(xxitem_pat_w{xtarg}{i}, yytarg_pat_w));
                            end

                            xwr_rep = mean(t_rep_corr); xxwr = mean(t_xxwr);
                            
                            %*************** same category for N+1 trial
                            % 1_category, 2_subcategory
                            % 1_same, 2_differ, 3_related
                            for xlevel = 1:2
                                
                                if xlevel == 1
                                    if (xcate == xxcate)
                                        xsame = 1;
                                    else
                                        xsame = 2;
                                    end
                                    
                                elseif xlevel == 2
                                    
                                    if (xcate == xxcate) && (xsubcate == xxsubcate)
                                        xsame = 1;
                                    elseif (xcate == xxcate) && (xsubcate ~= xxsubcate)
                                        xsame = 3;
                                    elseif (xcate ~= xxcate)
                                        xsame = 2;
                                    end
                                end
                                
                                timecourse.sameseq{xlevel}{xsame}.cond{xcond}.targ{xtarg}.tr{xtr}.corr_w = ...
                                    horzcat(timecourse.sameseq{xlevel}{xsame}.cond{xcond}.targ{xtarg}.tr{xtr}.corr_w, ...
                                    .5*(log(1+xwr_rep) - log(1-xwr_rep)));
                                
                                if xtr > length(xunit)
                                    it_tr = xtr - length(xunit);
                                    if it_tr <= dur_sync
                                        
                                        timecourse.sameseq{xlevel}{xsame}.sync.cond{xcond}.targ{xtarg}.tr{it_tr}.corr_w = ...
                                            horzcat(timecourse.sameseq{xlevel}{xsame}.sync.cond{xcond}.targ{xtarg}.tr{it_tr}.corr_w,...
                                            .5*(log(1+xwr_rep) - log(1-xwr_rep)));
                                        
                                        %*************** actual N+1 item
                                        timecourse.sameseq{xlevel}{xsame}.sync.next.cond{xcond}.targ{xtarg}.tr{it_tr}.corr_w = ...
                                            horzcat(timecourse.sameseq{xlevel}{xsame}.sync.next.cond{xcond}.targ{xtarg}.tr{it_tr}.corr_w,...
                                            .5*(log(1+xxwr) - log(1-xxwr)));
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end    
        fprintf('\n')
    end
    
    fprintf('\n')
    
    %% *************** condition: timecourse: mean per TR
    for xcond = 1:n_condition
        for xtarg = 1:n_targ
            for xtr = 1:n_trs
                for xlevel = 1:2
                    
                    if xlevel==1, it_sames = 1:2; else it_sames = 1:3; end
                    
                    for xsame = it_sames % 1_same, 2_differ, 3_related
                        
                        rsa.timecourse.sameseq{xlevel}{xsame}.cond{xcond}.targ{xtarg}.tr{xtr}.corr_w = ...
                            mean(timecourse.sameseq{xlevel}{xsame}.cond{xcond}.targ{xtarg}.tr{xtr}.corr_w);
                        
                        %*************** sync
                        if xtr <= dur_sync
                            rsa.timecourse.sameseq{xlevel}{xsame}.sync.cond{xcond}.targ{xtarg}.tr{xtr}.corr_w = ...
                                mean(timecourse.sameseq{xlevel}{xsame}.sync.cond{xcond}.targ{xtarg}.tr{xtr}.corr_w);
                            
                            % N+1
                            rsa.timecourse.sameseq{xlevel}{xsame}.sync.next.cond{xcond}.targ{xtarg}.tr{xtr}.corr_w = ...
                                mean(timecourse.sameseq{xlevel}{xsame}.sync.next.cond{xcond}.targ{xtarg}.tr{xtr}.corr_w);
                        end
                        
                        %*************** n trials
                        if (xtarg == 1) && (xtr == 1)
                            rsa.timecourse.sameseq{xlevel}{xsame}.cond{xcond}.n_trial = ...
                                length(timecourse.sameseq{xlevel}{xsame}.cond{xcond}.targ{xtarg}.tr{xtr}.corr_w);
                        end
                    end
                end
            end
        end
    end
    
    %% ============= REPLACE N+1 SAME/DIFF BASED ON NEW-ITEM
    % N(new_item)=N+1: timecourse.sameseq{xlevel}{xsame}.cond{xcond}.newitem{xnew_same}.targ{xtarg}.tr{xtr}.corr
    % cond2: repCat: xsame=2, xnew_same: same category: new_item=next_item
    
    it_level = 1;%1_cate
    it_conds = 2;%2_repCat
    it_same  = 2;%2_diff
    n_targ   = 3;
    
    for xnew_same = 1:2% 1_same, 2_differ
        for xtarg = 1:n_targ
            for xtr = 1:n_trs
                
                timecourse.sameseq{it_level}{it_same}.cond{it_conds}.newitem{xnew_same}.targ{xtarg}.tr{xtr}.corr_w = [];
                
                if xtr <= dur_sync
                    timecourse.sameseq{it_level}{it_same}.sync.cond{it_conds}.newitem{xnew_same}.targ{xtarg}.tr{xtr}.corr_w = [];
                    timecourse.sameseq{it_level}{it_same}.sync.next.cond{it_conds}.newitem{xnew_same}.targ{xtarg}.tr{xtr}.corr_w = [];
                end
            end
        end
        
        n_newsame{xnew_same} = 0;
    end
    
    for xrun = 1:n_runs
        
        fprintf('run: %d | trial: ', xrun)
        
        for xtrial = 1:n_trials(xrun)-1
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%
            %*************** N trial
            tunit = getDATA(xmatrix', xheader, ...
                {'run','trial','presentation'}, {xrun, xtrial, 1});
            xcond = unique(xmatrix(findCol(xheader, {'condition'}), tunit));
            
            if xcond==it_conds
                
                fprintf('%s.', num2str(xtrial))
                
                %*************** extract params
                clear item_id xitem_pat xitem_pat_w
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%
                %*************** from N trial
                xcate        = unique(xmatrix(findCol(xheader, {'category'}), tunit));
                xsubcate     = unique(xmatrix(findCol(xheader, {'subcategory'}), tunit));
                xitem        = unique(xmatrix(findCol(xheader, {'item'}), tunit));
                related_items   = (1:n_item * n_subcategory) + (n_item * n_subcategory * (xcate - 1));
                %*************** from N trial: new item
                xnew_cate    = getDATA(xmatrix', xheader, ...
                    {'run','trial','presentation'}, {xrun, xtrial, 2},...
                    findCol(xheader, {'new_category'}));
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%
                %*************** N+1 trial
                ttunit      = getDATA(xmatrix', xheader, ...
                    {'run','trial','presentation'}, {xrun, xtrial+1, 1});
                xxcate      = unique(xmatrix(findCol(xheader, {'category'}), ttunit));
                xxsubcate   = unique(xmatrix(findCol(xheader, {'subcategory'}), ttunit));
                xxitem      = unique(xmatrix(findCol(xheader, {'item'}), ttunit));
                ncate_items = (1:n_item * n_subcategory) + (n_item * n_subcategory * (xxcate - 1));
                
                %%%%%%%%%%%%%%%%%%%%%%%%%% N trial
                %*************** 1_target, 2_nontarget, 3_baseline
                item_id{1}  = xitem + (g_items * (xsubcate-1)) + (g_items * n_subcategory * (xcate-1));
                item_id{2}  = related_items(~ismember(related_items, item_id{1}));
                item_id{3}  = item_array(~ismember(item_array, related_items));
                
                %*************** template
                for xtarg = 1:3%1_target, 2_nontarget, 3_baseline
                    for i = 1:length(item_id{xtarg})
                        xit = item_id{xtarg}(i);
                        xitem_pat_w{xtarg}{i} = subj_rsa{1}.template.items{xit}.pat_w;
                    end
                end
                
                %%%%%%%%%%%%%%%%%%%%%%%%%% N+1 trial
                %*************** 1_target, 2_nontarget, 3_baseline
                nitem_id{1} = xxitem + (g_items * (xxsubcate-1)) + (g_items * n_subcategory * (xxcate-1));
                nitem_id{2} = ncate_items(~ismember(ncate_items, nitem_id{1}));
                nitem_id{3} = item_array(~ismember(item_array, ncate_items));
                
                %*************** template
                for xtarg = 1:3%1_target, 2_nontarget, 3_baseline
                    for i = 1:length(nitem_id{xtarg})
                        xxit = nitem_id{xtarg}(i);
                        xxitem_pat_w{xtarg}{i} = subj_rsa{1}.template.items{xxit}.pat_w;
                    end
                end
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%
                %*************** timecourse
                xunit    = find(getDATA(xmatrix', xheader, {'run', 'trial'}, {xrun, xtrial}));
                xunit_tc = xunit(1):(xunit(1) + n_trs - 1);
                xspike   = ~(xmatrix(findCol(xheader, {'spike'}), xunit_tc));
                
                if sum(xspike)~=0
                    for xtr = 1:length(xunit_tc)
                        if xspike(xtr)
                            for xtarg = 1:3%1_target, 2_nontarget, 3_baseline
                                t_rep_corr = []; t_xxwr = [];
                                
                                for i = 1:length(item_id{xtarg})
                                    %%%%%%%%%%%%%%%%%%%%%%%%%% N trial
                                    clear ytarg_pat ytarg_pat_w
                                    %*************** target
                                    xit         = item_id{xtarg}(i);
                                    ycate       = fix((xit-1)/(g_items * n_subcategory)) + 1;
                                    
                                    xbeta       = xsubj{1}{ycate}.patterns{xit + 1}.mat;
                                    ytarg_pat   = xsubj{2}{ycate}.patterns{1}.mat(:, xunit_tc(xtr));
                                    ytarg_pat_w = (ytarg_pat' .* xbeta')';
                                    
                                    %*************** similarity: z-transformed pearson correlation
                                    t_rep_corr = horzcat(t_rep_corr, corr2(xitem_pat_w{xtarg}{i}, ytarg_pat_w));
                                    
                                    %%%%%%%%%%%%%%%%%%%%%%%%%% N+1 trial
                                    clear ytarg_pat ytarg_pat_w
                                    %*************** target
                                    xxit         = nitem_id{xtarg}(i);
                                    yycate       = fix((xxit-1)/(g_items * n_subcategory)) + 1;
                                    
                                    xxbeta       = xsubj{1}{yycate}.patterns{xxit + 1}.mat;
                                    yytarg_pat   = xsubj{2}{yycate}.patterns{1}.mat(:, xunit_tc(xtr));
                                    yytarg_pat_w = (yytarg_pat' .* xxbeta')';
                                    
                                    %*************** similarity: z-transformed pearson correlation
                                    t_xxwr = horzcat(t_xxwr, corr2(xxitem_pat_w{xtarg}{i}, yytarg_pat_w));
                                end
                                
                                xwr_rep = mean(t_rep_corr); xxwr = mean(t_xxwr);
                                
                                %*************** same category for N+1 trial
                                % 1_category, 2_subcategory
                                % 1_same, 2_differ, 3_related
                                for xlevel = it_level
                                    if xlevel == 1
                                        if (xcate == xxcate)
                                            xsame = 1;
                                        else
                                            xsame = 2;
                                        end
                                        
                                    elseif xlevel == 2
                                        
                                        if (xcate == xxcate) && (xsubcate == xxsubcate)
                                            xsame = 1;
                                        elseif (xcate == xxcate) && (xsubcate ~= xxsubcate)
                                            xsame = 3;
                                        elseif (xcate ~= xxcate)
                                            xsame = 2;
                                        end
                                    end
                                    
                                    if xsame==it_same
                                        if xnew_cate==xxcate % same category
                                            xnew_same = 1;
                                        else
                                            xnew_same = 2;
                                        end
                                        
                                        n_newsame{xnew_same} = n_newsame{xnew_same} + 1;
                                        
                                        timecourse.sameseq{xlevel}{xsame}.cond{xcond}.newitem{xnew_same}.targ{xtarg}.tr{xtr}.corr_w = ...
                                            horzcat(timecourse.sameseq{xlevel}{xsame}.cond{xcond}.newitem{xnew_same}.targ{xtarg}.tr{xtr}.corr_w,...
                                            .5*(log(1+xwr_rep) - log(1-xwr_rep)));
                                        
                                        if xtr > length(xunit)
                                            it_tr = xtr - length(xunit);
                                            if it_tr <= dur_sync
                                                
                                                timecourse.sameseq{xlevel}{xsame}.sync.cond{xcond}.newitem{xnew_same}.targ{xtarg}.tr{it_tr}.corr_w = ...
                                                    horzcat(timecourse.sameseq{xlevel}{xsame}.sync.cond{xcond}.newitem{xnew_same}.targ{xtarg}.tr{it_tr}.corr_w,...
                                                    .5*(log(1+xwr_rep) - log(1-xwr_rep)));
                                                
                                                %*************** actual N+1 item
                                                timecourse.sameseq{xlevel}{xsame}.sync.next.cond{xcond}.newitem{xnew_same}.targ{xtarg}.tr{it_tr}.corr_w = ...
                                                    horzcat(timecourse.sameseq{xlevel}{xsame}.sync.next.cond{xcond}.newitem{xnew_same}.targ{xtarg}.tr{it_tr}.corr_w,...
                                                    .5*(log(1+xxwr) - log(1-xxwr)));
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
        fprintf('\n')
    end
    
    fprintf('\n')
    
    %% ============= PREDICTING NEXT ITEM (STRONGREST ITEM)
    %*************** for replace category/subcategory
    %*************** repCat (2): targ:1_item, 2_new_item, 3_related_item, 4_related new_item, 5_others
    %*************** repSub (3): targ:1_item, 2_new_item, 3_related, 4_others
    %*************** timecourse.cond{xcond}.item_targ{xtarg}.tr{xtr}.corr_w
    xph       = 2;
    it_conds  = 2:3;
    item_peak = args.item_peak;%for replace: 1_before, 2_after switch
    n_peak    = item_peak(2, 2) - item_peak(2, 1) + 1;
    
    for xcond = it_conds
        
        if xcond==2, n_targ = 5; else n_targ = 4; end
        
        for xtarg = 1:n_targ
            
            for xtr = 1:n_trs
                timecourse.cond{xcond}.item_targ{xtarg}.tr{xtr}.corr_w = [];
            end
        end
    end
    
    fprintf('(+) Replace: predicted next item\n')
    
    for xrun = 1:n_runs
        
        fprintf('run: %d | trial: ', xrun)
        
        for xtrial = 1:n_trials(xrun)
            
            fprintf('%s.', num2str(xtrial))
            
            %*************** extract params
            clear item_id xitem_pat xitem_pat_w
            
            tunit = getDATA(xmatrix', xheader, ...
                {'run', 'trial','presentation'}, {xrun, xtrial, 1});
            xcond = unique(xmatrix(findCol(xheader, {'condition'}), tunit));
            
            if (xcond == 2) || (xcond == 3)
                
                if xcond==2, n_targ = 5; else n_targ = 4; end %#ok<*SEPEX>
                
                %*************** target item
                xcate         = unique(xmatrix(findCol(xheader, {'category'}), tunit));
                xsubcate      = unique(xmatrix(findCol(xheader, {'subcategory'}), tunit));
                xitem         = unique(xmatrix(findCol(xheader, {'item'}), tunit));
                related_items = (1:n_item * n_subcategory) + (n_item * n_subcategory * (xcate - 1));
                
                %*************** new item
                new_unit      = getDATA(xmatrix', xheader, ...
                    {'run', 'trial','presentation'}, {xrun, xtrial, 2});
                xnew_cate     = unique(xmatrix(findCol(xheader, {'new_category'}), new_unit));
                xnew_subcate  = unique(xmatrix(findCol(xheader, {'new_subcategory'}), new_unit));
                
                if xcond == 2
                    %*********** all items in new category
                    new_item_array = (1:n_item * n_subcategory) + (n_item * n_subcategory * (xnew_cate - 1));
                else
                    %*********** all items in new subcategory
                    new_item_array = (1:n_item) + (n_item * (xnew_subcate - 1)) + ...
                        (n_item * n_subcategory * (xnew_cate - 1));
                end
                
                %*************** timecourse
                xunit = find(getDATA(xmatrix', xheader, {'run', 'trial'}, {xrun, xtrial}));
                xpeak = (xunit(1):xunit(1) + 5) + (item_peak(2, 1) - 1);
                
                %*************** predict new_item
                clear xnew_evidence
                for it = 1:length(new_item_array)
                    clear xtemp xwr twr
                    xnew  = new_item_array(it);
                    xtemp = subj_rsa{1}.template.items{xnew}.pat_w;
                    
                    for it_tr = 1:length(xpeak)
                        xtr = xpeak(it_tr);
                        
                        xbeta      = xsubj{1}{xnew_cate}.patterns{xnew + 1}.mat;
                        xnew_pat   = xsubj{2}{xnew_cate}.patterns{1}.mat(:, xtr);
                        xnew_pat_w = (xnew_pat' .* xbeta')';
                        
                        %*************** similarity: z-transformed pearson correlation
                        twr(it_tr) = corr2(xtemp, xnew_pat_w);
                    end
                    
                    xwr_rep = mean(twr);
                    
                    xnew_evidence(it) = .5 * (log(1+xwr_rep) - log(1-xwr_rep));
                end
                
                [~, xwhich] = max(xnew_evidence);
                
                %*************** repCat (2): targ:1_item, 2_new_item, 3_related_item, 4_related new_item, 5_others
                %*************** repSub (3): targ:1_item, 2_new_item, 3_related, 4_others
                item_id{1} = xitem + (g_items * (xsubcate-1)) + (g_items * n_subcategory * (xcate-1));
                item_id{2} = new_item_array(xwhich);
                
                if xcond == 2
                    %*********** 3_related_item, 4_related new_item, 5_others
                    item_id{3} = related_items(~ismember(related_items, item_id{1}));
                    item_id{4} = new_item_array(~ismember(new_item_array, item_id{2}));
                    item_id{5} = item_array(~ismember(item_array, [related_items new_item_array]));
                elseif xcond == 3
                    %*********** 3_related, 4_others
                    item_id{3} = related_items(~ismember(related_items, [item_id{1} item_id{2}]));
                    item_id{4} = item_array(~ismember(item_array, related_items));
                end
                
                %*************** template
                for xtarg = 1:n_targ
                    for i = 1:length(item_id{xtarg})
                        xit = item_id{xtarg}(i);
                        xitem_pat_w{xtarg}{i} = subj_rsa{1}.template.items{xit}.pat_w;
                    end
                end
                
                %*************** timecourse
                xunit_tc = xunit(1):(xunit(1) + n_trs - 1);
                xspike   = ~(xmatrix(findCol(xheader, {'spike'}), xunit_tc));
                
                if sum(xspike)~=0
                    for xtr = 1:length(xunit_tc)
                        if xspike(xtr)
                            for xtarg = 1:n_targ
                                
                                t_rep_corr = [];
                                
                                for i = 1:length(item_id{xtarg})
                                    clear ytarg_pat ytarg_pat_w
                                    %*************** target
                                    xit         = item_id{xtarg}(i);
                                    ycate       = fix((xit-1)/(g_items * n_subcategory)) + 1;
                                    
                                    if (xtarg==1) && (xcate~=ycate)
                                        fprintf('warning: not matching category for target item\n')
                                    end
                                    
                                    xbeta       = xsubj{1}{ycate}.patterns{xit + 1}.mat;
                                    ytarg_pat   = xsubj{2}{ycate}.patterns{1}.mat(:, xunit_tc(xtr));
                                    ytarg_pat_w = (ytarg_pat' .* xbeta')';
                                    
                                    %*************** similarity: z-transformed pearson correlation
                                    t_rep_corr = horzcat(t_rep_corr, corr2(xitem_pat_w{xtarg}{i}, ytarg_pat_w));
                                end
                                
                                xwr_rep = mean(t_rep_corr);
                                
                                timecourse.cond{xcond}.item_targ{xtarg}.tr{xtr}.corr_w = ...
                                    horzcat(timecourse.cond{xcond}.item_targ{xtarg}.tr{xtr}.corr_w, ...
                                    .5*(log(1+xwr_rep) - log(1-xwr_rep)));
                            end
                        end
                    end
                end
            end
        end
        
        fprintf('\n')
    end
    
    %% ============= REACTIVATION OF N ITEM IN N+1 TRIAL: replace (temporary removal): if target is reactivated in N+1 trial
    xph = 2;
    % conds: 'maintain','repCat','repSubcate','suppress','clear'
    % collect N-based RSA on N+1
    
    % Ma,Su,Cl: react.cond{xcond}.same{xsame}.targ{xtarg}.tr{xtr}.corr_w
    %  (same:N == N+1 category)
    %    1_N-target(1,fa), 3_N+1-target(1,fa)
    %    4_Ns-relatives(16,fa),
    %    6_others(36,fr sc)
    %  (diff:N ~= N+1 category)
    %    1_N-target(1,fa), 3_N+1-target(1,fr)
    %    4_N-relative(17,fa), 5_N+1-relative(17,fr),
    %    6_others(18,sc)
    
    % RepCat.2: react.cond{xcond}.same{xsame}.newsame{xnewsame}.targ{xtarg}.tr{xtr}.corr_w
    %  (same:N == N+1 category: N-new ~= N+1)
    %    1_N-target(1,fa), 2_N-newtarg(18,fr), 3_N+1-target(1,fa),
    %    4_Ns-relatives(16,fa)
    %    6_others(18,sc)
    %  (diff:N ~= N+1 category):
    %    (same:N-new == N+1)
    %      1_N-target(1,fa), 2_N-newtarg(17,fr), 3_N+1-target(1,fr)
    %      4_N-relative(17,fa)
    %      6_others(18,sc)
    %    (diff:N-new ~= N+1)
    %      1_N-target(1,fa), 2_N-newtarg(18,fr), 3_N+1-target(1,sc),
    %      4_N-relative(17,fa), 5_N+1-relative(17,sc)
    
    % RepSubcat.3: react.cond{xcond}.same{xsame}.targ{xtarg}.tr{xtr}.corr_w
    %  (same:N == N+1 category: N-new == N+1)
    %    1_N-target(1,fa), 2_N-newtarg(16,fa), 3_N+1-target(1,fa),
    %    6_others(36,fr sc)
    %  (diff:N ~= N+1 category: N-new ~= N+1)
    %    1_N-target(1,fa), 2_N-newtarg(17,fa), 3_N+1-target(1,fr),
    %    5_N+1-relative(17,fr)
    %    6_others(18,sc)
    
    n_targs = 6;
    
    for xcond = 1:n_condition
        for xtarg = 1:n_targs
            for xsame = 1:2
                for xtr = 1:n_trs
                    if (xcond==2) && (xsame == 2)
                        for xnewsame = 1:2
                            react.cond{xcond}.same{xsame}.newsame{xnewsame}.targ{xtarg}.tr{xtr}.corr_w = [];
                        end
                    else
                        react.cond{xcond}.same{xsame}.targ{xtarg}.tr{xtr}.corr_w = [];
                    end
                end
            end
        end
    end
    
    for xrun = 1:n_runs
        
        fprintf('run: %d | trial: ', xrun)
        
        for xtrial = 1:n_trials(xrun)-1
            clear xunit xcond xcate xsubcate xitem xitem_id xrelatives id_off
            
            fprintf('%s.', num2str(xtrial))
            
            for xtarg = 1:n_targs, xitem_id{xtarg} = []; end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%
            %*************** units:1_N, 2_N-new, 3_N+1
            xunit{1} = getDATA(xmatrix', xheader, ...
                {'run','trial','presentation'}, {xrun, xtrial, 1});
            xunit{2} = getDATA(xmatrix', xheader, ...
                {'run','trial','presentation'}, {xrun, xtrial, 2});
            xunit{3} = getDATA(xmatrix', xheader, ...
                {'run','trial','presentation'}, {xrun, xtrial+1, 1});
            
            xcond = unique(xmatrix(findCol(xheader, {'condition'}), xunit{1}));
            
            %*************** item info
            for xtarg = [1 3] %1_N, 3_N+1 item
                xcate{xtarg}      = unique(xmatrix(findCol(xheader, {'category'}), xunit{xtarg}));
                xsubcate{xtarg}   = unique(xmatrix(findCol(xheader, {'subcategory'}), xunit{xtarg}));
                xitem{xtarg}      = unique(xmatrix(findCol(xheader, {'item'}), xunit{xtarg}));
                
                id_off{xtarg}     = g_items * n_subcategory * (xcate{xtarg}-1);
                xrelatives{xtarg} = (1:(g_items * n_subcategory)) + id_off{xtarg};
                xitem_id{xtarg}   = xitem{xtarg} + (g_items * (xsubcate{xtarg}-1)) + (g_items * n_subcategory * (xcate{xtarg}-1));
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%
            %*************** N vs. N+1 category
            if xcate{1} == xcate{3}, xsame = 1;
            else                     xsame = 2; end
            
            if (xcond==2) || (xcond==3)
                %*************** from N trial: new category
                xtarg = 2;% N-new category
                xcate{xtarg}      = unique(xmatrix(findCol(xheader, {'new_category'}), xunit{xtarg}));
                xsubcate{xtarg}   = unique(xmatrix(findCol(xheader, {'new_subcategory'}), xunit{xtarg}));
                id_off{xtarg}     = g_items * n_subcategory * (xcate{xtarg}-1);
                xrelatives{xtarg} = (1:(g_items * n_subcategory)) + id_off{xtarg};
                
                %*************** N-new(2) vs. N+1(3) category
                if xcate{2} == xcate{3}, xnewsame = 1;
                else                     xnewsame = 2; end
                
                if (xcond==2)
                    if (xsame == 1) %(same:N == N+1 category: N-new ~= N+1)
                        % 1_N-target(1,fa), 2_N-newtarg(18,fr), 3_N+1-target(1,fa),
                        % 4_Ns-relatives(16,fa)
                        % 6_others(18,sc)
                        xitem_id{2} = xrelatives{2};
                        xitem_id{4} = xrelatives{1}(~ismember(xrelatives{1}, [xitem_id{1} xitem_id{3}]));
                        xitem_id{6} = item_array(~ismember(item_array, [xrelatives{1} xrelatives{2}]));
                    else %(diff:N ~= N+1 category)
                        if xnewsame == 1 %(same:N-new == N+1)
                            % 1_N-target(1,fa), 2_N-newtarg(17,fr), 3_N+1-target(1,fr)
                            % 4_N-relative(17,fa)
                            % 6_others(18,sc)
                            xitem_id{2} = xrelatives{2}(~ismember(xrelatives{2}, xitem_id{3}));
                            xitem_id{4} = xrelatives{1}(~ismember(xrelatives{1}, xitem_id{1}));
                            xitem_id{6} = item_array(~ismember(item_array, [xrelatives{1} xrelatives{2}]));
                        else %(diff:N-new ~= N+1) ********
                            % 1_N-target(1,fa), 2_N-newtarg(18,fr), 3_N+1-target(1,sc),
                            % 4_N-relative(17,fa), 5_N+1-relative(17,sc)
                            xitem_id{2} = xrelatives{2};
                            xitem_id{4} = xrelatives{1}(~ismember(xrelatives{1}, xitem_id{1}));
                            xitem_id{5} = xrelatives{3}(~ismember(xrelatives{3}, xitem_id{3}));
                        end
                    end
                else
                    if (xsame == 1) %(same:N == N+1 category: N-new == N+1)
                        % 1_N-target(1,fa), 2_N-newtarg(16,fa), 3_N+1-target(1,fa),
                        % 6_others(36,fr sc)
                        xitem_id{2} = xrelatives{1}(~ismember(xrelatives{1}, [xitem_id{1} xitem_id{3}]));
                        xitem_id{6} = item_array(~ismember(item_array, xrelatives{1}));
                    else %(diff:N ~= N+1 category: N-new ~= N+1)
                        % 1_N-target(1,fa), 2_N-newtarg(17,fa), 3_N+1-target(1,fr),
                        % 5_N+1-relative(17,fr)
                        % 6_others(18,sc)
                        xitem_id{2} = xrelatives{1}(~ismember(xrelatives{1}, xitem_id{1}));
                        xitem_id{5} = xrelatives{3}(~ismember(xrelatives{3}, xitem_id{3}));
                        xitem_id{6} = item_array(~ismember(item_array, [xrelatives{1} xrelatives{3}]));
                    end
                end
            else % Ma,Su,Cl
                if (xsame == 1) %(same:N == N+1 category)
                    % 1_N-target(1,fa), 3_N+1-target(1,fa)
                    % 4_Ns-relatives(16,fa),
                    % 6_others(36,fr sc)
                    xitem_id{4} = xrelatives{1}(~ismember(xrelatives{1}, [xitem_id{1} xitem_id{3}]));
                    xitem_id{6} = item_array(~ismember(item_array, xrelatives{1}));
                else %(diff:N ~= N+1 category) ********
                    % 1_N-target(1,fa), 3_N+1-target(1,fr)
                    % 4_N-relative(17,fa), 5_N+1-relative(17,fr),
                    % 6_others(18,sc)
                    xitem_id{4} = xrelatives{1}(~ismember(xrelatives{1}, xitem_id{1}));
                    xitem_id{5} = xrelatives{3}(~ismember(xrelatives{3}, xitem_id{3}));
                    xitem_id{6} = item_array(~ismember(item_array, [xrelatives{1} xrelatives{3}]));
                end
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %*************** check item ids
            t_array = [];
            for xtarg = 1:n_targs
                if ~isempty(xitem_id{xtarg})
                    t_array = horzcat(t_array, xitem_id{xtarg});
                end
            end
            t_array = sort(unique(t_array));
            
            if ~isequal(item_array, t_array)
                it_run = unique(xmatrix(findCol(xheader, {'it_run'}), xunit{1}));
                fprintf('warning: no full item_id: run %d, trial %d, cond %d\n', it_run, xtrial, xcond)
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %*************** template:
            % 1_N, 2_N-new, 3_N+1, 4_N-relative, 5_N+1-relative, 6_others
            for xtarg = 1:n_targs
                if ~isempty(xitem_id{xtarg})
                    for i = 1:length(xitem_id{xtarg})
                        xit = xitem_id{xtarg}(i);
                        xitem_pat_w{xtarg}{i} = subj_rsa{1}.template.items{xit}.pat_w;
                    end
                end
            end
            
            %*************** N+1 timecourse
            xonset   = find(xunit{3}, 1);
            xunit_tc = xonset:(xonset + n_trs - 1);
            xspike   = ~(xmatrix(findCol(xheader, {'spike'}), xunit_tc));
            
            if sum(xspike)~=0
                for xtr = 1:length(xunit_tc)
                    if xspike(xtr)
                        for xtarg = 1:n_targs
                            if ~isempty(xitem_id{xtarg})
                                t_rep_corr = [];
                                
                                for i = 1:length(xitem_id{xtarg})
                                    %%%%%%%%%%%%%%%%%%%%%%%%%% N trial
                                    clear ytarg_pat ytarg_pat_w
                                    %*************** target
                                    xit         = xitem_id{xtarg}(i);
                                    ycate       = fix((xit-1)/(g_items * n_subcategory)) + 1;
                                    xbeta       = xsubj{1}{ycate}.patterns{xit + 1}.mat;
                                    ytarg_pat   = xsubj{2}{ycate}.patterns{1}.mat(:, xunit_tc(xtr));
                                    ytarg_pat_w = (ytarg_pat' .* xbeta')';
                                    
                                    %*************** similarity: z-transformed pearson correlation
                                    t_rep_corr = horzcat(t_rep_corr, corr2(xitem_pat_w{xtarg}{i}, ytarg_pat_w));
                                end
                                
                                xwr_rep = mean(t_rep_corr);
                                
                                if (xcond==2) && (xsame == 2)
                                    react.cond{xcond}.same{xsame}.newsame{xnewsame}.targ{xtarg}.tr{xtr}.corr_w = ...
                                        horzcat(react.cond{xcond}.same{xsame}.newsame{xnewsame}.targ{xtarg}.tr{xtr}.corr_w, ...
                                        .5*(log(1+xwr_rep) - log(1-xwr_rep)));
                                else
                                    react.cond{xcond}.same{xsame}.targ{xtarg}.tr{xtr}.corr_w = ...
                                        horzcat(react.cond{xcond}.same{xsame}.targ{xtarg}.tr{xtr}.corr_w, ...
                                        .5*(log(1+xwr_rep) - log(1-xwr_rep)));
                                end
                            end
                        end
                    end
                end
            end
        end
        fprintf('\n')
    end
    
    fprintf('\n')
    
    %% ============= ADDITIONAL
    %*************** STUDY: 1st presentation of target decoding
    % template from repetition / per run
    % subj_rsa{1}.template.items{1}.repeat{xrep}
    % subj_rsa{1}.template.items{1}.run{loc_run}
    xpercept_win = args.percept_win;
    item_ids = unique(xmatrix(findCol(xheader, {'image_id'}), :));
    item_ids = item_ids(~ismember(item_ids, 0));
    
    for xtarg = 1:3
        for yrun = loc_run_array
            for xcate = 1:n_category
                for xsubcate = 1:n_subcategory    
                    decoding.cate{xcate}.subcate{xsubcate}.targ{xtarg}.rep{yrun}.corr_w = [];
                    decoding.cate{xcate}.subcate{xsubcate}.targ{xtarg}.run{yrun}.corr_w = [];
                end
            end
            
            decoding.targ{xtarg}.rep{yrun}.corr_w = [];
            decoding.targ{xtarg}.run{yrun}.corr_w = [];
        end
    end
    
    for s_item = item_ids %select only the first presented item
        %*************** extract params
        clear item_id 
        
        s_run      = min(getDATA(xmatrix', xheader, ...
            {'image_id'}, {s_item}, findCol(xheader, {'run'})));
        s_trial    = min(getDATA(xmatrix', xheader, ...
            {'image_id','run'}, {s_item, s_run}, findCol(xheader, {'trial'})));
        
        tunit      = getDATA(xmatrix', xheader, ...
            {'image_id','run','trial'}, {s_item, s_run, s_trial});
        xcate      = unique(xmatrix(findCol(xheader, {'category'}), tunit));
        xsubcate   = unique(xmatrix(findCol(xheader, {'subcategory'}), tunit));
        xitem      = unique(xmatrix(findCol(xheader, {'item'}), tunit));
        related_items = (1:n_item * n_subcategory) + (n_item * n_subcategory * (xcate - 1));
        
        %*************** study: 1_target, 2_nontarget, 3_baseline
        item_id{1} = xitem + (g_items * (xsubcate-1)) + (g_items * n_subcategory * (xcate-1));
        item_id{2} = related_items(~ismember(related_items, item_id{1}));
        item_id{3} = item_array(~ismember(item_array, related_items));
        
        %*************** template
        clear xitem_pat_w
        for xtarg = 1:3%1_target, 2_nontarget, 3_baseline
            for i = 1:length(item_id{xtarg})
                xit = item_id{xtarg}(i);
                for yrun = loc_run_array
                    xx = subj_rsa{1}.template.items{xit}.rep{yrun};
                    if ~isempty(xx)
                        xitem_pat_w.targ{xtarg}.item{i}.rep{yrun} = xx.pat_w;
                    else
                        xitem_pat_w.targ{xtarg}.item{i}.rep{yrun} = [];
                    end
                    
                    xx = subj_rsa{1}.template.items{xit}.run{yrun};
                    if ~isempty(xx)
                        xitem_pat_w.targ{xtarg}.item{i}.run{yrun} = xx.pat_w;
                    else
                        xitem_pat_w.targ{xtarg}.item{i}.run{yrun} = [];
                    end
                end
            end
        end
        
        %*************** timecourse
        xunit     = find(getDATA(xmatrix', xheader, {'run', 'trial'}, {s_run, s_trial}));
        xunit_tc  = xunit(xpercept_win);
        xspike    = ~(xmatrix(findCol(xheader, {'spike'}), xunit_tc));
        
        if sum(xspike)~=0
            for xtarg = 1:3%1_target, 2_nontarget, 3_baseline
                for yrun = loc_run_array
                    
                    t_rep_corr = []; t_run_corr = [];
                    
                    for xtr = 1:length(xunit_tc)
                        if xspike(xtr)
                            for i = 1:length(item_id{xtarg})
                                clear ytarg_pat ytarg_pat_w
                                %*************** target
                                xit         = item_id{xtarg}(i);
                                ycate       = fix((xit-1)/(g_items * n_subcategory)) + 1;
                                
                                if (xtarg==1) && (xcate~=ycate)
                                    fprintf('warning: not matching category for target item\n')
                                end
                                
                                xbeta       = xsubj{1}{ycate}.patterns{xit + 1}.mat;
                                ytarg_pat   = xsubj{2}{ycate}.patterns{1}.mat(:, xunit_tc(xtr));
                                ytarg_pat_w = (ytarg_pat' .* xbeta')';
                                
                                %*************** similarity: z-transformed pearson correlation
                                clear xtemp_pat_w
                                xtemp_pat_w = xitem_pat_w.targ{xtarg}.item{i}.rep{yrun};
                                
                                if ~isempty(xtemp_pat_w)
                                    xr         = corr2(xtemp_pat_w, ytarg_pat_w);
                                    t_rep_corr = horzcat(t_rep_corr, .5*(log(1+xr) - log(1-xr)));
                                end
                                
                                clear xtemp_pat_w
                                xtemp_pat_w = xitem_pat_w.targ{xtarg}.item{i}.run{yrun};
                                    
                                if ~isempty(xtemp_pat_w)
                                    xr         = corr2(xtemp_pat_w, ytarg_pat_w);
                                    t_run_corr = horzcat(t_run_corr, .5*(log(1+xr) - log(1-xr)));
                                end
                            end
                        end
                    end
                    
                    %*************** fisher's z-correlation
                    clear xwr_rep xwr_run
                    if ~isempty(t_rep_corr)
                        decoding.cate{xcate}.subcate{xsubcate}.targ{xtarg}.rep{yrun}.corr_w = ...
                            horzcat(decoding.cate{xcate}.subcate{xsubcate}.targ{xtarg}.rep{yrun}.corr_w,...
                             mean(t_rep_corr));
                        
                        decoding.targ{xtarg}.rep{yrun}.corr_w = ...
                            horzcat(decoding.targ{xtarg}.rep{yrun}.corr_w,...
                            mean(t_rep_corr));
                    end
                    
                    if ~isempty(t_run_corr)
                        decoding.cate{xcate}.subcate{xsubcate}.targ{xtarg}.run{yrun}.corr_w = ...
                            horzcat(decoding.cate{xcate}.subcate{xsubcate}.targ{xtarg}.run{yrun}.corr_w,...
                            mean(t_run_corr));
                        
                        decoding.targ{xtarg}.run{yrun}.corr_w = ...
                            horzcat(decoding.targ{xtarg}.run{yrun}.corr_w,...
                            mean(t_run_corr));
                    end
                end
            end
        end
    end
    
% else
%     %*************** load existing file
%     xrsa = load(rsa_fname);%,'subj_rsa'
%     
%     subj_rsa{xph}.rsa        = xrsa.subj_rsa{xph}.rsa;
%     subj_rsa{xph}.timecourse = xrsa.subj_rsa{xph}.timecourse;
%     subj_rsa{xph}.react      = xrsa.subj_rsa{xph}.react;
%     subj_rsa{xph}.decoding   = decoding;

end

%% ============= SAVE PATTERNS
%*************** reset subj_rsa
if args.rsa_pattern(2)
    xph = 2;
    subj_rsa{xph}.decoding   = decoding;
    subj_rsa{xph}.rsa        = rsa;
    subj_rsa{xph}.timecourse = timecourse;
    subj_rsa{xph}.react      = react;
end

save(rsa_fname,'subj_rsa','-v7.3');
fsize = dir(rsa_fname);
fprintf('saved: file size: %s GB\n', num2str(fsize.bytes/(10^9)));

toc

%% *************** copy pattern to grp_pattern directory
% dst_fname = fullfile(dirs.rsa.group.pattern, sprintf('%s_%s.mat', basename, args.subject_id)); 
% copyfile(rsa_fname, dst_fname);

