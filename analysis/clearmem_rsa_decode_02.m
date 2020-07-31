function[] = clearmem_rsa_decode_02(args, dirs)
%*************** create template by averaging the extracted pattern
%*************** category level

%% ============= SETUP DIRECTORY
spm_dir         = dirs.rsa.roi.spm;
output_dir      = dirs.rsa.pattern;

%% ============= SETUP PARAMETERS
n_category      = args.index{1}.param.n_category;
category_names  = {'face','fruit','scene'};
conds_names     = args.index{2}.param.conds_names;
n_condition     = length(conds_names);
n_trials        = args.index{2}.param.n_trials;%of study
cate_members    = 1:n_category;
n_trs           = args.tc_tr;

spm_p_thresh    = args.spm_p_thresh;
spm_v_extent    = args.spm_v_extent;
spm_thresh_type = args.spm_thresh_type;

rsa_mask        = args.rsa_mask;

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
tic

%% ============= 01: LOCALIZER EPI PATTERNS
%*************** load mask + read in epis
%*************** zsoring epis
%*************** option: read epis in mask | wholebrain
%*************** define: args.wholebrain '0 | 1'
xph     = 1; args.xphase = xph;
xmatrix = args.regs{xph}.matrix;
xheader = args.regs{xph}.header;
n_runs  = args.index{xph}.param.n_runs;

for xmaskcate = 1:n_category
    %% ============= INITIATE SUBJ. READ MASK FROM TSPM
    xmask_name = sprintf('spmT_000%d_%s_thresh%s_ext%s_%s_mask',...
        xmaskcate, spm_thresh_type, num2str(spm_p_thresh), num2str(spm_v_extent), rsa_mask);
        
    ph1.basename = sprintf('%s_%s_%s_zscored_%s', args.phase_name{xph}, category_names{xmaskcate}, xmask_name, args.epi_name);
    fname        = sprintf('%s/ph1_%s.mat', dirs.rsa.scratch, ph1.basename);
        
    if args.load_pattern_rsa{xph}
        clear subj nm
        fprintf('\n(+) initiating subject structure: %s, %s \n', category_names{xmaskcate}, args.subject_id);
        subj = init_subj(args.experiment, args.subject_id);%identifier of the subj
        
        fprintf('\n#####################################################################\n');
        fprintf('* STEP_01: %s: load mask + read/zscore epis patterns\n', args.phase_name{xph});
        
        %*************** MASK
        nm.mask = xmask_name;
        fprintf('\n(+) load mask: %s.%s\n', nm.mask, args.epiext);
        
        xmask    = fullfile(spm_dir, sprintf('%s.%s', nm.mask, args.epiext));
        subj     = load_spm_mask(subj, nm.mask, xmask);
        mask_obj = get_object(subj, 'mask', nm.mask); %'mask object' contents
        subj.masks{1}.header.description = ...
            sprintf('% mask from TSPM', category_names{xmaskcate});
        
        %% ============= EXTRACT PATTERNS: LOCALIZER
        %*************** EPI PATTERN: subj.patterns{1}
        fprintf('\n(+) load epi data under mask %s with %s voxels\n', ...
            nm.mask, num2str(mask_obj.nvox));
        
        nm.pats = sprintf('%s_%s_patterns', args.phase_name{xph}, category_names{xmaskcate});
        
        for xrun = 1:n_runs
            args.epi_names{xph}{xrun} = ...
                fullfile(dirs.runs{xph}{xrun}, sprintf('%s.%s.gz', args.epi_name, args.epiext));
        end
        
        fprintf('\n... loading epi patterns\n');
        subj = load_spm_pattern_gz(subj, nm.pats, nm.mask, args.epi_names{xph});
        
        %% ============= RESET PATTERNS FOR TRIM-TRS
        %*************** trim the first 10TRs of each run if untrimed
        if (~args.trim_trs)% 0_non_trimmed, 1_trimmed
            
            %*************** reset patterns
            fprintf('\n  ... reset patterns: trim the first %sTRs: from subj.patterns{1}.mat\n', ...
                num2str(args.xtrim));
            
            xmatrix = args.regs{xph}.matrix;%3180
            xheader = args.regs{xph}.header;
            xunit   = xmatrix(findCol(xheader, {'it_volume'}), :);
            
            selected_pattern = subj.patterns{1}.mat(:, xunit);
            
            subj    = set_mat(subj, 'pattern', subj.patterns{1}.name, selected_pattern);
            summarize(subj, 'objtype', 'pattern');
        end
        
        %% ============= ZSCORING: subj.patterns{2}
        %*************** the standard deviation of the timecourse is one.
        fprintf('\n(+) z-scoring data: voxel-wise\n');
        
        %*************** SELECTORS: subj.selectors{1}
        nm.runs   = sprintf('%s_runs', args.phase_name{xph});% selectors object name
        subj      = initset_object(subj,'selector', nm.runs, args.regs{xph}.selectors);
        subj      = zscore_runs(subj, nm.pats, nm.runs);
        
        nm.pats_z = sprintf('%s_z', nm.pats);% subj.patterns{2}.name;
        
        summarize(subj,'objtype','pattern')
        
        %% ============= EXTRACTING WEIGHT(BETA) FROM SPM
        subj = load_spm_pattern(subj, 'spm_beta', nm.mask, ...
            fullfile(dirs.rsa.roi.spm, sprintf('spmT_000%d.%s', xmaskcate+1, args.epiext)));
        subj.patterns{3}.header.description = ...
            sprintf('% beta from TSPM', category_names{xmaskcate});
        
        %*************** save phase 1.
        ph1.args = args; ph1.nm = nm; ph1.subj = subj;
        save(fname,'ph1','-v7.3');
        
    else
        fprintf('  ... load patterns from ph1\n'); load(fname); 
        subj = ph1.subj; nm = ph1.nm;
    end
    
    summarize(subj)
    
    %% ============= EXTRACTING WEIGHT(BETA) FROM SPM
    xbeta = subj.patterns{3}.mat;
    patterns{xph}.maskcate{xmaskcate}.beta    = xbeta; %#ok<*AGROW>
    patterns{xph}.maskcate{xmaskcate}.allPats = subj.patterns{2}.mat;
    
    %% ============= EXTRACTING PATTERNS FROM LOCALIZER
    %*************** TR shifted
    for xcate = 1:n_category
        for xrun = 1:n_runs
            clear xpat xpat_w
            
            fprintf('localizer: mask cate: %s, targ: %s, run: %d\n', ...
                category_names{xmaskcate}, category_names{xcate}, xrun); 
            
            %*************** localizer zscored raw patterns
            xunit    = find(getDATA(xmatrix', xheader, {'run','category','spike'}, {xrun, xcate, 0}));
            xunit_sh = xunit + args.shift_TRs;
            
            if xunit_sh(end) > size(subj.patterns{2}.mat, 2)
                endTR    = size(subj.patterns{2}.mat, 2);
                endunit  = find(xunit_sh == endTR);
                xunit_sh = xunit_sh(1:endunit);
            end
            
            if ~isempty(xunit_sh)
                fprintf('xunit_sh_end: %s\n',  num2str(xunit_sh(end)));
                
                xpat = subj.patterns{2}.mat(:, xunit_sh);
                
                patterns{xph}.maskcate{xmaskcate}.targcate{xcate}.meanPat(:,xrun) = mean(xpat,2);
                
                fprintf('xpat size: %s\n', num2str(size(xpat,2)));
                
                %*************** weighted with beta
                for i = 1:size(xpat,2)
                    xxpat       = xpat(:,i)';
                    xpat_w(:,i) = xxpat.*xbeta';
                end                
               
                patterns{xph}.maskcate{xmaskcate}.targcate{xcate}.weighted.meanPat(:,xrun) = mean(xpat_w,2);
            end
        end
        
        patterns{xph}.maskcate{xmaskcate}.targcate{xcate}.template = ...
            mean(patterns{xph}.maskcate{xmaskcate}.targcate{xcate}.meanPat, 2);
        patterns{xph}.maskcate{xmaskcate}.targcate{xcate}.weighted.template = ...
            mean(patterns{xph}.maskcate{xmaskcate}.targcate{xcate}.weighted.meanPat, 2);
    end
    
    patterns{xph}.maskcate{xmaskcate}.nvoxels  = size(subj.patterns{2}.mat, 1);
    t_vox(xmaskcate) = patterns{xph}.maskcate{xmaskcate}.nvoxels;
    
end

patterns{xph}.info.nVox   = t_vox;
patterns{xph}.info.nTRs   = size(subj.patterns{2}.mat, 2);
patterns{xph}.info.matrix = xmatrix;
patterns{xph}.info.header = xheader;

%% ============= 02: STUDY EPI PATTERNS
%*************** load mask + read in epis
%*************** zsoring epis
%*************** option: read epis in mask | wholebrain
%*************** define: args.wholebrain '0 | 1'
xph       = 2; args.xphase = xph;

xparam    = args.index{xph}.param;
dur_stim  = xparam.dur_stim;
dur_manip = xparam.dur_manipulation;
dur_iti   = max(xparam.dur_iti);
dur_sync  = n_trs - dur_stim - dur_manip - dur_iti;%next trial 

xmatrix   = args.regs{xph}.matrix;
xheader   = args.regs{xph}.header;
n_runs    = args.index{xph}.param.n_runs;

for xmaskcate = 1:n_category
    %% ============= INITIATE SUBJ. READ MASK FROM TSPM
    xmask_name = sprintf('spmT_000%d_%s_thresh%s_ext%s_%s_mask',...
        xmaskcate, spm_thresh_type, num2str(spm_p_thresh), num2str(spm_v_extent), rsa_mask);
        
    ph2.basename = sprintf('%s_%s_%s_zscored_%s', args.phase_name{xph}, category_names{xmaskcate}, xmask_name, args.epi_name);
    fname        = sprintf('%s/ph2_%s.mat', dirs.rsa.scratch, ph2.basename);
        
    if args.load_pattern_rsa{xph}
        clear subj nm
        fprintf('\n(+) initiating subject structure: %s, %s \n', category_names{xmaskcate}, args.subject_id);
        subj = init_subj(args.experiment, args.subject_id);%identifier of the subj
        
        fprintf('\n#####################################################################\n');
        fprintf('* STEP_02: %s: load mask + read/zscore epis patterns\n', args.phase_name{xph});
        
        %*************** MASK
        nm.mask = xmask_name;
        fprintf('\n(+) load mask: %s.%s\n', nm.mask, args.epiext);
        
        xmask    = fullfile(spm_dir, sprintf('%s.%s', nm.mask, args.epiext));
        subj     = load_spm_mask(subj, nm.mask, xmask);
        mask_obj = get_object(subj, 'mask', nm.mask); %'mask object' contents
        subj.masks{1}.header.description = ...
            sprintf('% mask from TSPM', category_names{xmaskcate});
        
        %% ============= EXTRACT PATTERNS: LOCALIZER
        %*************** EPI PATTERN: subj.patterns{1}
        fprintf('\n(+) load epi data under mask %s with %s voxels\n', ...
            nm.mask, num2str(mask_obj.nvox));
        
        nm.pats = sprintf('%s_%s_patterns', args.phase_name{xph}, category_names{xmaskcate});
        
        for xrun = 1:n_runs
            args.epi_names{xph}{xrun} = ...
                fullfile(dirs.runs{xph}{xrun}, sprintf('%s.%s.gz', args.epi_name, args.epiext));
        end
        
        fprintf('\n... loading epi patterns\n');
        subj = load_spm_pattern_gz(subj, nm.pats, nm.mask, args.epi_names{xph});
        
        %% ============= RESET PATTERNS FOR TRIM-TRS
        %*************** trim the first 10TRs of each run if untrimed
        if (~args.trim_trs)% 0_non_trimmed, 1_trimmed
            
            %*************** reset patterns
            fprintf('\n  ... reset patterns: trim the first %sTRs: from subj.patterns{1}.mat\n', ...
                num2str(args.xtrim));
            
            xmatrix = args.regs{xph}.matrix;%3180
            xheader = args.regs{xph}.header;
            xunit   = xmatrix(findCol(xheader, {'it_volume'}), :);
            
            selected_pattern = subj.patterns{1}.mat(:, xunit);
            
            subj    = set_mat(subj, 'pattern', subj.patterns{1}.name, selected_pattern);
            summarize(subj, 'objtype', 'pattern');
        end
        
        %% ============= ZSCORING: subj.patterns{2}
        %*************** the standard deviation of the timecourse is one.
        fprintf('\n(+) z-scoring data: voxel-wise\n');
        
        %*************** SELECTORS: subj.selectors{1}
        nm.runs   = sprintf('%s_runs', args.phase_name{xph});% selectors object name
        subj      = initset_object(subj,'selector', nm.runs, args.regs{xph}.selectors);
        subj      = zscore_runs(subj, nm.pats, nm.runs);
        
        nm.pats_z = sprintf('%s_z', nm.pats);% subj.patterns{2}.name;
        
        summarize(subj,'objtype','pattern')
        
        %*************** save phase 1.
        ph2.args = args; ph2.nm = nm; ph2.subj = subj;
        save(fname,'ph2','-v7.3');
        
    else
        fprintf('  ... load patterns from ph2\n'); load(fname); 
        subj = ph2.subj; nm = ph2.nm;
    end
    
    summarize(subj)
    
    %% ============= EXTRACTING PATTERNS FROM STUDY
    %*************** study zscored raw patterns: all TRS
    %*************** no TR shift for timecourse RSA
    %*************** TR shift for accuracy
    patterns{xph}.maskcate{xmaskcate}.allPats = subj.patterns{2}.mat;
    
    %*************** category based: verify accuracy
    for xcate = 1:n_category
        for xrun = 1:n_runs
            clear xpat 
            
            fprintf('study mask cate: %s, targ: %s, run: %d\n', ...
                category_names{xmaskcate}, category_names{xcate}, xrun); 
            
            %*************** localizer zscored raw patterns
            xunit    = find(getDATA(xmatrix', xheader, ...
                {'run','category','presentation','spike'}, {xrun, xcate, 1, 0}));
            xunit_sh = xunit + args.shift_TRs;
            xpat     = subj.patterns{2}.mat(:, xunit_sh);

            patterns{xph}.maskcate{xmaskcate}.targcate{xcate}.meanPat(:,xrun) = mean(xpat,2);
        end
        
        patterns{xph}.maskcate{xmaskcate}.targcate{xcate}.template = ...
            mean(patterns{xph}.maskcate{xmaskcate}.targcate{xcate}.meanPat, 2);
    end
    
    %% *************** info
    patterns{xph}.maskcate{xmaskcate}.nvoxels = size(subj.patterns{2}.mat, 1);
    t_vox(xmaskcate) = patterns{xph}.maskcate{xmaskcate}.nvoxels;
    
end

patterns{xph}.info.nVox   = t_vox;
patterns{xph}.info.nTRs   = size(subj.patterns{2}.mat, 2);
patterns{xph}.info.matrix = xmatrix;
patterns{xph}.info.header = xheader; %#ok<*NASGU>

%% *************** STUDY: timecourse: condition
% patterns{xph}.maskcate{xmaskcate}.allPats
% patterns{xph}.timecourse.cond{xcond}.target{xtarg}.cate{xcate}.trPats{xtr} 
% xtarg: 1_target, 2_newtarg/nontarg, 3_nontarg

for xcond = 1:n_condition
    for xtarg = 1:3
        for xcate = 1:n_category
            for xtr = 1:n_trs
                if xcond~=2
                    patterns{xph}.timecourse.cond{xcond}.targcate{xcate}.targ{xtarg}.trPats{xtr} = [];
                else
                    for xnewtarg = 1:n_category
                        patterns{xph}.timecourse.cond{xcond}.targcate{xcate}.newtarg{xnewtarg}.targ{xtarg}.trPats{xtr} = [];                        
                    end
                end
            end
            
            for xtr = 1:dur_sync
                if xcond~=2
                    patterns{xph}.timecourse.sync.cond{xcond}.targcate{xcate}.targ{xtarg}.trPats{xtr} = [];
                else
                    for xnewtarg = 1:n_category
                        patterns{xph}.timecourse.sync.cond{xcond}.targcate{xcate}.newtarg{xnewtarg}.targ{xtarg}.trPats{xtr} = [];                        
                    end
                end
            end
        end
    end
end

%% *************** timecourse
% xtarg: 1_target, 2_newtarg/nontarg, 3_nontarg

xph = 2;
for xrun = 1:n_runs
    for xtrial = 1:n_trials(xrun)
        xunit = find(getDATA(xmatrix', xheader, {'run', 'trial'}, {xrun, xtrial}));
        xcond = unique(xmatrix(findCol(xheader, {'condition'}), xunit));
        xcate = unique(getDATA(xmatrix', xheader, ...
            {'run', 'trial','presentation'}, {xrun, xtrial, 1}, ...
            findCol(xheader, {'category'})));
        
        if xcond == 2
            %*************** new category for cond 2/3
            xnewcate = unique(getDATA(xmatrix', xheader, ...
                {'run', 'trial','presentation'}, {xrun, xtrial, 2}, ...
                findCol(xheader, {'new_category'})));
            
            xnoncate = cate_members(~ismember(cate_members, [xcate xnewcate]));
            
        else
            %*************** non target category
            xnoncate = cate_members(~ismember(cate_members, xcate));
        end
        
        xunit_tc = xunit(1):(xunit(1) + n_trs - 1);
        xspike   = ~(xmatrix(findCol(xheader, {'spike'}), xunit_tc));
        
        if sum(xspike)~=0
            for xtr = 1:length(xunit_tc)
                if xspike(xtr)
                    clear xpat
                    %*************** target
                    xpat{1} = patterns{xph}.maskcate{xcate}.allPats(:, xunit_tc(xtr));
                        
                    if xcond~=2
                        %*************** non target
                        for xx = 1:length(xnoncate)
                            xpat{xx+1} = patterns{xph}.maskcate{xnoncate(xx)}.allPats(:, xunit_tc(xtr)); 
                        end
                        
                        for xtarg = 1:3
                            patterns{xph}.timecourse.cond{xcond}.targcate{xcate}.targ{xtarg}.trPats{xtr} = ...
                                horzcat(patterns{xph}.timecourse.cond{xcond}.targcate{xcate}.targ{xtarg}.trPats{xtr},...
                                xpat{xtarg});
                            
                            if xtr > length(xunit)
                                it_tr = xtr - length(xunit);
                                if it_tr <= dur_sync
                                    patterns{xph}.timecourse.sync.cond{xcond}.targcate{xcate}.targ{xtarg}.trPats{it_tr} = ...
                                        horzcat(patterns{xph}.timecourse.sync.cond{xcond}.targcate{xcate}.targ{xtarg}.trPats{it_tr},...
                                        xpat{xtarg});
                                end
                            end
                        end
                    else
                        %*************** new category
                        xpat{2} = patterns{xph}.maskcate{xnewcate}.allPats(:, xunit_tc(xtr)); 
                        %*************** non target
                        xpat{3} = patterns{xph}.maskcate{xnoncate}.allPats(:, xunit_tc(xtr)); 
                        
                        for xtarg = 1:3
                            patterns{xph}.timecourse.cond{xcond}.targcate{xcate}.newtarg{xnewcate}.targ{xtarg}.trPats{xtr} = ...
                                horzcat(patterns{xph}.timecourse.cond{xcond}.targcate{xcate}.newtarg{xnewcate}.targ{xtarg}.trPats{xtr},...
                                xpat{xtarg});
                            
                            if xtr > length(xunit)
                                it_tr = xtr - length(xunit);
                                if it_tr <= dur_sync
                                    patterns{xph}.timecourse.sync.cond{xcond}.targcate{xcate}.newtarg{xnewcate}.targ{xtarg}.trPats{it_tr} = ...
                                        horzcat(patterns{xph}.timecourse.sync.cond{xcond}.targcate{xcate}.newtarg{xnewcate}.targ{xtarg}.trPats{it_tr},...
                                        xpat{xtarg});
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

%% *************** condition: timecourse: mean per TR
for xcond = 1:n_condition
    for xcate = 1:n_category
        for xtarg = 1:3
            for xtr = 1:n_trs
                if xcond~=2
                    patterns{xph}.timecourse.cond{xcond}.targcate{xcate}.targ{xtarg}.mean_trPats{xtr} = ...
                        mean(patterns{xph}.timecourse.cond{xcond}.targcate{xcate}.targ{xtarg}.trPats{xtr}, 2);
                else
                    for xnewcate = 1:n_category
                        patterns{xph}.timecourse.cond{xcond}.targcate{xcate}.newtarg{xnewcate}.targ{xtarg}.mean_trPats{xtr} = ...
                            mean(patterns{xph}.timecourse.cond{xcond}.targcate{xcate}.newtarg{xnewcate}.targ{xtarg}.trPats{xtr}, 2);
                    end
                end
            end
            
            for xtr = 1:dur_sync
                if xcond~=2
                    patterns{xph}.timecourse.sync.cond{xcond}.targcate{xcate}.targ{xtarg}.mean_trPats{xtr} = ...
                        mean(patterns{xph}.timecourse.sync.cond{xcond}.targcate{xcate}.targ{xtarg}.trPats{xtr}, 2);
                else
                    for xnewcate = 1:n_category
                        patterns{xph}.timecourse.sync.cond{xcond}.targcate{xcate}.newtarg{xnewcate}.targ{xtarg}.mean_trPats{xtr} = ...
                            mean(patterns{xph}.timecourse.sync.cond{xcond}.targcate{xcate}.newtarg{xnewcate}.targ{xtarg}.trPats{xtr}, 2);
                    end
                end
            end
        end
    end
end

%% ============= SAVE PATTERNS
basename = sprintf('patterns_%s_spmT_%s_thresh%s_ext%s_%s_mask', ...
    args.level, spm_thresh_type, num2str(spm_p_thresh), num2str(spm_v_extent), rsa_mask);
fname = fullfile(dirs.rsa.pattern, sprintf('%s.mat', basename)); 

save(fname,'patterns','-v7.3');

fsize = dir(fname);

fprintf('saved: file size: %s GB\n', num2str(fsize.bytes/(10^9)));

toc
        
%% *************** copy pattern to grp_pattern directory
% dst_fname = fullfile(dirs.rsa.group.pattern, sprintf('%s_%s.mat', basename, args.subject_id)); 
% copyfile(fname, dst_fname);
% 
% fsize = dir(fname);
% fprintf('file size: %s GB\n', num2str(fsize.bytes/(10^9)));


