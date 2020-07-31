function[subj, nm] = mvpa_ph01_patterns_gz(args, subj, nm, dirs)
    
    %---------------------------------------------------------------------------
    %*************** load mask + read in epis.
    %*************** zsoring epis
    %*************** options: read epis in mask | wholebrain
    %---------------------------------------------------------------------------    
    
    xph = args.xphase;
    
    %% ============= IN MASK
    if args.xphase == 2
        mask_dir  = dirs.mvpa.selected_mask{xph}; 
    else
        mask_dir = dirs.mask; 
%         if args.wholebrain, mask_dir = dirs.epi_mid;
%         else                mask_dir = dirs.mask; end
    end
    
    %*************** MASK
    nm.mask{xph} = args.mask_name;
    fprintf('\n(+) load mask: %s.%s\n', nm.mask{xph}, args.epiext);
    
    xmask_gz = fullfile(mask_dir, sprintf('%s.%s.gz', args.mask_name, args.epiext));

    subj     = load_spm_mask_gz(subj, nm.mask{xph}, xmask_gz);
    mask_obj = get_object(subj, 'mask', nm.mask{xph}); %'mask object' contents
    subj.masks{1}.header.description = 'raw pattern under the mask';
    
    %% ============= EPI PATTERN: subj.patterns{1}
    fprintf('\n(+) load epi data under mask %s with %s voxels\n', ...
        args.mask_name, num2str(mask_obj.nvox));
    
    nm.pats{xph} = sprintf('%s_patterns', args.phase_name{xph});

    %*************** shift epi
    if strcmp(args.regress_type, 'shift')
        for xrun = 1:length(dirs.runs{xph})
            args.epi_names{xph}{xrun} = ...
                fullfile(dirs.runs{xph}{xrun}, sprintf('%s.%s.gz', args.epi_name, args.epiext));
        end
                
        fprintf('\n... loading epi patterns: n_runs: %d\n', length(dirs.runs{xph}));
        
        %*************** read pattern from operation
        if args.read_mat{xph}
            
            xwhole_mask = find(args.whole_pattern.masks{1}.mat);
            xmask       = find(subj.masks{2}.mat);
            new_mat     = ismember(xwhole_mask, xmask);
            xpattern    = args.whole_pattern.patterns{1}.mat(new_mat, :);
            
            subj = initset_object(subj,'pattern', nm.pats{xph}, xpattern, 'masked_by', nm.mask{2});
            
            subj.patterns{4}.header = args.whole_pattern.patterns{1}.header;
            
        else
            subj = load_spm_pattern_gz(subj, nm.pats{xph}, nm.mask{xph}, args.epi_names{xph});
        end
         
    elseif strcmp(args.regress_type, 'beta')
        args.epi_names = fullfile(dirs.beta.home, ...
                         sprintf('LSS_beta_localizer_%s.nii.gz', args.subject_id));
        
        fprintf('... loading beta patterns\n');
        subj = load_spm_pattern_gz(subj, nm.pats{xph}, nm.mask{xph}, args.epi_names);
    end
    
    %% ============= RESET PATTERNS
    %*************** trim the first 10TRs of each run if untrimed 
    if (~args.trim_trs) && (~args.read_mat{xph})% 0_non_trimmed, 1_trimmed
        
        if xph == 2, xcell = 4; else xcell = 1; end
        
        %*************** reset patterns
        fprintf('\n  ... reset patterns: trim the first %sTRs: from subj.patterns{%d}.mat\n', ...
            num2str(args.xtrim), xcell);

        xmatrix = args.regs{xph}.matrix;%6990
        xheader = args.regs{xph}.header;
        xunit   = xmatrix(findCol(xheader, {'it_volume'}), :);
        
        selected_pattern = subj.patterns{xcell}.mat(:, xunit);
        
        subj    = set_mat(subj, 'pattern', subj.patterns{xcell}.name, selected_pattern);
        summarize(subj, 'objtype', 'pattern');
    end
    
    %% ============= ZSCORING: subj.patterns{2}
    %*************** zsoring each voxel, include 'rest' in the zscoring
    %*************** the standard deviation of the timecourse is one.
    fprintf('\n(+) z-scoring data: voxel-wise\n');
    
    %*************** SELECTORS: subj.selectors{1}
    nm.runs{xph}   = sprintf('%s_runs', args.phase_name{xph});% selectors object name
    subj           = initset_object(subj,'selector', nm.runs{xph}, args.regs{xph}.selectors);
    subj           = zscore_runs(subj, nm.pats{xph}, nm.runs{xph});
    
    nm.pats_z{xph} = sprintf('%s_z', nm.pats{xph});% subj.patterns{2}.name;

    summarize(subj,'objtype','pattern')

end