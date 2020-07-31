function[subj, nm] = mvpa_ph01_patterns(args, subj, nm, dirs)
    
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
        if args.wholebrain, mask_dir = dirs.epi_mid;
        else                mask_dir = dirs.mask; end
    end
    
    %*************** MASK
    nm.mask{xph} = args.mask_name;
    fprintf('\n(+) load mask: %s.%s\n', args.mask_name, args.epiext);
    
    xmask = fullfile(mask_dir, sprintf('%s.%s', args.mask_name, args.epiext));
    if ~exist(xmask,'file'), gunzip(sprintf('%s.gz', xmask), mask_dir); end
        
    subj         = load_spm_mask(subj, nm.mask{xph}, xmask);
    mask_obj     = get_object(subj, 'mask', nm.mask{xph}); %'mask object' contents
    subj.masks{1}.header.description = 'raw pattern under the mask';
    delete(xmask);
    
    %% ============= EPI PATTERN: subj.patterns{1}
    fprintf('\n(+) load epi data under mask %s with %s voxels\n', ...
        args.mask_name, num2str(mask_obj.nvox));
    
    nm.pats{xph} = sprintf('%s_patterns', args.phase_name{xph});

    %*************** shift epi
    if strcmp(args.regress_type, 'shift')
        for xrun = 1:args.index{xph}.param.n_runs
            
            xepi = fullfile(dirs.runs{xph}{xrun}, sprintf('%s.%s.gz', args.epi_name, args.epiext));
            
            %*************** unzipped epis
            tmp_dir = sprintf('%s/tmp_%s_%s', dirs.runs{xph}{xrun}, args.level, args.rest);
            if ~isdir(tmp_dir), mkdir(tmp_dir); end
            
            args.epi_names{xph}{xrun} = ...
                fullfile(tmp_dir, sprintf('%s.%s', args.epi_name, args.epiext));
            
            %*************** unzip nii.gz
            fprintf('... unzip epi: %s\n', args.epi_names{xph}{xrun});
            gunzip(xepi, tmp_dir);
            
        end
        
        fprintf('\n');
        
        fprintf('... loading epi patterns\n');
        subj     = load_spm_pattern(subj, nm.pats{xph}, nm.mask{xph}, args.epi_names{xph});
        
        for xrun = 1:numel(args.epi_names{xph})
            tmp_dir = sprintf('%s/tmp_%s_%s', dirs.runs{xph}{xrun}, args.level, args.rest);
            rmdir(tmp_dir, 's')
        end
        
    elseif strcmp(args.regress_type, 'beta')
        args.epi_names = fullfile(dirs.beta.home, ...
                         sprintf('LSS_beta_localizer_%s.nii', args.subject_id));
        %*************** unzip nii.gz
        if ~exist(args.epi_names, 'file')
            fprintf('... unzip beta nii: %s\n\n', args.epi_names);
            gunzip(sprintf('%s.gz', args.epi_names), dirs.beta.home);
        end
        
        fprintf('... loading beta patterns\n');
        subj     = load_spm_pattern(subj, nm.pats{xph}, nm.mask{xph}, args.epi_names);
        
        delete(args.epi_names);
    end
    
    %% ============= RESET PATTERNS
    %*************** trim the first 10TRs of each run if untrimed 
    if ~(args.trim_trs)% 0_non_trimmed, 1_trimmed
        
        if xph == 2, xcell = 4; else xcell = 1; end
        
        %*************** reset patterns
        fprintf('\n  ... reset patterns: trim the first %sTRs: from subj.patterns{%d}.mat\n', ...
            num2str(args.xtrim), xcell);

        xmatrix = args.regs{xph}.matrix;
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