function[] = clearmem_operation_cs_subj(args, dirs)
% creating regressors for cs mvpa in python
xout_dir = dirs.param;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%============== concatenated nii.gz
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% first TRs were trimmed

for xph = [1 3]
    
    xphase_name = args.phase_name{xph};
    
    %% ============= regressors
    if xph==3
        args.index{xph} = create_design_index_study(args, dirs);
        
        if args.four_oper_regress
            args.regs{xph} = create_mvpa_regressor_operation_four(args, dirs);
        else
            args.regs{xph} = create_mvpa_regressor_operation(args, dirs);
        end
    else
        args.index{xph} = create_design_index_localizer(args, dirs);
        args.regs{xph}  = create_mvpa_regressor_localizer(args, dirs);
    end
    
    xidx      = args.index{xph}.matrix(findCol(args.index{xph}.header, {'vol_index'}), :);
    xsel_vols = find(xidx);
    
    xrun_id   = args.index{xph}.matrix(findCol(args.index{xph}.header, {'it_run'}), xsel_vols);
    xvol_idx  = xidx(xsel_vols);
%     xsels     = args.regs{xph}.selectors(xsel_vols);
    xregs     = args.regs{xph}.regressors_index(xsel_vols);
    xspikes   = args.index{xph}.matrix(findCol(args.index{xph}.header, {'spike'}), xsel_vols);
    
    xregress  = [xvol_idx; xrun_id; xregs; xspikes];
%     xtable    = array2table(xregress', 'VariableNames', {'vol_idx','run_idx','regressor','spike'});
%     writetable(xtable, fullfile(xout_dir, sprintf('%s_cs_regressor.csv', xphase_name)));
%     dlmwrite(fullfile(xout_dir, sprintf('%s_cs_regressor.txt', xphase_name)), xregress');
    
    save(fullfile(xout_dir, sprintf('%s_cs_regressor.mat', xphase_name)), 'xregress');
    clear xregress
    
    %% ============= EXTRACTING PATTERNS
%     %*************** setup names
%     nm.mask  = args.mask_name;
%     nm.pat   = sprintf('%s_mni_patterns', xphase_name);
%     nm.pat_z = sprintf('%s_z', nm.pat);% subj.patterns{2}.name;
%     nm.sel   = sprintf('%s_runs', xphase_name);% selectors object name
%     
%     %*************** initiate subj
%     subj = init_subj(args.experiment, args.subject_id);%identifier of the subj
%     
%     %*************** load mask
%     if args.wholebrain, mask_dir = dirs.mni_mask; end
%     
%     fprintf('\n(+) load mask: %s.%s\n', nm.mask, args.epiext);
%     
%     xmask_gz = fullfile(mask_dir, sprintf('%s.%s.gz', args.mask_name, args.epiext));
%     subj     = load_spm_mask_gz(subj, nm.mask, xmask_gz);
%     
%     %% ============= EPI PATTERN: subj.patterns{1}
%     fprintf('\n(+) load epi data under mask %s\n', args.mask_name);
%     
%     %*************** epi
%     for xrun = 1:length(dirs.runs{xph})
%         args.epi_names{xph}{xrun} = ...
%             fullfile(dirs.runs{xph}{xrun}, sprintf('%s.%s.gz', args.epi_name, args.epiext));
%     end
%     
%     fprintf('\n... loading epi patterns\n');
%     
%     subj = load_spm_pattern_gz(subj, nm.pat, nm.mask, args.epi_names{xph});
%     
%     %*************** remove original pattern: we will use z-scored patterns
%     xpat     = get_mat(subj, 'pattern', nm.pat);
%     xsel_pat = xpat(:, xsel_vols);
%     subj     = remove_mat(subj,'pattern', nm.pat);
%     subj     = set_mat(subj, 'pattern', nm.pat, xsel_pat);
%     
%     %% ============= ZSCORING: subj.patterns{2}
%     %*************** zsoring each voxel, include 'rest' in the zscoring
%     fprintf('\n(+) z-scoring data: voxel-wise\n');
%     
%     %*************** SELECTORS: subj.selectors{1}
%     subj = initset_object(subj,'selector', nm.sel, xsels);
%     subj = zscore_runs(subj, nm.pat, nm.sel);
%     
%     summarize(subj)
    
    %% ============= write nii 
%     xnii_name = sprintf('zpat_%s_trim_%s_%s', ...
%         xphase_name, args.epi_name, nm.mask);
%     
%     write_to_spm(subj, 'pattern', nm.pat_z, ...
%         'output_filename', fullfile(dirs.mvpa.cs.subj, xnii_name));
%     
%     fprintf('... zipping mask.nii.gz\n');
%     gzip(fullfile(dirs.mvpa.cs.subj, sprintf('%s.nii', xnii_name)), dirs.bold);
%     
%     fprintf('... deleting mask.nii\n');
%     delete(fullfile(dirs.mvpa.cs.subj, sprintf('%s.nii', xnii_name)));
%     clear subj
    
end

end