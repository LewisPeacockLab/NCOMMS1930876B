function[subj, args, nm] = mvpa_ph03_featselection(args, subj, nm, dirs)

%---------------------------------------------------------------------------
%*************** creating cross-validation indices
%*************** feature selection anova.
%*************** options: shiftTRs, peakWindow, featSelThresh
%---------------------------------------------------------------------------

xph = args.xphase;

%% ============= CREATING THE CROSS-VALIDATION INDICES
%*************** n-minus-one cross-validation schema for (within localizer only)
if args.cross_valid
    fprintf('\n(+) creating the cross-validation indices\n');
    
    if strcmp(args.regress_type, 'shift')
        nm.run_xvalid{xph} = sprintf('%s_norest_xval', nm.runs{xph});
        subj               = create_xvalid_indices(subj, nm.runs{xph}, ...
                             'actives_selname', nm.conds_sh_norest{xph}, ...
                             'new_selstem', nm.run_xvalid{xph});
                         
    elseif strcmp(args.regress_type, 'beta')
        nm.run_xvalid{xph} = sprintf('%s_xval', nm.runs{xph});
        subj               = create_xvalid_indices(subj, nm.runs{xph}, ...
                             'new_selstem', nm.run_xvalid{xph});
    end
else
    %*************** set regressors with 1 for all localizer runs
    nm.train_selector  = sprintf('%s_norest_train', nm.runs{xph});
    subj               = init_object(subj,'selector', nm.train_selector);
    train_selector_sh  = get_mat(subj,'selector', nm.conds_sh_norest{xph});
    subj               = set_mat(subj,'selector', nm.train_selector, train_selector_sh);
end

summarize(subj,'objtype','selector')

%% ============= FEATURE SELECT ANOVA
%*************** voxel-wise ANOVA: p value for conditions in each voxel over the time
fprintf('\n(+) feature selection using ANOVA\n');

if args.featSel
    %% *************** setup names
    if args.featVox
        nm.pats_z_stats{xph} = sprintf('%s_stat', nm.pats_z{xph});
        nm.mask_selec{xph}   = sprintf('%s_%s_%svox', nm.pats_z_stats{xph}, ...
            args.epi_name, num2str(args.fsVosNum));%for ANOVA-ed mask
    else
        nm.mask_selec{xph}   = sprintf('%s_%s_thresh%s', nm.pats_z{xph}, ...
            args.epi_name, args.featSelThresh);%for ANOVA-ed mask
    end
    
    %*************** reset mask name
    if args.class_selecting
        nm.mask_selec{xph} = sprintf('cate%s_%s', sprintf('%d',args.selected_category), nm.mask_selec{xph});
    end
    
    if strcmp(args.regress_type,'shift')
        if (args.xphase==3) % operation
            mask_basename = sprintf('%s_%s_sh%d_blk_%s', args.mask_name, ...
                nm.mask_selec{xph}, args.shift_TRs, args.rest);
        else
            mask_basename = sprintf('%s_%s_%s_sh%d_blk_%s', args.mask_name, ...
                nm.mask_selec{xph}, args.level, args.shift_TRs, args.rest);
        end
    elseif strcmp(args.regress_type,'beta')
        mask_basename = sprintf('%s_%s_%s', args.mask_name, nm.mask_selec{xph}, args.level);
    end
    
    if args.cross_valid% localizer
        n_validation = numel(unique(args.regs{xph}.selectors));
        for i=1:n_validation
            nm.mask_selected{xph}{i}      = sprintf('%s_%d', nm.mask_selec{xph}, i);
            nm.mask_selected_file{xph}{i} = sprintf('%s_%d', mask_basename, i);
        end
    else% decoding
        n_validation                  = 1;%1 cross-validation indices for training run (to decode testing run)
        nm.mask_selected{xph}{1}      = nm.mask_selec{xph};
        nm.mask_selected_file{xph}{1} = mask_basename;
    end
    
    %% *************** check if feature selected mask already exists
    disp('... starting feature selection ANOVA on data')
    
    if     strcmp(args.regress_type, 'shift'), nm_regs = nm.conds_sh{xph};
    elseif strcmp(args.regress_type, 'beta'),  nm_regs = nm.conds{xph}; end
    
    if args.featVox
        if args.cross_valid
            %%%%%%% use selector with cross-validation indices_tag #
            subj = feature_select(subj, nm.pats_z{xph}, nm_regs, nm.run_xvalid{xph}, ...
                'new_map_patname', nm.pats_z_stats{xph},...
                'thresh', []);
        else
            %%%%%%% use train selector: all 1
            subj = peek_feature_select(subj, nm.pats_z{xph}, nm_regs, nm.train_selector, ...
                'new_map_patname', nm.pats_z_stats{xph},...
                'thresh', []);
        end
        
        %*************** find N best_voxels
        subj = create_sorted_mask(subj, nm.pats_z_stats{xph}, ...
            nm.mask_selec{xph}, args.fsVosNum,...
            'descending',true);
        
    else
        if args.cross_valid
            %%%%%%% use selector with cross-validation indices_tag #
            subj = feature_select(subj, nm.pats_z{xph}, nm_regs, nm.run_xvalid{xph}, ...
                'thresh', str2double(args.featSelThresh));
            
        else
            %%%%%%% use train selector: all 1
            subj = peek_feature_select(subj, nm.pats_z{xph}, nm_regs, nm.train_selector, ...
                'thresh', str2double(args.featSelThresh));
        end
    end
    
    %*************** write feature-selection mask to disk
    for i=1:n_validation
        if strcmp(args.regress_type, 'shift')
            selsected_mask_dir = dirs.mvpa.selected_mask{xph};
        elseif strcmp(args.regress_type, 'beta')
            selsected_mask_dir = dirs.rsa.selected_mask{xph};
        end
        
        fprintf('... writing feature-selection mask: %s/%s\n',...
            selsected_mask_dir, nm.mask_selected_file{xph}{i});
        
        write_to_spm(subj,'mask', subj.masks{i+1}.name, ...
            'output_filename', fullfile(selsected_mask_dir, ...
            nm.mask_selected_file{xph}{i}));
        
        xselected_mask = fullfile(dirs.mvpa.selected_mask{xph}, ...
            sprintf('%s.%s',nm.mask_selected_file{xph}{i}, args.epiext));
        
        fprintf('... zipping mask.nii.gz');
        gzip(xselected_mask, dirs.mvpa.selected_mask{xph});
        
        fprintf('... deleting mask.nii');
        delete(xselected_mask);
    end
    
    args.anova = 1;
    
    summarize(subj,'display_groups', false)
    summarize(subj,'objtype','mask')
    
    %% *************** display size of peek mask
    mask_obj      = get_object(subj,'mask', nm.mask{xph});
    
    for i=1:n_validation
        feat_mask_obj = get_object(subj,'mask', subj.masks{i+1}.name);
        peeksize      = feat_mask_obj.nvox;
        
        fprintf('\n... phase_%s feature selection: %d voxels (%.1f %% of %d voxels)\n',...
            args.phase_name{xph}, peeksize, 100*peeksize/mask_obj.nvox, mask_obj.nvox);
    end
    
else
    nm.mask_selected{xph}{1} = nm.mask{xph};
end

end