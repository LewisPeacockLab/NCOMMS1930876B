function[] = clearmem_localizer_importance_map(args, dirs)

%% ============= UNPACK ARGS.
xph               = args.xphase;
mask_name         = args.mask_name;
args.regress_type = args.train_regress_type;
n_regs            = length(args.regs{xph}.regressor_name);
subject_list      = args.subject_list;
n_subs            = length(args.g_sub);
imp_type          = {'pos','neg','overlap'};
reg_list          = 1:3;
xgrp_subs         = args.g_sub;%args.filtered_subs;

%*************** output directory
xpeak_dir         = fullfile(dirs.mvpa.group.imp_map{xph}, ...
    sprintf('top_%s_%s', num2str(args.peak_thresh * 100), args.rest));

if ~isdir(xpeak_dir), mkdir(xpeak_dir); end

%*************** output basename
basename          = args.analysis_basename;

%*************** reference volume
xnorm_refer_epi   = fullfile(dirs.mvpa.imp_map{xph},'norm_impmap_localizer_subcategory_bold_avg_mcf_brain_mask_rest_mcduff_mean_cond1_actors_pos.nii.gz');

%*************** set FSL environment
setenv('FSLDIR', dirs.fsl);  % this to tell where FSL folder is
setenv('FSLOUTPUTTYPE', 'NIFTI_GZ'); % this to tell what the output type would be
        
%% ============= SETUP FILE NAMES
%*************** ph1. base filename
ph1.basename = sprintf('%s_%s_zscored_%s', args.phase_name{xph}, args.mask_name, args.epi_name); 

%*************** ph2. base filename
if strcmp(args.regress_type, 'shift')
    ph2.basename = sprintf('%s_%s_%s%dtr_blk_%s',...
                   ph1.basename, args.level, args.regress_type, ...
                   args.shift_TRs, args.rest);
elseif strcmp(args.regress_type, 'beta')
    ph2.basename = sprintf('%s_%s_%s', ...
                   ph1.basename, args.level, args.regress_type);
end

%*************** reset ph2. filename
if args.class_selecting
    ph2.basename = sprintf('cate%s_%s', sprintf('%d',args.selected_category), ph2.basename);
    args.penalty_check{xph} = 0;
end

%*************** ph3. base filenames
if args.featVox
    ph3.basename = sprintf('%s_featsel_%svox', ph2.basename, num2str(args.fsVosNum));
else
    ph3.basename = sprintf('%s_featsel_thresh%s', ph2.basename, num2str(args.featSelThresh));
end

%*************** ph4. base filenames
class_basename  = sprintf('classified_%s_%s', ph3.basename, args.classifier);

args.grp_results_name = class_basename;

%*************** reset group name
if args.class_selecting
    grp_name = sprintf('grp_imp_map_%s_cate_%s', ...
        args.imp_type, basename);
else
    grp_name = sprintf('grp_imp_map_%s_%s', args.imp_type, basename);
end

grp_fname = sprintf('%s/%s.mat', dirs.mvpa.group.imp_map{xph}, grp_name);

%% ============= 2ND LEVEL SAVE
if args.group_mat{xph}
    
    if exist(grp_fname, 'file')
        fprintf('(+) load 2nd level normalized importance map\n');
        load(grp_fname); 
    end %grp_norm_pattern
    
    %% ============= 1ST LEVEL SUBJECT DATA
    for xsub = args.g_sub
        %*************** setup subject & directories
        args.subject_id = subject_list(xsub).name;
        dirs            = setup_directory(dirs, args);
        
        if ~isdir(dirs.mvpa.imp_map{xph}), mkdir(dirs.mvpa.imp_map{xph}); end
        
        %% ============= load penalty_check
        %*************** ph2. base filename
        if strcmp(args.regress_type, 'shift')
            penalty_ph2.basename = sprintf('%s_%s_%s%dtr_blk_%s',...
                ph1.basename, args.level, args.regress_type, ...
                args.shift_TRs, args.rest);
        elseif strcmp(args.regress_type, 'beta')
            penalty_ph2.basename = sprintf('%s_%s_%s', ...
                ph1.basename, args.level, args.regress_type);
        end
        
        %*************** reset ph2. filename
        if args.class_selecting
            penalty_ph2.basename = sprintf('cate%s_%s', sprintf('%d',args.selected_category), penalty_ph2.basename);
        end
        
        %*************** ph3. base filename
        if args.featVox
            penalty_ph3.basename = sprintf('%s_featsel_%svox', penalty_ph2.basename, num2str(args.fsVosNum));
        else
            penalty_ph3.basename = sprintf('%s_featsel_thresh%s', penalty_ph2.basename, num2str(args.featSelThresh));
        end
        
        %*************** basename phase 4.
        penalty_basename  = sprintf('classified_%s_%s', penalty_ph3.basename, args.classifier);
        
        %*************** load penalty
        loc_pen = load(sprintf('%s/penalty_check_%s.mat', dirs.mvpa.output{xph}, penalty_basename));%'pen_check'
        [xacc, whichmax] = max(loc_pen.pen_check.performance); %#ok<*NODEF>
        max_penalty      = loc_pen.pen_check.penalty(whichmax);
        args.penalty     = max_penalty;
        
        fprintf('... max_penalty: %s: acc: %1.4f\n', num2str(max_penalty), xacc);
        
        %*************** definde penalty + classification
        args.xpenalty = min(max_penalty);
        fprintf('... %s max_penalty: %s\n', args.subject_id, num2str(min(max_penalty)));
        
        %*************** basename phase 4.
        ph4.basename  = sprintf('%s_penalty%s', class_basename, num2str(args.xpenalty));
        
        %*************** load phase 4.
        fname = sprintf('%s/%s.mat', dirs.mvpa.output{xph}, ph4.basename);
        load(fname);%'ph4'
        
        fprintf('\n... loaded classification results (penalty: %s) of %s: %s\n', ...
            num2str(args.xpenalty), args.subject_id, fname);
        
        %% ============= 1ST-LEVEL IMPORTANCE MAP
        xdir = fullfile(dirs.epi_mid, 'warp_param');
        if ~isdir(xdir), mkdir(xdir); end
        
        xname = sprintf('subj_imp_map_%s_%s', args.imp_type, basename);
        
        if args.class_selecting
            xname = sprintf('subj_imp_map_%s_cate%s_%s', ...
                args.imp_type, sprintf('%d',args.selected_category), basename);
        end
        
        if strcmp(args.rest, 'rest')
            xname = sprintf('%s_%s', xname, args.rest);
        end
        
        fname = sprintf('%s/%s.mat', dirs.mvpa.imp_map{xph}, xname);
        
        %*************** running 1st level importance map
        if args.subj_impmap{xph}
            clear subj
            
            fprintf('(+) 1st level importance map: %s\n', args.subject_id);
            
            [subj] = create_importance_maps_logreg(ph4.subj, ph4.results, args, dirs);
            
            save(fname, 'subj','-v7.3');
        else
            load(fname);% subj
        end
        
        %% ============= SAVE STANDARDIZED PATTERNS
        %%************** grp_norm_pattern{xsub}{xreg}{pos/neg}
        n_ori_pats = length(ph4.subj.patterns);
        n_iters    = length(ph4.results.iterations);
        n_patterns = n_ori_pats + (n_iters * n_regs);

        for xreg = 1:n_regs
            for i = 1:2%pos,neg
                xunit = n_patterns + i + (2 * (xreg-1));
                grp_norm_pattern{xsub}{xreg}{i} = subj.patterns{xunit}.norm.pat_mat; %#ok<*AGROW,*NASGU>
            end
        end    
    end
    
    %% ============= 2ND LEVEL SAVE
    fprintf('(+) save normalized importance map pattern in group structure\n\n');
   
    save(grp_fname, 'grp_norm_pattern','-v7.3');    
    
else
    %% ============= LOAD GROUP NORMED IMPORTANCE MAP
    %% *************** setup subject & directories
    args.subject_id = subject_list(1).name;
    dirs            = setup_directory(dirs, args);
    
    fprintf('(+) load normalized importance map pattern from group structure\n\n');

    load(grp_fname);%'grp_norm_pattern'
    
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%============= 2ND LEVEL IMPORTANCE MAP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create_group_norm_importance_map

if args.class_selecting
    impmap_group = sprintf('grp_impmap_%s_%s_cate%s_%s_%s_%s_mean_cond', ...
        args.phase_name{args.xphase}, args.level, ...
        sprintf('%d',args.selected_category), args.mask_name, args.rest, args.imp_type);
else
    impmap_group = sprintf('grp_impmap_%s_%s_%s_%s_%s_mean_cond', ...
        args.phase_name{args.xphase}, args.level, args.mask_name, args.rest, args.imp_type);
end

xname = sprintf('subj_grp_imp_map_%s_%s', args.imp_type, basename);
if args.class_selecting
    xname = sprintf('subj_grp_imp_map_%s_cate%s_%s', ...
        args.imp_type, sprintf('%d',args.selected_category), basename);
end

grp_subj_fname = sprintf('%s/%s.mat', dirs.mvpa.group.imp_map{xph}, xname);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%============= 2ND LEVEL MEAN IMPORTANCE MAP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if args.group_impmap{xph}    
    %% ============= initiate group subject structure
    %*************** standard mni whole brain structure
    
    fprintf('(+) 2nd level group mean importance map\n');  
    
    clear grp_subj
    
    grp_subj     = init_subj(sprintf('%s_impmap', args.experiment), 'group');%identifier of the subj
    grp_subj     = load_spm_mask_gz(grp_subj, 'mni_mask', args.mni_brain);
    xmask        = get_mat(grp_subj, 'mask', 'mni_mask');
    
    xmask_cord   = find(xmask);
    n_masked_vox = size(xmask_cord, 1);
    
    %% ============== mean patterns
    for xreg = 1:n_regs
        %*************** build in pattern structure in grp_subj
        
        xcond_name = args.regs{args.xphase}.regressor_name{xreg};
        
        %% ============= positive | negative
        for i = 1:2%pos,neg
            clear t_pattern xpattern xgrp_mean_imp_map xvol_pat
            
            fprintf('... cond%d/%d_%s...\n', xreg, n_regs, imp_type{i});  
            
            mean_impmap_name = sprintf('%s%d_%s_%s', impmap_group, xreg, xcond_name, imp_type{i});
            
            %*************** get patterns
            xpattern = zeros(n_masked_vox, n_subs);
            
            for xsub = xgrp_subs %args.g_sub
                t_pattern = grp_norm_pattern{xsub}{xreg}{i}; %#ok<*AGROW,*NASGU>
                xpattern(:, xsub) = t_pattern(xmask_cord);
            end
            
            xgrp_mean_impmap = zeros(n_masked_vox, 1);
            for xvox = 1:n_masked_vox
                tpat = xpattern(xvox, :);
                if sum(tpat)~=0
                    xgrp_mean_impmap(xvox, 1) = mean(tpat(tpat~=0));
                end
            end
            
            %*************** new volume
            xvol_pat = zeros(size(xmask));
            xvol_pat(xmask_cord) = xgrp_mean_impmap;
            
            %*************** grp_subj structure
            grp_subj = init_object(grp_subj,'pattern', mean_impmap_name);
            grp_subj = set_mat(grp_subj,'pattern', mean_impmap_name, xgrp_mean_impmap);
            
            grp_subj = set_objfield(grp_subj, 'pattern', mean_impmap_name, 'masked_by', args.mni_brain);
            grp_subj = set_objfield(grp_subj, 'pattern', mean_impmap_name, 'group_name', impmap_group);
            grp_subj = set_objfield(grp_subj, 'pattern', mean_impmap_name, 'pat_mat', xvol_pat, 'ignore_absence', true);
            
            mean_new_filename = fullfile(dirs.mvpa.group.imp_map{xph}, ...
                sprintf('%s.nii', mean_impmap_name));
            
            %*************** create nii.gz
            [~, ~, refer_vol] = icatb_read_gzip_nii(xnorm_refer_epi);
            
            xcur_vol       = refer_vol;%grp_subj.masks{1}.header.vol;
            xcur_vol.fname = mean_new_filename;
            spm_write_vol(xcur_vol, xvol_pat);
            
            gzip(mean_new_filename, dirs.mvpa.group.imp_map{xph});
            delete(mean_new_filename);
            
            %*************** change orientation
            system(sprintf('%s/bin/fslorient -copysform2qform %s', ...
                dirs.fsl, sprintf('%s.gz', mean_new_filename)));
            
            %% ============= threshold the results
            if args.peak_thresh==0.1 % top 10%
                fprintf('... thresholding top %s: cond%d/%d_%s...\n', ...
                    num2str(args.peak_thresh * 100), xreg, n_regs, imp_type{i});
            else
                fprintf('... thresholding std %s SD: cond%d/%d_%s...\n', ...
                    num2str(args.peak_thresh), xreg, n_regs, imp_type{i});
            end
            
            xmean = mean(xgrp_mean_impmap);
            xsd   = std(xgrp_mean_impmap);
            
            if i==1
                if args.peak_thresh==0.1 % top 10%
                    sorted_Y   = sort(xgrp_mean_impmap, 'descend');
                    xcriterion = sorted_Y(round(size(xgrp_mean_impmap, 1) * args.peak_thresh));
                    
                else % above mean + 2 sd
                    xcriterion = xmean + (xsd * args.peak_thresh);
                end
            end
            
            %*************** cutoff table
            xunit = xgrp_mean_impmap >= xcriterion;
            
            xheader = {'condition','map','mean','sd','cutoff',...
                'total_voxels','select_voxels','percent','max','min'};
            xarray{i + (3 * (xreg-1)), 1}  = xcond_name;
            xarray{i + (3 * (xreg-1)), 2}  = imp_type{i};
            xarray{i + (3 * (xreg-1)), 3}  = xmean;
            xarray{i + (3 * (xreg-1)), 4}  = xsd;
            xarray{i + (3 * (xreg-1)), 5}  = xcriterion;
            xarray{i + (3 * (xreg-1)), 6}  = size(xgrp_mean_impmap,1);
            xarray{i + (3 * (xreg-1)), 7}  = size(find(xunit),1);
            xarray{i + (3 * (xreg-1)), 8}  = (size(find(xunit),1)/size(xgrp_mean_impmap,1)) * 100;
            xarray{i + (3 * (xreg-1)), 9}  = max(xgrp_mean_impmap(xunit));
            xarray{i + (3 * (xreg-1)), 10} = min(xgrp_mean_impmap(xunit));
            
            %*************** new volume
            xpeak_grp_mean_impmap = zeros(size(xgrp_mean_impmap));
            xpeak_grp_mean_impmap(xunit) = xgrp_mean_impmap(xunit);
            
            xvol_pat_peak = zeros(size(xmask));
            xvol_pat_peak(xmask_cord) = xpeak_grp_mean_impmap;
            
            grp_subj = set_objfield(grp_subj, 'pattern', mean_impmap_name, 'peak_pat_mat', xpeak_grp_mean_impmap, 'ignore_absence', true);
            grp_subj = set_objfield(grp_subj, 'pattern', mean_impmap_name, 'peak_pat_vol', xvol_pat_peak, 'ignore_absence', true);
            grp_subj = set_objfield(grp_subj, 'pattern', mean_impmap_name, 'peak_top', args.peak_thresh, 'ignore_absence', true);
            
            if args.peak_thresh==0.1 % top 10%
                peak_mean_new_filename = fullfile(xpeak_dir, ...
                    sprintf('peak_top%s_%s.nii', num2str(args.peak_thresh*100), mean_impmap_name));
            else 
                peak_mean_new_filename = fullfile(xpeak_dir, ...
                    sprintf('peak_sd%s_%s.nii', num2str(args.peak_thresh), mean_impmap_name));
            end
            
            xcur_vol.fname = peak_mean_new_filename;
            spm_write_vol(xcur_vol, xvol_pat_peak);
            
            gzip(peak_mean_new_filename, xpeak_dir);
            delete(peak_mean_new_filename);
            
            %*************** change orientation
            system(sprintf('%s/bin/fslorient -copysform2qform %s', ...
                dirs.fsl, sprintf('%s.gz', peak_mean_new_filename)));
            
        end
        
        %% ============= overlapped positive & negative
        
        fprintf('... overlapping b/w pos & neg: cond%d/%d...\n', xreg, n_regs);  
            
        xpos_mat   = get_objfield(grp_subj, 'pattern', sprintf('%s%d_%s_%s', impmap_group, xreg, xcond_name, imp_type{1}), 'peak_pat_mat');
        xneg_mat   = get_objfield(grp_subj, 'pattern', sprintf('%s%d_%s_%s', impmap_group, xreg, xcond_name, imp_type{2}), 'peak_pat_mat');
        
        xpos_index = find(xpos_mat);
        xneg_index = find(xneg_mat);
        
        overlap_index = intersect(xpos_index, xneg_index);
        
        %*************** cutoff table
        i = 3;
        xarray{i + (3 * (xreg-1)), 1} = xcond_name;
        xarray{i + (3 * (xreg-1)), 2} = imp_type{i};
        xarray{i + (3 * (xreg-1)), 3} = xmean;
        xarray{i + (3 * (xreg-1)), 4} = xsd;
        xarray{i + (3 * (xreg-1)), 5} = xcriterion;
        xarray{i + (3 * (xreg-1)), 6} = size(xgrp_mean_impmap,1);
        xarray{i + (3 * (xreg-1)), 7} = size(overlap_index,1);
        xarray{i + (3 * (xreg-1)), 8} = size(overlap_index,1)/size(xgrp_mean_impmap,1);
        
        %*************** overlaping volume
        xgrp_overlap_impmap = zeros(n_masked_vox, 1);
        xgrp_overlap_impmap(overlap_index) = 1;%(xpos_mat(overlap_index) - xneg_mat(overlap_index));
        
        %*************** new volume
        xvol_pat_overlap = zeros(size(xmask));
        xvol_pat_overlap(xmask_cord) = xgrp_overlap_impmap;
        
        mean_impmap_name = sprintf('%s%d_%s_%s', impmap_group, xreg, xcond_name, imp_type{1});%in pos
        
        grp_subj = set_objfield(grp_subj, 'pattern', mean_impmap_name, 'peak_overlap_mat', xgrp_overlap_impmap, 'ignore_absence', true);
        grp_subj = set_objfield(grp_subj, 'pattern', mean_impmap_name, 'peak_overlap_vol', xvol_pat_overlap, 'ignore_absence', true);
        
        %*************** save new volume
        overlap_impmap_name = sprintf('%s%d_%s_%s', impmap_group, xreg, xcond_name, imp_type{3});
        
        if args.peak_thresh==0.1
            peak_overlap_new_filename = fullfile(xpeak_dir, ...
                sprintf('peak_top%s_%s.nii', num2str(args.peak_thresh * 100), overlap_impmap_name));
        else
            peak_overlap_new_filename = fullfile(xpeak_dir, ...
                sprintf('peak_sd%s_%s.nii', num2str(args.peak_thresh), overlap_impmap_name));
        end
        
        xcur_vol.fname = peak_overlap_new_filename;
        spm_write_vol(xcur_vol, xvol_pat_overlap);
        
        gzip(peak_overlap_new_filename, xpeak_dir);
        delete(peak_overlap_new_filename);
        
        %*************** change orientation
        system(sprintf('%s/bin/fslorient -copysform2qform %s', ...
            dirs.fsl, sprintf('%s.gz', peak_overlap_new_filename)));
        
    end
    
    fprintf('\n\n');
    
    %% ============= SAVE grp_subj
    
    save(grp_subj_fname, 'grp_subj', '-v7.3');

    %% *************** voxel selection table
    xtable = cell2table(xarray, 'VariableNames', xheader);
    
    if args.peak_thresh==0.1
        writetable(xtable, sprintf('%s/voxel_select_top%s_%s.csv', ...
            dirs.mvpa.group.imp_map{xph}, num2str(args.peak_thresh * 100), basename));
    else
        writetable(xtable, sprintf('%s/voxel_select_sd%s_%s.csv', ...
            dirs.mvpa.group.imp_map{xph}, num2str(args.peak_thresh), basename));
    end
    
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%============= 2ND LEVEL DIFFERENCES OF IMPORTANCE MAP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% not including rest category for both "rest" and "norest" classifiers

if args.group_impmap_diff{xph}
    
    fprintf('(+) 2nd level group difference of importance map\n'); 
    
    %*************** load grp_subj structure
    if ~(args.group_impmap{xph}), load(grp_subj_fname); end %'grp_subj' 
    
    %*************** get mask coordination
    xmask        = get_mat(grp_subj, 'mask', 'mni_mask');
    
    xmask_cord   = find(xmask);
    n_masked_vox = size(xmask_cord, 1);
    
    %% ============= get mat from grp_subj
    clear xgrp_mean_impmap
    
    for xreg = 1:n_regs
        for i = 1:2%pos,neg
            
            fprintf('cond%d/%d_%s...\n', xreg, n_regs, imp_type{i});
            
            xcond_name = args.regs{args.xphase}.regressor_name{xreg};
            
            mean_impmap_name = sprintf('%s%d_%s_%s', impmap_group, xreg, xcond_name, imp_type{i});
            
            xgrp_mean_impmap{xreg}{i} = get_mat(grp_subj,'pattern', mean_impmap_name);
            
        end
    end
    
    fprintf('\n');
        
    %% ============= differences [1 -1/2 -1/2]
    
    diff_impmap_group = sprintf('diff_%s', impmap_group);
            
    for xreg = 1:n_regs
        for i = 1:2%pos,neg
            
            clear xgrp_diff_impmap nontarg xvol_pat xcur_vol sorted_Y xvol_pat_peak
            
            fprintf('... cond%d/%d_%s...\n', xreg, n_regs, imp_type{i});  
            
            xcond_name = args.regs{args.xphase}.regressor_name{xreg};
            
            %*************** mean of the non target importance map
            nontarg = reg_list(~ismember(reg_list, xreg));
            
            nontarg_impmap = zeros(size(xgrp_mean_impmap{1}{i}));
            for xntarg = 1:length(nontarg)
                nontarg_impmap = nontarg_impmap + xgrp_mean_impmap{nontarg(xntarg)}{i};
            end
            
            xgrp_diff_impmap = xgrp_mean_impmap{xreg}{i} - (nontarg_impmap/length(nontarg));
            
            %*************** new volume
            xvol_pat = zeros(size(xmask));
            xvol_pat(xmask_cord) = xgrp_diff_impmap;
            
            %*************** setup grp_subj structures
            diff_impmap_name  = sprintf('%s%d_%s_%s', diff_impmap_group, xreg, xcond_name, imp_type{i});
            
            grp_subj = init_object(grp_subj,'pattern', diff_impmap_name);
            grp_subj = set_mat(grp_subj,'pattern', diff_impmap_name, xgrp_diff_impmap);
            
            grp_subj = set_objfield(grp_subj, 'pattern', diff_impmap_name, 'masked_by', args.mni_brain);
            grp_subj = set_objfield(grp_subj, 'pattern', diff_impmap_name, 'group_name', diff_impmap_group);
            grp_subj = set_objfield(grp_subj, 'pattern', diff_impmap_name, 'pat_vol', xvol_pat, 'ignore_absence', true);
            
            diff_mean_filename = fullfile(dirs.mvpa.group.imp_map{xph}, sprintf('%s.nii', diff_impmap_name));
            
            %*************** create nii.gz
            [~, ~, refer_vol] = icatb_read_gzip_nii(xnorm_refer_epi);
            
            xcur_vol       = refer_vol;%grp_subj.masks{1}.header.vol;
            xcur_vol.fname = diff_mean_filename;
            spm_write_vol(xcur_vol, xvol_pat);
            
            gzip(diff_mean_filename, dirs.mvpa.group.imp_map{xph});
            delete(diff_mean_filename);
            
            %*************** change orientation
            system(sprintf('%s/bin/fslorient -copysform2qform %s', ...
                dirs.fsl, sprintf('%s.gz', diff_mean_filename)));
            
            %% ============= threshold the results
            
            if args.peak_thresh==0.1 % top 10%
                fprintf('... thresholding top %s: cond%d/%d_%s...\n', ...
                    num2str(args.peak_thresh * 100), xreg, n_regs, imp_type{i});
            else
                fprintf('... thresholding std %s SD: cond%d/%d_%s...\n', ...
                    num2str(args.peak_thresh), xreg, n_regs, imp_type{i});
            end
            
            if i==1
                if args.peak_thresh==0.1 % top 10%
                    sorted_Y   = sort(xgrp_diff_impmap, 'descend');
                    xcriterion = sorted_Y(round(size(xgrp_diff_impmap, 1) * args.peak_thresh));
                    
                else % above mean + 2 sd
                    xmean      = mean(xgrp_diff_impmap);
                    xsd        = std(xgrp_diff_impmap);
                    xcriterion = xmean + (xsd * args.peak_thresh);
                end
            end
            
            xunit      = xgrp_diff_impmap >= xcriterion;
            
            xpeak_grp_diff_impmap = zeros(size(xgrp_diff_impmap));
            xpeak_grp_diff_impmap(xunit) = xgrp_diff_impmap(xunit);
            
            %*************** new volume
            xvol_pat_peak = zeros(size(xmask));
            xvol_pat_peak(xmask_cord) = xpeak_grp_diff_impmap;
            
            grp_subj = set_objfield(grp_subj, 'pattern', diff_impmap_name, 'peak_pat_mat', xpeak_grp_diff_impmap, 'ignore_absence', true);
            grp_subj = set_objfield(grp_subj, 'pattern', diff_impmap_name, 'peak_pat_vol', xvol_pat_peak, 'ignore_absence', true);
            grp_subj = set_objfield(grp_subj, 'pattern', diff_impmap_name, 'peak_top', args.peak_thresh, 'ignore_absence', true);
            
            if args.peak_thresh==0.1 % top 10%
                peak_diff_new_filename = fullfile(xpeak_dir, ...
                    sprintf('peak_top%s_%s.nii', num2str(args.peak_thresh*100), diff_impmap_name));
            else 
                peak_diff_new_filename = fullfile(xpeak_dir, ...
                    sprintf('peak_sd%s_%s.nii', num2str(args.peak_thresh), diff_impmap_name));
            end
            
            xcur_vol.fname = peak_diff_new_filename;
            spm_write_vol(xcur_vol, xvol_pat_peak);
            
            gzip(peak_diff_new_filename, xpeak_dir);
            delete(peak_diff_new_filename);
            
            %*************** change orientation
            system(sprintf('%s/bin/fslorient -copysform2qform %s', ...
                dirs.fsl, sprintf('%s.gz', peak_diff_new_filename)));
            
        end
        
        fprintf('\n\n');
        
        %% ============= overlapped positive & negative
        
        fprintf('... overlapping b/w pos & neg: cond%d/%d...\n', xreg, n_regs);  
        
        xpos_mat   = get_objfield(grp_subj, 'pattern', sprintf('%s%d_%s_%s', diff_impmap_group, xreg, xcond_name, imp_type{1}), 'peak_pat_mat');
        xneg_mat   = get_objfield(grp_subj, 'pattern', sprintf('%s%d_%s_%s', diff_impmap_group, xreg, xcond_name, imp_type{2}), 'peak_pat_mat');
        
        xpos_index = find(xpos_mat);
        xneg_index = find(xneg_mat);
        
        overlap_index = intersect(xpos_index, xneg_index);
        
        %*************** overlaping volume
        xgrp_overlap_impmap = zeros(n_masked_vox, 1);
        xgrp_overlap_impmap(overlap_index) = 1;%(xpos_mat(overlap_index) - xneg_mat(overlap_index));
        
        %*************** new volume
        xvol_pat_overlap = zeros(size(xmask));
        xvol_pat_overlap(xmask_cord) = xgrp_overlap_impmap;
        
        diff_impmap_name  = sprintf('%s%d_%s_%s', diff_impmap_group, xreg, xcond_name, imp_type{1});%in pos
        
        grp_subj = set_objfield(grp_subj, 'pattern', diff_impmap_name, 'peak_overlap_mat', xgrp_overlap_impmap, 'ignore_absence', true);
        grp_subj = set_objfield(grp_subj, 'pattern', diff_impmap_name, 'peak_overlap_vol', xvol_pat_overlap, 'ignore_absence', true);
        
        %*************** save new volume
        overlap_impmap_name = sprintf('%s%d_%s_%s', diff_impmap_group, xreg, xcond_name, imp_type{3});
        
        if args.peak_thresh==0.1
            peak_overlap_new_filename = fullfile(xpeak_dir, ...
                sprintf('peak_top%s_%s.nii', num2str(args.peak_thresh * 100), overlap_impmap_name));
        else
            peak_overlap_new_filename = fullfile(xpeak_dir, ...
                sprintf('peak_sd%s_%s.nii', num2str(args.peak_thresh), overlap_impmap_name));
        end
        
        xcur_vol.fname = peak_overlap_new_filename;
        spm_write_vol(xcur_vol, xvol_pat_overlap);
        
        gzip(peak_overlap_new_filename, xpeak_dir);
        delete(peak_overlap_new_filename);
        
        %*************** change orientation
        system(sprintf('%s/bin/fslorient -copysform2qform %s', ...
            dirs.fsl, sprintf('%s.gz', peak_overlap_new_filename)));
        
    end
    
    %% *************** re-save grp_subj structure
    save(grp_subj_fname, 'grp_subj', '-v7.3');
    
end

end