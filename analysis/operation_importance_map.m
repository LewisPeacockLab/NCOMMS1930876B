function[] = clearmem_operation_importance_map(args, dirs)
% step 6: importance map 1st + 2nd
% group-level
%   1. average all z-scored individual maps into group map
%   2. cutoff above 1 standard deviation

%% ============= UNPACK ARGS.
xph               = args.xphase;
n_regs            = length(args.regs{xph}.regressor_name);
subject_list      = args.subject_list;
imp_type          = {'pos','neg','comb'};
reg_list          = 1:n_regs;
xsub_groups       = args.filtered_subs;
n_imp_type        = length(imp_type);

%*************** output directory
xpeak_dir         = fullfile(dirs.mvpa.group.imp_map{xph}, ...
    sprintf('top_%s_%s', num2str(args.peak_thresh * 100), args.rest));

if ~isdir(xpeak_dir), mkdir(xpeak_dir); end %#ok<*ISDIR>

%*************** output basename
basename          = args.analysis_basename;

%*************** reference volume
xnorm_refer_epi   = fullfile(dirs.mvpa.imp_map{xph},'norm_impmap_operation_category_bold_avg_mcf_brain_mask_norest_mcduff_mean_cond1_maintain_pos.nii.gz');

%*************** set FSL environment
setenv('FSLDIR', dirs.fsl);  % this to tell where FSL folder is
setenv('FSLOUTPUTTYPE', 'NIFTI_GZ'); % this to tell what the output type would be

%% ============= SETUP FILE NAMES
ph1.basename = sprintf('%s_%s_%s', args.phase_name{xph}, args.mask_name, args.epi_name);
ph2.basename = sprintf('%s_%s_tr%s_blk_%s',...
    ph1.basename, args.regress_type, ...
    num2str(args.shift_TRs), args.rest);
ph3.basename = sprintf('%s_featsel_thresh%s', ph2.basename, num2str(args.featSelThresh));
ph4.basename = sprintf('%s_decoding_setup', ph3.basename);

%*************** basename phase 5.
class_basename   = sprintf('classified_%s_%s', ph4.basename, args.classifier);
penalty_basename = sprintf('classified_%s_%s', ph3.basename, args.classifier);

args.grp_results_name = class_basename;

%*************** reset group name
grp_name  = sprintf('grp_imp_map_%s_%s', args.imp_type, basename);
grp_fname = sprintf('%s/%s.mat', dirs.mvpa.group.imp_map{xph}, grp_name);

%% ============= 2ND LEVEL SAVE
%%============= 1ST LEVEL SUBJECT DATA
for xsub = xsub_groups
    %*************** setup subject & directories
    args.subject_id = subject_list(xsub).name;
    dirs            = setup_directory(dirs, args);
    
    if ~isdir(dirs.mvpa.imp_map{xph}), mkdir(dirs.mvpa.imp_map{xph}); end
    
    fprintf('... sub: %s\n', args.subject_id);
    
    %% ============= load penalty_check
    load(sprintf('%s/penalty_check_%s.mat', dirs.mvpa.output{xph}, penalty_basename));%#ok<*LOAD> %'pen_check'
    [xacc, whichmax] = max(pen_check.performance); %#ok<*NODEF>
    max_penalty      = pen_check.penalty(whichmax);
    args.xpenalty    = max_penalty;
    
    fprintf('... max_penalty: %s: acc: %1.4f\n', num2str(max_penalty), xacc);
    
    %*************** load phase 4
    fname        = sprintf('%s/ph4_%s.mat', dirs.mvpa.scratch{xph}, ph4.basename);
    load(fname);%'ph4'
    
    %*************** load phase 5:results
    ph5.basename = sprintf('%s_penalty%s', class_basename, num2str(args.xpenalty));
    fname        = sprintf('%s/%s.mat', dirs.mvpa.output{xph}, ph5.basename);
    load(fname);%'ph5'  
    
    fprintf('\n... loaded classification results (penalty: %s) of %s: %s\n', ...
        num2str(args.xpenalty), args.subject_id, fname);
    
    %% ============= 1ST-LEVEL IMPORTANCE MAP
    xdir = fullfile(dirs.epi_mid, 'warp_param');
    if ~isdir(xdir), mkdir(xdir); end
    
    xname = sprintf('subj_imp_map_%s_%s', args.imp_type, basename);
    fname = sprintf('%s/%s.mat', dirs.mvpa.imp_map{xph}, xname);
    
    %*************** loading 1st level importance map
    fprintf('(+) loading 1st level importance map: %s\n', args.subject_id);
    
    load(fname);% subj
    
    %% ============= SAVE STANDARDIZED PATTERNS
    %%************** grp_norm_pattern{xsub}{xreg}{pos/neg}
    n_ori_pats = length(ph4.subj.patterns);
    n_iters    = length(ph5.results.iterations);
    n_patterns = n_ori_pats;% + (n_iters * n_regs);
    
    for xreg = 1:n_regs
        for i = 1:n_imp_type%pos,neg,comb
            xunit = n_patterns + i + (n_imp_type * (xreg-1));
            grp_norm_pattern{xsub}{xreg}{i} = subj.patterns{xunit}.norm.pat_mat; %#ok<*AGROW,*NASGU>
        end
    end
end

%% ============= 2ND LEVEL SAVE
fprintf('(+) save normalized importance map pattern in group structure\n\n');

save(grp_fname, 'grp_norm_pattern','-v7.3');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%============= 2ND LEVEL IMPORTANCE MAP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create_group_norm_importance_map

impmap_group   = sprintf('grp_impmap_%s_%s_%s_%s_%s_mean_cond', ...
    args.phase_name{args.xphase}, args.level, args.mask_name, args.rest, args.imp_type);
xname          = sprintf('subj_grp_imp_map_%s_%s', args.imp_type, basename);
grp_subj_fname = sprintf('%s/%s.mat', dirs.mvpa.group.imp_map{xph}, xname);

%%============= initiate group subject structure
%*************** standard mni whole brain structure

fprintf('(+) 2nd level group mean importance map\n');

clear grp_subj

xregressor_name = {'maintain','repCat','target','global'};

grp_subj     = init_subj(sprintf('%s_impmap', args.experiment), 'group');%identifier of the subj
grp_subj     = load_spm_mask_gz(grp_subj, 'mni_mask', args.mni_brain);
xmask        = get_mat(grp_subj, 'mask', 'mni_mask');

xmask_cord   = find(xmask);
n_masked_vox = size(xmask_cord, 1);

%% ============== mean patterns
for xreg = 1:n_regs
    
    clear xcriterion xgrp_mean_impmap xpeak_grp_mean_impmap
    %*************** build in pattern structure in grp_subj
    
    xcond_name = xregressor_name{xreg};
    
    it_imp = 1:n_imp_type;
    
    %% ============= positive | negative | combined
    for i = it_imp
        clear t_pattern xpattern xvol_pat
        
        fprintf('... cond%d/%d_%s...\n', xreg, n_regs, imp_type{i});
        
        mean_impmap_name = sprintf('%s%d_%s_%s', impmap_group, xreg, xcond_name, imp_type{i});
        
        %*************** get patterns
        xpattern = zeros(n_masked_vox, length(xsub_groups));
        
        for it_sub = 1:length(xsub_groups)
            xsub = xsub_groups(it_sub);
            
            t_pattern = grp_norm_pattern{xsub}{xreg}{i}; %#ok<*AGROW,*NASGU>
            xpattern(:, it_sub) = t_pattern(xmask_cord);
        end
        
        xgrp_mean_impmap{i} = zeros(n_masked_vox, 1);
        for xvox = 1:n_masked_vox
            tpat = xpattern(xvox, :);
            if sum(tpat)~=0
                xgrp_mean_impmap{i}(xvox, 1) = mean(tpat(tpat~=0));
            end
        end
        
        %*************** new volume
        xvol_pat = zeros(size(xmask));
        xvol_pat(xmask_cord) = xgrp_mean_impmap{i};
        
        %*************** grp_subj structure
        grp_subj = init_object(grp_subj,'pattern', mean_impmap_name);
        grp_subj = set_mat(grp_subj,'pattern', mean_impmap_name, xgrp_mean_impmap{i});
        
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
    end
    
    %% ============= criterion of the results
    %*************** get criterion from pos + neg distribution
    comb_mean = [];
    for i=1:2
        comb_mean = vertcat(comb_mean, xgrp_mean_impmap{i});
    end
    
    xmean = mean(comb_mean);
    xsd   = std(comb_mean)/sqrt(length(xsub_groups));
    
    %*************** cutoff table
    sorted_Y   = sort(comb_mean, 'descend');
    xcriterion = sorted_Y(round(size(comb_mean, 1) * (args.peak_thresh/2)));
    
    %% ============= threshold the results
    for i = it_imp
        
        mean_impmap_name = sprintf('%s%d_%s_%s', impmap_group, xreg, xcond_name, imp_type{i});
        
        fprintf('... thresholding top %s: cond%d/%d_%s...\n', ...
            num2str(args.peak_thresh * 100), xreg, n_regs, imp_type{i});
        
        %%%%%%%%%%%%%%%%% INTENSITY
        clear xunit
        xunit = xgrp_mean_impmap{i} >= xcriterion;
        
        xheader = {'condition','map','mean','sd','cutoff',...
            'total_voxels','select_voxels','percent','max','min'};
        xarray{i + (n_imp_type * (xreg-1)), 1}  = xcond_name;
        xarray{i + (n_imp_type * (xreg-1)), 2}  = imp_type{i};
        xarray{i + (n_imp_type * (xreg-1)), 3}  = xmean;
        xarray{i + (n_imp_type * (xreg-1)), 4}  = xsd;
        xarray{i + (n_imp_type * (xreg-1)), 5}  = xcriterion;
        xarray{i + (n_imp_type * (xreg-1)), 6}  = size(xgrp_mean_impmap{i},1);
        xarray{i + (n_imp_type * (xreg-1)), 7}  = size(find(xunit),1);
        xarray{i + (n_imp_type * (xreg-1)), 8}  = (size(find(xunit),1)/size(xgrp_mean_impmap{i},1)) * 100;
        xarray{i + (n_imp_type * (xreg-1)), 9}  = max(xgrp_mean_impmap{i}(xunit));
        xarray{i + (n_imp_type * (xreg-1)), 10} = min(xgrp_mean_impmap{i}(xunit));
        
        %*************** new volume
        xpeak_grp_mean_impmap{i} = zeros(size(xgrp_mean_impmap{i}));
        xpeak_grp_mean_impmap{i}(xunit) = abs(xgrp_mean_impmap{i}(xunit));
        
        xvol_pat_peak = zeros(size(xmask));
        xvol_pat_peak(xmask_cord) = xpeak_grp_mean_impmap{i};
        
        grp_subj = set_objfield(grp_subj, 'pattern', mean_impmap_name, 'peak_pat_mat', xpeak_grp_mean_impmap{i}, 'ignore_absence', true);
        grp_subj = set_objfield(grp_subj, 'pattern', mean_impmap_name, 'peak_pat_vol', xvol_pat_peak, 'ignore_absence', true);
        grp_subj = set_objfield(grp_subj, 'pattern', mean_impmap_name, 'peak_top', args.peak_thresh, 'ignore_absence', true);
        
        if args.peak_thresh < 1 % top 10%
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
end

fprintf('\n\n');

%% ============= SAVE grp_subj
xname          = sprintf('subj_grp_imp_map_%s_%s', args.imp_type, basename);
grp_subj_fname = sprintf('%s/%s.mat', dirs.mvpa.group.imp_map{xph}, xname);

save(grp_subj_fname, 'grp_subj', '-v7.3');

%% *************** voxel selection table
xtable = cell2table(xarray, 'VariableNames', xheader);
writetable(xtable, sprintf('%s/voxel_select_top%s_%s.csv', ...
    dirs.mvpa.group.imp_map{xph}, num2str(args.peak_thresh * 100), basename));

%% ============= CLUSTER CORRECTION/SMOOTHING
clear tsubj xarray

xclust_vox = 10;
xsmooth    = 12;
xheader    = {'condition','map','mean','sd',...
    'total_voxels','select_voxels','percent','max','min'};

if args.peak_thresh < 1 % top 10%
    xmaps = dir(fullfile(xpeak_dir, ...
        sprintf('peak_top%s_grp_impmap_operation_%s_*', ...
        num2str(args.peak_thresh*100), args.level)));
else
    xmaps = dir(fullfile(xpeak_dir, ...
        sprintf('peak_sd%s_grp_impmap_operation_%s_*',...
        num2str(args.peak_thresh), args.level)));
end

tsubj = init_subj(args.experiment, 'clearmem_impmap');%identifier of the subj
tsubj = load_spm_mask_gz(tsubj, 'MNI152_T1_3mm', args.mni_brain);

for xmap = 1:length(xmaps)
    clear xstr ximp xcmap xmat
    xx = fullfile(xpeak_dir, xmaps(xmap).name);
    
    %% ============== CLUSTER CORRECTION
    
    fprintf('(+) cluster correction: %s\n', num2str(xclust_vox));
    
    xcoutput = fullfile(xpeak_dir, sprintf('c%s_%s', ...
        num2str(xclust_vox), xmaps(xmap).name));
    
    if exist(xcoutput, 'file'), delete(xcoutput); end
    
    system(sprintf('%s/abin/3dmerge -dxyz=1 -1clust 1 %s -prefix %s %s', ...
        dirs.home, num2str(xclust_vox), xcoutput, xx));
    
    %% %%%%%%%%%%%%%%% INTENSITY
    xcmap = sprintf('c%s_%s', num2str(xclust_vox), xmaps(xmap).name);
    tsubj = load_spm_pattern_gz(tsubj, xcmap, 'MNI152_T1_3mm', xcoutput);
    xmat  = get_mat(tsubj,'pattern', xcmap);
    
    xstr  = strsplit(xmaps(xmap).name,'_');
    ximp  = strsplit(xstr{17},'.');
    
    xarray{xmap, 1} = xstr{16};
    xarray{xmap, 2} = ximp{1};
    xarray{xmap, 3} = mean(xmat);
    xarray{xmap, 4} = std(xmat);
    xarray{xmap, 5} = size(xmat,1);
    xarray{xmap, 6} = size(find(xmat),1);
    xarray{xmap, 7} = (size(find(xmat),1)/size(xmat,1)) * 100;
    xarray{xmap, 8} = max(xmat(xmat~=0));
    xarray{xmap, 9} = min(xmat(xmat~=0));
    
    %% ============== SMOOTHING
    xsoutput = fullfile(xpeak_dir, sprintf('s%s_c%s_%s', ...
        num2str(xsmooth), num2str(xclust_vox), xmaps(xmap).name));
    
    if exist(xsoutput, 'file'), delete(xsoutput); end
    
    system(sprintf('%s/abin/3dBlurToFWHM -FWHM %s -prefix %s -input %s', ...
        dirs.home, num2str(xsmooth), xsoutput, xcoutput));
    
end

%% *************** voxel selection table
xtable = cell2table(xarray, 'VariableNames', xheader);
writetable(xtable, sprintf('%s/voxel_select_top%s_c%s_%s.csv', ...
    xpeak_dir, num2str(args.peak_thresh * 100), num2str(xclust_vox), basename));

end