function[] = create_spm_mask(args, dirs)

subject_id = args.subject_id;

%% ============= MASK CHECK
%*************** set FSL environment
setenv('FSLDIR', dirs.fsl);  % this to tell where FSL folder is
setenv('FSLOUTPUTTYPE', 'NIFTI_GZ'); % this to tell what the output type would be

for i = 1:2
    xmask = fullfile(dirs.mask, sprintf('%s.nii.gz', args.mask_name{i}));
    
    if ~exist(xmask, 'file')
        %*************** MNI to subject functional space
        fprintf('... converting MNI mask %s to subject space \n', args.mask_name{i});
        xmni_mask = fullfile(dirs.mni_mask, sprintf('%s.nii.gz', args.mask_name{i}));
        
        %*************** make_mask.sh
        % # $1: FSL_DIR
        % # $2: subject's directory
        % # $3: MNI mask
        system(sprintf('./make_mask_mni2sub.sh %s %s %s', ...
            dirs.fsl, dirs.subj_home, xmni_mask))
    end
end

%% ============= COMBINE MASK WITH COMMON VOXELS
clear subj
subj = init_subj(args.experiment, subject_id);%identifier of the subj

%*************** explicit mask (epi whole brain mask)
explicit_mask = fullfile(dirs.rsa.roi.home, 'spm_category', 'mask.nii');
if ~exist(explicit_mask,'file')
    gunzip(sprintf('%s.gz', explicit_mask), dirs.epi_mid);
end

subj = load_spm_mask(subj, 'wholebrain', explicit_mask);
whole_mask_mat = get_mat(subj, 'mask', 'wholebrain');
whole_mask_coord = find(whole_mask_mat==1);

%*************** filter all rois with wholebrain mask
for i = 1:length(args.mask_name)
    roi_mask = fullfile(dirs.mask, sprintf('%s.nii', args.mask_name{i}));
    
    if ~exist(roi_mask,'file')
        gunzip(sprintf('%s.gz', roi_mask), dirs.mask);
    end
    
    subj = load_spm_mask(subj, args.mask_name{i}, roi_mask);
    
    roi_mask_mat  = get_mat(subj, 'mask', args.mask_name{i});
    roi_mask_coord = find(roi_mask_mat==1);

    %*************** generate new mask with overlapping voxels
    xfiltered_coord = whole_mask_coord(ismember(whole_mask_coord, roi_mask_coord));
    new_mask = fullfile(dirs.mask, sprintf('%s_filtered.nii', args.mask_name{i}));
    
    refer_mask = spm_vol(explicit_mask);
    xmask      = spm_read_vols(refer_mask); xmask(:,:,:) = 0; xmask(xfiltered_coord)=1;
    refer_mask.fname = new_mask;
    spm_write_vol(refer_mask, xmask);
    
    gzip(new_mask, dirs.mask);
    
    delete(roi_mask);
    delete(new_mask);
end
end