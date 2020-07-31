function[] = clearmem_rsa_decode_01_tmp(args, dirs)
%*************** find ROI with GLM analysis
% only for result: saving TSPM nii for category mask

xph = 1;

%*************** SPM / MVPA PRINCETON TOOLBOX
if strcmp(args.xcluster, 'local')
    rmpath ~/Box/github/lewpealab/mvpa/
    rmpath ~/Box/github/lewpealab/mvpa/afni_matlab/
    
    addpath /Applications/spm12/
    
elseif strcmp(args.xcluster, 'blanca')
    rmpath /work/ics/data/projects/banichlab/studies/wmem/fmri/mvpa/utaustin/cu_src/mvpa/
    rmpath /work/ics/data/projects/banichlab/studies/wmem/fmri/mvpa/utaustin/cu_src/mvpa/afni_matlab/
    
    addpath /projects/ics/software/spm/spm12_v6470/
end

spm('Defaults','fMRI'); spm_jobman('initcfg');

%% ============= SETUP DIRECTORY
runs_dir       = dirs.runs{xph};
whole_mask_dir = dirs.epi_mid; 
mask_dir       = dirs.mask;
output_dir     = dirs.rsa.roi.spm;
regress_dir    = dirs.rsa.roi.regressor;

xmask = fullfile(dirs.mask, sprintf('%s.nii.gz', args.mask_name));

if ~exist(xmask, 'file')
    %*************** set FSL environment
    setenv('FSLDIR', dirs.fsl);  % this to tell where FSL folder is
    setenv('FSLOUTPUTTYPE', 'NIFTI_GZ'); % this to tell what the output type would be
    
    %*************** nii to roi.mat
    xmni_mask = fullfile(dirs.mni_mask, sprintf('%s.nii.gz', args.mask_name));
        
    fprintf('... converting mask %s \n', args.mask_name);
    
    %*************** make_mask.sh
    % # $1: FSL_DIR
    % # $2: subject's directory
    % # $3: MNI mask
    
    system(sprintf('./make_mask_mni2sub.sh %s %s %s', ...
        dirs.fsl, dirs.subj_home, xmni_mask))                
end

%% ============= LOGGING
fprintf('(+) using scratch dir: %s\n', output_dir);
%*************** turn on diary to capture analysis output
diary off;
diary(sprintf('%s/diary.txt', output_dir));
fprintf('running code: %s at %s\n', mfilename, datestr(now, 0))
fprintf('#####################################################################\n\n');

%% ============= SETUP ARGS.
fprintf('\n(+) loading args...\n');
disp(args);

% n_volumes       = args.index{xph}.param.n_volumes;
% n_runs          = args.index{xph}.param.n_runs;

%*************** PARAMETERS
TR              = 0.46;
beta_n_slice    = 16;
beta_ref_slice  = beta_n_slice/2;
beta_thresh     = 0.8;%-1000;
 
%*************** RESULT PARAMETERS
% t threshold of uncorrected p < 0.05
spm_p_thresh    = args.spm_p_thresh;
spm_v_extent    = args.spm_v_extent;
spm_thresh_type = args.spm_thresh_type;

% explicit mask: whole brain
whole_name      = args.whole_mask;
xmask_name      = args.rsa_mask;
exp_mask        = fullfile(whole_mask_dir, sprintf('%s.%s', whole_name, args.epiext));
xmask           = fullfile(mask_dir,       sprintf('%s.%s', xmask_name, args.epiext));

%% ============= MASK CHECK
setenv('FSLDIR', dirs.fsl);  % this to tell where FSL folder is
setenv('FSLOUTPUTTYPE', 'NIFTI_GZ'); % this to tell what the output type would be

%*************** anatomical ROIs
%*************** nii to roi.mat
it_roi_nii  = sprintf('%s.gz', xmask);

if ~exist(it_roi_nii, 'file')
    
    xmni_mask = fullfile(dirs.mni_mask, sprintf('%s.nii.gz', args.mask_name));
    
    fprintf('... converting mask %s \n', sprintf('%s.nii.gz', args.mask_name));
    
    %*************** make_mask.sh
    % # $1: FSL_DIR
    % # $2: subject's directory
    % # $3: MNI mask
    
    system(sprintf('./make_mask_mni2sub.sh %s %s %s', ...
        dirs.fsl, dirs.subj_home, xmni_mask))
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%============= GLM FOR ROI: SPM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%*************** explict mask
fprintf('\n(+) checking explicit mask: %s.%s\n', args.whole_mask, args.epiext);

if ~exist(exp_mask, 'file')
    gunzip(sprintf('%s.gz', exp_mask), whole_mask_dir);
end

%*************** vvs mask
if ~exist(xmask, 'file')
    gunzip(sprintf('%s.gz', xmask), mask_dir);
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%============= GLM FOR ROI: SPM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%*************** explict mask

%% ============= RUN SPM
fprintf('\n#####################################################################\n');
fprintf('* Starting SPM GLM Beta estimates subject: %s, %s\n', ...
    args.subject_id, args.phase_name{xph});

%% ============= SETUP BATCH JOB
clear matlabbatch

if ~(args.new_roi)
    %% ============= RESULTS FOR VVS
    matlabbatch{1}.spm.stats.results.spmmat                    = {fullfile(dirs.rsa.roi.spm, 'SPM.mat')};
    matlabbatch{1}.spm.stats.results.conspec.titlestr          = '';
    matlabbatch{1}.spm.stats.results.conspec.contrasts         = Inf;
    matlabbatch{1}.spm.stats.results.conspec.threshdesc        = spm_thresh_type;
    matlabbatch{1}.spm.stats.results.conspec.thresh            = spm_p_thresh;
    matlabbatch{1}.spm.stats.results.conspec.extent            = spm_v_extent;
    matlabbatch{1}.spm.stats.results.conspec.conjunction       = 1;
    matlabbatch{1}.spm.stats.results.conspec.mask.image.name   = {sprintf('%s,1', xmask)};
    matlabbatch{1}.spm.stats.results.conspec.mask.image.mtype  = 0;
    matlabbatch{1}.spm.stats.results.units                     = 1;
    matlabbatch{1}.spm.stats.results.export{1}.ps              = true;
    matlabbatch{1}.spm.stats.results.export{2}.binary.basename = ...
        sprintf('%s_thresh%s_ext%s_%s_mask', spm_thresh_type, num2str(spm_p_thresh), num2str(spm_v_extent), xmask_name);
    
    %% ============= RESULTS FOR WHOLE BRAIN
    matlabbatch{2}.spm.stats.results.spmmat                    = {fullfile(dirs.rsa.roi.spm, 'SPM.mat')};
    matlabbatch{2}.spm.stats.results.conspec.titlestr          = '';
    matlabbatch{2}.spm.stats.results.conspec.contrasts         = Inf;
    matlabbatch{2}.spm.stats.results.conspec.threshdesc        = spm_thresh_type;
    matlabbatch{2}.spm.stats.results.conspec.thresh            = spm_p_thresh;
    matlabbatch{2}.spm.stats.results.conspec.extent            = spm_v_extent;
    matlabbatch{2}.spm.stats.results.conspec.conjunction       = 1;
    matlabbatch{2}.spm.stats.results.conspec.mask.image.name   = {sprintf('%s,1', exp_mask)};
    matlabbatch{2}.spm.stats.results.conspec.mask.image.mtype  = 0;
    matlabbatch{2}.spm.stats.results.units                     = 1;
    matlabbatch{2}.spm.stats.results.export{1}.ps              = true;
    matlabbatch{2}.spm.stats.results.export{2}.binary.basename = ...
        sprintf('%s_thresh%s_ext%s_%s_mask',spm_thresh_type,num2str(spm_p_thresh),num2str(spm_v_extent), whole_name);
else
    %% ============= RESULTS FOR HIPPOCAMPUS
    matlabbatch{1}.spm.stats.results.spmmat                    = {fullfile(dirs.rsa.roi.spm, 'SPM.mat')};
    matlabbatch{1}.spm.stats.results.conspec.titlestr          = '';
    matlabbatch{1}.spm.stats.results.conspec.contrasts         = Inf;
    matlabbatch{1}.spm.stats.results.conspec.threshdesc        = spm_thresh_type;
    matlabbatch{1}.spm.stats.results.conspec.thresh            = spm_p_thresh;
    matlabbatch{1}.spm.stats.results.conspec.extent            = spm_v_extent;
    matlabbatch{1}.spm.stats.results.conspec.conjunction       = 1;
    matlabbatch{1}.spm.stats.results.conspec.mask.image.name   = {sprintf('%s,1', xmask)};
    matlabbatch{1}.spm.stats.results.conspec.mask.image.mtype  = 0;
    matlabbatch{1}.spm.stats.results.units                     = 1;
    matlabbatch{1}.spm.stats.results.export{1}.ps              = true;
    matlabbatch{1}.spm.stats.results.export{2}.binary.basename = ...
        sprintf('%s_thresh%s_ext%s_%s_mask', spm_thresh_type, num2str(spm_p_thresh), num2str(spm_v_extent), xmask_name);
    
end

%% ============= RUN JOB BATCH
tic;
fprintf('\n(+) run job batch at %s:\n', datestr(now, 0));
spm_jobman('run', matlabbatch);

%% ============= ORGANIZE DIRECTORIES
fprintf('Subject: %s end at %s\n\n', args.subject_id, datestr(now, 0));

t = toc; %#ok<*SAGROW,*AGROW,*NASGU>

fprintf('Subject: %s end at %s\n', args.subject_id, datestr(now, 0));
fprintf('... took %4.4f seconds\n', t);

%% ============= DELETE NII FILES
%*************** explict mask
fprintf('\n(+) deleting mask.nii: %s.%s\n', args.whole_mask, args.epiext);

if exist(exp_mask, 'file'), delete(exp_mask); end

%*************** vvs mask
if exist(xmask, 'file'), delete(xmask); end

%% ============= CD TO HOME
cd(args.analysis)

