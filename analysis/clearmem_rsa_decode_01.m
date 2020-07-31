function[] = clearmem_rsa_decode_01(args, dirs)
%*************** find ROI with GLM analysis

xph = 1;

%*************** SPM / MVPA PRINCETON TOOLBOX
if strcmp(args.xcluster, 'local')
    rmpath ~/Documents/github/lewpealab/mvpa/
    rmpath ~/github/lewpealab/mvpa/afni_matlab/
    
    addpath /Applications/spm12/
    
elseif strcmp(args.xcluster, 'blanca')
    rmpath('/work/ics/data/projects/banichlab/studies/wmem/fmri/mvpa/utaustin/cu_src/mvpa/')
    
    addpath(genpath('/projects/ics/software/spm/spm12_v6470'))
    addpath('/projects/ics/software/spm/spm12_v6470/compat')
end

spm('Defaults','fMRI'); spm_jobman('initcfg');

%% ============= SETUP DIRECTORY
runs_dir       = dirs.runs{xph};
whole_mask_dir = dirs.mask;%dirs.epi_mid; 
mask_dir       = dirs.mask;
output_dir     = dirs.rsa.roi.spm;
regress_dir    = dirs.rsa.roi.regressor;

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

n_volumes = args.index{xph}.param.n_volumes;
n_runs    = args.index{xph}.param.n_runs;
xnames    = args.regs{xph}.names;
n_class   = length(xnames);

if args.item_within
    n_item_within = 18;
end

%*************** PARAMETERS
TR              = 0.46;
beta_n_slice    = 16;
beta_ref_slice  = beta_n_slice/2;
beta_thresh     = 0.8;%-1000;
if args.item_within
    beta_thresh = -1000;
end
 
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

%*************** functionals
for xrun = 1:n_runs
    xepi = fullfile(runs_dir{xrun}, ...
        sprintf('%s.%s', args.epi_name, args.epiext));
    
    fprintf('\n(+) checking functionals: %s\n', xepi);
    
    if ~exist(xepi, 'file')
        gunzip(sprintf('%s.gz',xepi), ...
            fullfile(runs_dir{xrun}));
    end
end

%% ============= RUN SPM
fprintf('\n#####################################################################\n');
fprintf('* Starting SPM GLM Beta estimates subject: %s, %s\n', ...
    args.subject_id, args.phase_name{xph});

%% ============= REGRESSORS
%*************** define file names
cond_regress  = fullfile(regress_dir, sprintf('localizer_spm_regressor_%s.mat', args.level));
multi_regress = fullfile(regress_dir, 'localizer_multi_regressor.txt');

fprintf('\n(+) defining multiple conditons: %s\n', cond_regress);
fprintf('\n(+) defining multiple regressors: %s\n', multi_regress);

%% ============= SETUP BATCH JOB
clear matlabbatch

fprintf('\n(+) setup batch job:\n');
matlabbatch{1}.spm.stats.fmri_spec.dir            = {output_dir};
matlabbatch{1}.spm.stats.fmri_spec.timing.units   = 'scans';
matlabbatch{1}.spm.stats.fmri_spec.timing.RT      = TR;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t  = beta_n_slice;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = beta_ref_slice;

%% ============= EPIS DATA
fprintf('\n(+) reading functionals:\n');

for xrun = 1:n_runs
    for xvol = 1:n_volumes(xrun)
        matlabbatch{1}.spm.stats.fmri_spec.sess.scans{xvol + sum(n_volumes(1:xrun-1)), 1} = ...
            fullfile(runs_dir{xrun}, sprintf('%s.%s,%s', args.epi_name, args.epiext, num2str(xvol)));
    end
end

%% ============= DESIGN
%*************** fmri model specification
fprintf('\n(+) fmri model specification:\n');

matlabbatch{1}.spm.stats.fmri_spec.sess.multi       = {cond_regress};
matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg   = {multi_regress};
matlabbatch{1}.spm.stats.fmri_spec.sess.hpf         = 128;
matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
matlabbatch{1}.spm.stats.fmri_spec.volt             = 1;%do not model interaction between trial types
matlabbatch{1}.spm.stats.fmri_spec.global           = 'None';%global normalization
matlabbatch{1}.spm.stats.fmri_spec.mthresh          = beta_thresh;
matlabbatch{1}.spm.stats.fmri_spec.mask             = {sprintf('%s,1',exp_mask)};%explicit mask
matlabbatch{1}.spm.stats.fmri_spec.cvi              = 'AR(1)';

%*************** model estimation
fprintf('\n(+) model estimation:\n');

matlabbatch{2}.spm.stats.fmri_est.spmmat(1)         = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals   = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical  = 1;

%% ============= CONTRAST
matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));

if strcmp(args.level, 'category')
    matlabbatch{3}.spm.stats.con.consess{1}.fcon.name    = 'face vs. fruit vs. scene';
    matlabbatch{3}.spm.stats.con.consess{1}.fcon.weights = [1 0 0; 0 1 0; 0 0 1];
    matlabbatch{3}.spm.stats.con.consess{1}.fcon.sessrep = 'none';
    matlabbatch{3}.spm.stats.con.consess{2}.tcon.name    = 'face>others';
    matlabbatch{3}.spm.stats.con.consess{2}.tcon.weights = [1 -0.5 -0.5];
    matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
    matlabbatch{3}.spm.stats.con.consess{3}.tcon.name    = 'fruit>others';
    matlabbatch{3}.spm.stats.con.consess{3}.tcon.weights = [-0.5 1 -0.5];
    matlabbatch{3}.spm.stats.con.consess{3}.tcon.sessrep = 'none';
    matlabbatch{3}.spm.stats.con.consess{4}.tcon.name    = 'scene>others';
    matlabbatch{3}.spm.stats.con.consess{4}.tcon.weights = [-0.5 -0.5 1];
    matlabbatch{3}.spm.stats.con.consess{4}.tcon.sessrep = 'none';

elseif strcmp(args.level, 'item')
    for xitem = 1:n_class
        if args.item_within
            xcate   = fix((xitem-1)/n_item_within)+1;
            xunit   = (1:n_item_within) + (n_item_within * (xcate-1));
            xweight = zeros(1, n_class);
            xweight(xunit) = -(1/(n_item_within-1));
            xweight(xitem) = 1;
        else
            xweight = ones(1, n_class) * -(1/(n_class-1));
            xweight(xitem) = 1;
        end
        
        matlabbatch{3}.spm.stats.con.consess{xitem}.tcon.name    = ...
            sprintf('%s>others', xnames{xitem});
        matlabbatch{3}.spm.stats.con.consess{xitem}.tcon.weights = xweight;
        matlabbatch{3}.spm.stats.con.consess{xitem}.tcon.sessrep = 'none';
    end    
end

matlabbatch{3}.spm.stats.con.delete = 0;

%% ============= RESULTS FOR VVS
matlabbatch{4}.spm.stats.results.spmmat                    = {fullfile(dirs.rsa.roi.spm, 'SPM.mat')};
matlabbatch{4}.spm.stats.results.conspec.titlestr          = '';
matlabbatch{4}.spm.stats.results.conspec.contrasts         = Inf;
matlabbatch{4}.spm.stats.results.conspec.threshdesc        = spm_thresh_type;
matlabbatch{4}.spm.stats.results.conspec.thresh            = spm_p_thresh;
matlabbatch{4}.spm.stats.results.conspec.extent            = spm_v_extent;
matlabbatch{4}.spm.stats.results.conspec.conjunction       = 1;
matlabbatch{4}.spm.stats.results.conspec.mask.image.name   = {sprintf('%s,1', xmask)};
matlabbatch{4}.spm.stats.results.conspec.mask.image.mtype  = 0;
matlabbatch{4}.spm.stats.results.units                     = 1;
matlabbatch{4}.spm.stats.results.export{1}.ps              = true;
matlabbatch{4}.spm.stats.results.export{2}.binary.basename = ...
    sprintf('%s_thresh%s_ext%s_%s_mask',spm_thresh_type,num2str(spm_p_thresh),num2str(spm_v_extent), xmask_name);

%% ============= RESULTS FOR WHOLE BRAIN
matlabbatch{5}.spm.stats.results.spmmat                    = {fullfile(dirs.rsa.roi.spm, 'SPM.mat')};
matlabbatch{5}.spm.stats.results.conspec.titlestr          = '';
matlabbatch{5}.spm.stats.results.conspec.contrasts         = Inf;
matlabbatch{5}.spm.stats.results.conspec.threshdesc        = spm_thresh_type;
matlabbatch{5}.spm.stats.results.conspec.thresh            = spm_p_thresh;
matlabbatch{5}.spm.stats.results.conspec.extent            = spm_v_extent;
matlabbatch{5}.spm.stats.results.conspec.conjunction       = 1;
matlabbatch{5}.spm.stats.results.conspec.mask.image.name   = {sprintf('%s,1', exp_mask)};
matlabbatch{5}.spm.stats.results.conspec.mask.image.mtype  = 0;
matlabbatch{5}.spm.stats.results.units                     = 1;
matlabbatch{5}.spm.stats.results.export{1}.ps              = true;
matlabbatch{5}.spm.stats.results.export{2}.binary.basename = ...
    sprintf('%s_thresh%s_ext%s_%s_mask',spm_thresh_type,num2str(spm_p_thresh),num2str(spm_v_extent), whole_name);

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

%*************** functionals
for xrun = 1:n_runs
    fprintf('\n(+) deleting functionals.nii: %s\n', runs_dir{xrun});
    
    xepi = fullfile(runs_dir{xrun}, ...
        sprintf('%s.%s', args.epi_name, args.epiext));
    
    if exist(xepi, 'file'), delete(xepi); end
end

%% ============= CD TO HOME
cd(args.analysis)

