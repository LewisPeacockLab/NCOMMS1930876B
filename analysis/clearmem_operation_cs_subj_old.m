function[] = clearmem_operation_cs_subj(args, dirs)

xph = args.xphase;

%*************** regressors
args.index{xph} = create_design_index_study(args, dirs);

if args.four_oper_regress
    args.regs{xph} = create_mvpa_regressor_operation_four(args, dirs);
else
    args.regs{xph} = create_mvpa_regressor_operation(args, dirs);
end

xregs = args.regs{xph}.regressors;
xsels = args.regs{xph}.selectors;

%% ============= SETUP NAMES
nm.mask       = args.mask_name;
nm.pat        = sprintf('%s_mni_patterns', args.phase_name{xph});
nm.pat_z      = sprintf('%s_z', nm.pat);% subj.patterns{2}.name;
nm.reg        = sprintf('%s_regs', args.phase_name{xph});
nm.reg_sh     = sprintf('%s_sh%d', nm.reg, args.shift_TRs);% subj.regressors{2}.name
nm.sel        = sprintf('%s_runs', args.phase_name{xph});% selectors object name
nm.sel_norest = sprintf('%s_norest', nm.reg_sh);% norest selector

%% *************** EXTRACTING PATTERNS
subj = init_subj(args.experiment, args.subject_id);%identifier of the subj

%% ============= IN MASK
if args.wholebrain, mask_dir = dirs.mni_mask; end

%*************** MASK
fprintf('\n(+) load mask: %s.%s\n', nm.mask, args.epiext);

xmask_gz = fullfile(mask_dir, sprintf('%s.%s.gz', args.mask_name, args.epiext));
subj     = load_spm_mask_gz(subj, nm.mask, xmask_gz);

%% ============= EPI PATTERN: subj.patterns{1}
fprintf('\n(+) load epi data under mask %s\n', args.mask_name);

%*************** epi
for xrun = 1:length(dirs.runs{xph})
    args.epi_names{xph}{xrun} = ...
        fullfile(dirs.runs{xph}{xrun}, sprintf('%s.%s.gz', args.epi_name, args.epiext));
end

fprintf('\n... loading epi patterns\n');

subj = load_spm_pattern_gz(subj, nm.pat, nm.mask, args.epi_names{xph});

%% ============= ZSCORING: subj.patterns{2}
%*************** zsoring each voxel, include 'rest' in the zscoring
%*************** the standard deviation of the timecourse is one.
fprintf('\n(+) z-scoring data: voxel-wise\n');

%*************** SELECTORS: subj.selectors{1}
subj = initset_object(subj,'selector', nm.sel, xsels);
subj = zscore_runs(subj, nm.pat, nm.sel);

%*************** remove original pattern: we will use z-scored patterns
subj = remove_object(subj, 'pattern', 'operation_mni_patterns');

%*************** regressor + shift
subj = initset_object(subj,'regressors', nm.reg, xregs);
subj = shift_regressors(subj, nm.reg, nm.sel, args.shift_TRs);% shift the regressor

summarize(subj)

%% ============= RESET SUBJ
% ph1.basename = sprintf('%s_%s_%s', args.phase_name{xph}, args.mask_name, args.epi_name); 
% fname = sprintf('%s/ph1_%s_test.mat', dirs.mvpa.scratch{xph}, ph1.basename);
% if ~isdir(dirs.mvpa.scratch{xph}), mkdir(dirs.mvpa.scratch{xph}); end
% load(fname);%,'subj'

%*************** selector + remove rest time points
subj = create_norest_sel(subj, nm.reg_sh);% creating norest selector object

%*************** cut pats only using timepoints
xpat_sh        = get_mat(subj, 'pattern', nm.pat_z);
xsel_norest    = get_mat(subj, 'selector', nm.sel_norest);
xunit          = find(xsel_norest);

xpat_sh_norest = xpat_sh(:, xunit);
subj           = remove_mat(subj,'pattern', nm.pat_z);
subj           = set_mat(subj,'pattern', nm.pat_z, xpat_sh_norest);

%*************** cut regs
xreg_sh        = get_mat(subj, 'regressors', nm.reg_sh);
xreg_sh_norest = xreg_sh(:, xunit);
subj           = remove_mat(subj,'regressors', nm.reg_sh);
subj           = set_mat(subj,'regressors', nm.reg_sh, xreg_sh_norest);

%*************** cut sels
xsel_sh_norest = xsel_norest(:, xunit);
subj           = set_mat(subj,'selector', nm.sel_norest, xsel_sh_norest);

%*************** removing old objects
subj = remove_object(subj, 'regressors', nm.reg);
subj = remove_object(subj, 'selector', nm.sel);
summarize(subj)

%% ============= save ph1
ph1.basename = sprintf('%s_%s_%s', args.phase_name{xph}, args.mask_name, args.epi_name); 
fname = sprintf('%s/ph1_%s.mat', dirs.mvpa.scratch{xph}, ph1.basename);
if ~isdir(dirs.mvpa.scratch{xph}), mkdir(dirs.mvpa.scratch{xph}); end

ph1.args = args; 
ph1.subj = subj;
ph1.nm   = nm;

save(fname,'ph1','-v7.3');

%% *************** write nii 
% ph1.basename = sprintf('%s_%s_%s', args.phase_name{xph}, args.mask_name, args.epi_name); 
% fname = sprintf('%s/ph1_%s.mat', dirs.mvpa.scratch{xph}, ph1.basename);
% if ~isdir(dirs.mvpa.scratch{xph}), mkdir(dirs.mvpa.scratch{xph}); end
% load(fname)
% subj = ph1.subj;

xnii_name = sprintf('zpat_operation_norest_sh%s_%s_%s_%s', ...
    num2str(args.shift_TRs), args.epi_name, nm.mask, args.subject_id);
write_to_spm(subj, 'mask', nm.mask, ...
    'output_filename', fullfile(dirs.mvpa.cs.subj, xnii_name));

fprintf('... zipping mask.nii.gz\n');
gzip(fullfile(dirs.mvpa.cs.subj, sprintf('%s.nii', xnii_name)), dirs.bold);

fprintf('... deleting mask.nii\n');
delete(fullfile(dirs.mvpa.cs.subj, sprintf('%s.nii', xnii_name)));

end