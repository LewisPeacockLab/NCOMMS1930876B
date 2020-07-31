function[] = clearmem_operation_cs_01(args, dirs)

xph = args.xphase;

%% ============= INITIATE THE SUPER SUBJEC
fprintf('#############################################\n')
fprintf('[1]: concatenate data + feature selection\n')
fprintf('#############################################\n\n')

subj = init_subj(args.experiment, 'cross_subjects');%identifier of the subj

%% ============= LOAD 1ST LEVEL
xbasename = sprintf('%s_%s_%s', args.phase_name{xph}, args.mask_name, args.epi_name);

for xit = 1:length(args.filtered_subs)
    xsub = args.filtered_subs(xit);
    
    clear ph1 fname
    %*************** setup subject & directories
    args.subject_num = xsub;
    args.subject_id  = args.subject_list(xsub).name;
    sub_dirs         = dirs;
    sub_dirs         = setup_directory(sub_dirs, args);
    
    fprintf('... loading ph1: %s\n', args.subject_id);
    
    %*************** loading ph1
    fname = sprintf('%s/ph1_%s.mat', sub_dirs.mvpa.scratch{xph}, xbasename);
    load(fname);%'ph1' ph1.subj
    
    xreg_sh     = get_mat(ph1.subj, 'regressors', ph1.nm.reg_sh);
    xsel_norest = get_mat(ph1.subj, 'selector', ph1.nm.sel_norest);
    
    %*************** concatenate masked z_patterns
    if xit==1
        xregs = xreg_sh;
        xsels = xsel_norest;
        xsubj_sel = ones(1, size(xregs, 2)) * xit;
        
        xwhole_pat_z = get_masked_pattern(ph1.subj, ph1.nm.pat_z, args.mask_name);
    else
        xregs = horzcat(xregs, xreg_sh); %#ok<*AGROW>
        xsels = horzcat(xsels, xsel_norest);
        xsubj_sel = horzcat(xsubj_sel, ones(1, size(xregs, 2)) * xit);
        
        xwhole_pat_z = horzcat(xwhole_pat_z, ...
            get_masked_pattern(ph1.subj, ph1.nm.pat_z, args.mask_name));
    end
end
clear ph1

%% ============= PH1: concatenate cross-subject pattern + regressor + selector
nm.pat_z      = 'cross_subj_operation_patterns_z';
nm.reg        = 'cross_subj_operation_regs';
nm.reg_sh     = sprintf('%s_sh%d', nm.reg, args.shift_TRs);% subj.regressors{2}.name
nm.sel        = 'cross_subj_operation_selector';
nm.sel_norest = sprintf('%s_norest', nm.reg_sh);% norest selector
nm.sel_xval   = sprintf('%s_norest_xval', nm.sel);

%*************** pattern
subj = init_object(subj,'pattern', nm.pat_z);
subj = set_mat(subj,'pattern', nm.pat_z, xwhole_pat_z);
subj.patterns{1}.masked_by = args.mask_name;
subj.patterns{1}.header = ph1.subj.patterns{1}.header;

%*************** mask
subj.masks = ph1.subj.masks;

%*************** shifted regressor
subj = initset_object(subj,'regressors', nm.reg_sh, xregs);

%*************** selector + remove rest time points
subj = initset_object(subj,'selector', nm.sel_norest, xsubj_sel);

%*************** write nii
xnii_name = sprintf('zpat_operation_norest_sh%s_%s_%s_n%s', ...
    num2str(args.shift_TRs), args.epi_name, nm.mask, num2str(args.n_sub));
write_to_spm(subj, 'mask', nm.mask, ...
    'output_filename', fullfile(dirs.mvpa.cs.home, xnii_name));

fprintf('... zipping mask.nii.gz');
gzip(fullfile(dirs.cs.home, sprintf('%s.nii', xnii_name)), dirs.cs.home);

fprintf('... deleting mask.nii');
delete(fullfile(dirs.cs.home, sprintf('%s.nii', xnii_name)));

%*************** saving: in scratch
xfname = fullfile(dirs.mvpa.cs.scratch, ...
    sprintf('ph1_cs_operation_%s_%s_%s_n%s.mat', ...
    args.cs_type, args.epi_name, args.mask_name, num2str(args.n_sub)));

cs_ph1.subj = subj; %#ok<*STRNU>
cs_ph1.nm   = nm;

fprintf('... saving ph1 concatenate data ...\n')
save(xfname, 'cs_ph1', '-v7.3');
clear cs_ph1

%% ============= PH2: FEATURE SELECTION
%%%%%%% n-minus-one cross-validation schema for (within localizer only)
subj = create_xvalid_indices(subj, nm.sel, ...
    'actives_selname', nm.sel_norest, 'new_selstem', nm.sel_xval);

%%%%%%% use selector with cross-validation indices_tag #
subj = feature_select(subj, nm.pat_z, nm.reg, nm.sel_xval, ...
    'thresh', str2double(args.featSelThresh));

%*************** saving
xfname = fullfile(dirs.mvpa.cs.scratch, ...
    sprintf('ph2_cs_operation_%s_%s_%s_fsthres%s_n%s.mat', ...
    args.cs_type, args.epi_name, args.mask_name, ...
    args.featSelThresh, num2str(args.n_sub)));

cs_ph2.subj = subj;
cs_ph2.nm   = nm;

fprintf('... saving ph2 concatenate data ...\n')
save(xfname, 'cs_ph2', '-v7.3');

end