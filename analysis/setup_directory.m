function[dirs] = setup_directory(dirs, args)

subject_id = args.subject_id;

%% ============= data directory
dirs.subj_home = fullfile(dirs.fmri_data, subject_id);
dirs.protocols = fullfile(dirs.fmri_data,'protocols');
dirs.mni_mask  = fullfile(dirs.fmri_data,'mni-masks');
dirs.ref_mask  = fullfile(dirs.fmri_data,'refer_mask');

dirs.anatomy   = fullfile(dirs.subj_home,'anatomy');
dirs.bold      = fullfile(dirs.subj_home,'bold');
dirs.epi_mid   = fullfile(dirs.bold, 'avg_func_ref');
dirs.mask      = fullfile(dirs.subj_home,'masks');
dirs.motion    = fullfile(dirs.bold, 'motion');
dirs.param     = fullfile(dirs.subj_home,'param');

run_ph = [1 2 2];

dirs.runs = {};

for xph = 1:length(run_ph)
    run_dirs = dir(fullfile(dirs.bold, sprintf('%s*%s*', ...
        args.experiment, args.phase_name{run_ph(xph)})));
    
    for i = 1:length(run_dirs)
        dirs.runs{xph}{i} = fullfile(dirs.bold, run_dirs(i).name);
    end
end

%% ============= mvpa output directory
if strcmp(args.cluster, 'local')
    xdir_fmri = dirs.fmri_data;
elseif strcmp(args.cluster, 'tacc')
    xdir_fmri = dirs.fmri_data;
elseif strcmp(args.cluster, 'blanca')
    xdir_fmri = dirs.fmri_data;
end

dirs.out.subj_home = fullfile(xdir_fmri, subject_id);
    
for xph = 1:length(args.phase_name)
    dirs.mvpa.home{xph}              = fullfile(dirs.out.subj_home, sprintf('mvpa_%s', args.phase_name{xph}));
    dirs.mvpa.group.home{xph}        = fullfile(xdir_fmri, sprintf('group_mvpa_%s', args.phase_name{xph}));
    dirs.mvpa.regressor{xph}         = fullfile(dirs.mvpa.home{xph},'regressor');
    dirs.mvpa.parse{xph}             = fullfile(dirs.mvpa.home{xph}, sprintf('parse_sh%s', num2str(args.shift_TRs)));
    dirs.mvpa.imp_map{xph}           = fullfile(dirs.mvpa.home{xph}, 'imp_map');
    dirs.mvpa.output{xph}            = fullfile(dirs.mvpa.home{xph},'output');
    dirs.mvpa.scratch{xph}           = fullfile(dirs.mvpa.home{xph},'scratch');
    dirs.mvpa.selected_mask{xph}     = fullfile(dirs.mvpa.home{xph},'selected_mask');
    
    if xph==3
        if args.stw{3}
            dirs.mvpa.home{xph}       = fullfile(dirs.out.subj_home, sprintf('mvpa_%s_stw', args.phase_name{xph}));
        end
        
        if args.four_oper_regress, n_oper = 4; else n_oper = 5; end
        
        dirs.mvpa.parse{xph}             = fullfile(dirs.mvpa.home{xph}, sprintf('parse_sh%s_%d', num2str(args.shift_TRs), n_oper));
        dirs.mvpa.imp_map{xph}           = fullfile(dirs.mvpa.home{xph}, sprintf('imp_map_%d', n_oper));
        dirs.mvpa.output{xph}            = fullfile(dirs.mvpa.home{xph}, sprintf('output_%d', n_oper));
        dirs.mvpa.scratch{xph}           = fullfile(dirs.mvpa.home{xph}, sprintf('scratch_%d', n_oper));
        
        dirs.mvpa.group.home{xph} = fullfile(xdir_fmri, sprintf('group_mvpa_%s_%d', args.phase_name{xph}, n_oper));%pass-filtered
    end
    
    if args.stw{3}
        dirs.mvpa.group.home{xph} = fullfile(xdir_fmri, sprintf('group_mvpa_%s_stw', args.phase_name{xph}));
    end
    
    dirs.mvpa.group.out{xph}         = fullfile(dirs.mvpa.group.home{xph}, sprintf('out_sh%s', num2str(args.shift_TRs)));
    dirs.mvpa.group.imp_map{xph}     = fullfile(dirs.mvpa.group.home{xph}, 'imp_map');
    dirs.mvpa.group.mvpa_result{xph} = fullfile(dirs.mvpa.group.home{xph}, 'mvpa_result');

    if xph==3
        dirs.mvpa.group.auc{xph}     = fullfile(dirs.mvpa.group.out{xph}, sprintf('auc_%s', args.rest));
        dirs.mvpa.group.auc_baseline{xph} = fullfile(dirs.mvpa.group.auc{xph}, 'baselines');
    else
        dirs.mvpa.group.auc{xph}     = fullfile(dirs.mvpa.group.out{xph}, sprintf('auc_%s_%s', args.level, args.rest));
        dirs.mvpa.group.auc_baseline{xph} = fullfile(dirs.mvpa.group.auc{xph}, 'baselines');
    end
end

%% ============= RSA output directory
if strcmp(args.xstep, 'rsa_WM')
    dirs.rsa.home          = fullfile(dirs.out.subj_home, 'rsa_study');
    dirs.rsa.roi.home      = fullfile(dirs.rsa.home, 'roi');
    dirs.rsa.roi.regressor = fullfile(dirs.rsa.roi.home, 'regressor');
    dirs.rsa.roi.spm       = fullfile(dirs.rsa.roi.home, sprintf('spm_%s', args.level));
    
    if strcmp(args.level, 'item') && (args.item_within)
        if  strcmp(args.epi_name, 'bold_mcf_brain_hpass_dt')
            dirs.rsa.roi.spm = fullfile(dirs.rsa.roi.home, 'spm_item_pf_within');
        elseif strcmp(args.epi_name, 'bold_mcf_brain')
            dirs.rsa.roi.spm = fullfile(dirs.rsa.roi.home, 'spm_item_within');
        end
    end
    
    dirs.rsa.scratch       = fullfile(dirs.rsa.home, 'scratch');
    dirs.rsa.pattern       = fullfile(dirs.rsa.home, 'pattern');
    dirs.rsa.parse         = fullfile(dirs.rsa.home, 'parse');
    
    if strcmp(args.epi_name, 'bold_mcf_brain_hpass_dt')
        dirs.rsa.group.home = fullfile(xdir_fmri, 'group_rsa_pf');
    elseif strcmp(args.epi_name, 'bold_mcf_brain')
        dirs.rsa.group.home = fullfile(xdir_fmri, 'group_rsa');
    end
    
    if args.item_within
        dirs.rsa.group.home = sprintf('%s_within', dirs.rsa.group.home);
    end
    
    dirs.rsa.group.pattern = fullfile(dirs.rsa.group.home, 'patterns');
    dirs.rsa.group.parse   = fullfile(dirs.rsa.group.home, sprintf('parse_%s', args.level));
    
    if strcmp(args.rsa_mask, 'bold_avg_mcf_brain_mask')
        dirs.rsa.group.parse = sprintf('%s_whole', dirs.rsa.group.parse);
    elseif strcmp(args.rsa_mask, 'hippocampus')
        dirs.rsa.group.parse = sprintf('%s_HC', dirs.rsa.group.parse);
    end
end

%% ============= cross-subject operation 
% dirs.mvpa.cs.home    = fullfile(xdir_fmri, sprintf('cs_mvpa_%s', args.phase_name{3}));
% dirs.mvpa.cs.scratch = fullfile(dirs.mvpa.cs.home, 'scratch');
% dirs.mvpa.cs.penalty = fullfile(dirs.mvpa.cs.scratch, 'penalty');
% dirs.mvpa.cs.out     = fullfile(dirs.mvpa.cs.home, sprintf('out_%s', args.cs_type));
% dirs.mvpa.cs.subj    = fullfile(dirs.mvpa.cs.home, subject_id); 

end