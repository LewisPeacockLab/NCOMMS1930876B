function[] = clearmem_mvpa_decode_04_concat(args, dirs)
% second level analysis

%% ============= UNPACK ARGS.
args.xphase       = 2;
xph               = args.xphase;
mask_name         = args.mask_name;
args.regress_type = args.test_regress_type;
subject_list      = args.subject_list;

%*************** output basename
basename          = args.analysis_basename;
        
%% ============= SETUP FILE NAMES
%*************** ph4. base file name
if strcmp(args.regress_type, 'shift')
    ph4.basename = sprintf('%s_sh%d_%s_fselected_%s_%s_%s_%s_zepi', ...
        args.phase_name{xph}, args.shift_TRs, args.rest, mask_name, ...
        args.featSelThresh, args.level, args.epi_name); 
elseif strcmp(args.regress_type, 'beta')
    ph4.basename = sprintf('%s_fselected_%s_%s_%s_%s_zepi', ...
        args.phase_name{xph}, mask_name, ...
        args.featSelThresh, args.level, args.beta_name); 
end

%% ============= 1ST LEVEL SUBJECT DATA
if args.group_mvpa{xph}
    for xsub = args.g_sub
        %*************** setup subject & directories
        args.subject_id = subject_list(xsub).name;
        dirs            = setup_directory(dirs, args);
        
        %% ============= LOAD LOCALIZER PENALTY
        %*************** load phase 6.
        fprintf(sprintf('... loading penalty from localizer: s_%s_%s\n', num2str(xsub), args.subject_id));
        
        %% ============= MVPAOUT
        %*************** concatenated parsed results
        cat_mvpaout_name  = sprintf('%s/mvpaout_decoding_concat_%s_%s_%s.mat', dirs.mvpa.parse{xph}, ...
            ph4.basename, args.regress_type, args.classifier);
        xparse            = load(cat_mvpaout_name);%'mvpaout'
        
        grp_mvpaout{xsub} = xparse.mvpaout; %#ok<*AGROW,*NASGU>  

    end
    
    %%*************** 2ND LEVEL SAVE
    fprintf('(+) save 2nd level mvpaout\n');
    
    fname = fullfile(dirs.mvpa.group.out{xph}, sprintf('grp_mvpaout_%s.mat', basename));
    save(fname, 'grp_mvpaout','-v7.3');
    
else
    %% ============= LOAD GROUP MVPAOUT
    %% *************** setup subject & directories
    args.subject_id = subject_list(1).name;
    dirs            = setup_directory(dirs, args);
    
    %%*************** 2ND LEVEL LOAD
    fprintf('(+) load 2nd level mvpaout\n');
    
    fname = fullfile(dirs.mvpa.group.out{xph}, sprintf('grp_mvpaout_%s.mat', basename));
    load(fname);%'grp_mvpaout'
    
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%============= 2ND LEVEL TIMECOURSE OF EVIDENCE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

grp_analysis_timecourse_workingmemory_concat(args, grp_mvpaout, dirs)

end