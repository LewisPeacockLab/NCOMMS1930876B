function[] = clearmem_mvpa_operation_04(args, dirs)
% second level analysis

%% ============= UNPACK ARGS.
xph               = args.xphase;
args.regress_type = args.test_regress_type;
subject_list      = args.subject_list;

%*************** output basename
basename = args.analysis_basename;

if args.fixed_penalty{xph}
    xsub_groups = args.filtered_subs;
else
    xsub_groups = args.g_sub;
end
        
%% ============= SETUP FILE NAMES
%*************** ph1. base filename
ph1.basename = sprintf('%s_%s_%s', args.phase_name{xph}, args.mask_name, args.epi_name); 

%*************** ph2. base filename
if strcmp(args.regress_type, 'shift')
    ph2.basename = sprintf('%s_%s_tr%s_blk_%s',...
                   ph1.basename, args.regress_type, ...
                   num2str(args.shift_TRs), args.rest);
elseif strcmp(args.regress_type, 'beta')
    ph2.basename = sprintf('%s_%s', ...
                   ph1.basename, args.regress_type);
end

%*************** ph3. base filename
if args.featVox
    ph3.basename = sprintf('%s_featsel_%svox', ph2.basename, num2str(args.fsVosNum));
else
    ph3.basename = sprintf('%s_featsel_thresh%s', ph2.basename, num2str(args.featSelThresh));
end

%*************** ph4. base filename
ph4.basename = sprintf('%s_decoding_setup', ph3.basename);

%*************** basename phase 5.
class_basename   = sprintf('classified_%s_%s', ph4.basename, args.classifier);
args.grp_results_name = class_basename;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%============= 2ND LEVEL PARSED MVPAOUT STRUCTURE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if args.group_mvpa{xph}
    
    %% *************** 2ND LEVEL LOAD
    fname = fullfile(dirs.mvpa.group.out{xph}, sprintf('grp_mvpaout_%s.mat', basename));
    if exist(fname, 'file'), load(fname); end %'grp_mvpaout'
    
    %% ============= 1ST LEVEL
    for xsub = xsub_groups
        %*************** setup subject & directories
        args.subject_id = subject_list(xsub).name;
        dirs            = setup_directory(dirs, args);
        
        fprintf('(+) 2nd level parsed mvpaout: s%s_%s\n', num2str(xsub), args.subject_id);
        
        %% ============= DEFINE PENALTY
        if ~(args.fixed_penalty{xph})
            penalty_basename = sprintf('classified_%s_%s', ph3.basename, args.classifier);
            load(sprintf('%s/penalty_check_%s.mat', dirs.mvpa.output{xph}, penalty_basename));%'penalty_check'
            [xacc, whichmax] = max(pen_check.performance); %#ok<*NODEF>
            max_penalty      = pen_check.penalty(whichmax);
            args.xpenalty    = max_penalty;
            
            fprintf('... max_penalty: %s: acc: %1.4f\n', num2str(args.xpenalty), xacc);
        end

        %*************** basename phase 5.
        args.results_name  = sprintf('%s_penalty%s', class_basename, num2str(args.xpenalty));
        
        %% ============= MVPAOUT
        %%*************** evidence
        mvpaout_name = sprintf('%s/mvpaout_%s.mat', dirs.mvpa.parse{xph}, args.results_name);
        xparse       = load(mvpaout_name);%'mvpaout'
        
        grp_mvpaout{xsub} = xparse.mvpaout; %#ok<*AGROW,*NASGU>

    end
    
    %%*************** 2ND LEVEL SAVE
    fname = fullfile(dirs.mvpa.group.out{xph}, sprintf('grp_mvpaout_%s.mat', basename));
    save(fname, 'grp_mvpaout', '-v7.3');
    
else
    %% ============= LOAD GROUP MVPAOUT
    %% *************** setup subject & directories
    fprintf('(+) load 2nd level mvpaout\n');
    
    args.subject_id = subject_list(1).name;
    dirs            = setup_directory(dirs, args);
    
    %%*************** 2ND LEVEL LOAD
    fname = fullfile(dirs.mvpa.group.out{xph}, sprintf('grp_mvpaout_%s.mat', basename));
    load(fname);%'grp_mvpaout'
    
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%============= 2ND LEVEL CLASSIFIED RESULTS STRUCTURE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if args.mvpa_results{xph}
    
    %%*************** 2ND LEVEL LOAD
    fname = fullfile(dirs.mvpa.group.mvpa_result{xph}, sprintf('mvpa_results_%s.mat', basename));
    if exist(fname, 'file'), load(fname); end %'mvpa_results'
    
    %% ============= 1ST LEVEL
    for xsub = xsub_groups
        %*************** setup subject & directories
        args.subject_id = subject_list(xsub).name;
        dirs            = setup_directory(dirs, args);
        
        fprintf('(+) 2nd level mvpaout: s%s_%s\n', num2str(xsub), args.subject_id);
        
        %% ============= DEFINE PENALTY
        if ~(args.fixed_penalty{xph})
            penalty_basename = sprintf('classified_%s_%s', ph3.basename, args.classifier);
            load(sprintf('%s/penalty_check_%s.mat', dirs.mvpa.output{xph}, penalty_basename));%'penalty_check'
            xpen = pen_check;
            
            [xacc, whichmax] = max(xpen.performance);
            max_penalty      = xpen.penalty(whichmax);
            fprintf('... max_penalty: %s: acc: %1.4f\n', num2str(max_penalty), xacc);
            
            args.xpenalty    = max_penalty;
        end
        
        %*************** basename phase 5.
        args.results_name  = sprintf('%s_penalty%s', class_basename, num2str(args.xpenalty));
        
        %% ============= MVPAOUT
        %%*************** evidence
        reuslts_fname = sprintf('%s/%s.mat', dirs.mvpa.output{xph}, args.results_name);
        load(reuslts_fname);%'ph5'
        
        mvpa_results{xsub}.results = ph5.results; %#ok<*AGROW,*NASGU>

    end
    
    %%*************** 2ND LEVEL SAVE
    fname = fullfile(dirs.mvpa.group.mvpa_result{xph}, sprintf('mvpa_results_%s.mat', basename));
    save(fname, 'mvpa_results', '-v7.3');
    
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%============= 2ND LEVEL TIMECOURSE OF EVIDENCE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAtlab2017
if strcmp(args.cluster,'local')
    grp_analysis_timecourse_operation(args, grp_mvpaout, dirs)
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%============= 2ND LEVEL TIMECOURSE OF ZSCORED EVIDENCE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAtlab2017
if strcmp(args.cluster,'local')
    grp_analysis_timecourse_operation_zscored(args, grp_mvpaout, dirs)
end
end