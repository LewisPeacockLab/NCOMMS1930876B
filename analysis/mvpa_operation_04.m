function[] = mvpa_operation_04(args, dirs)
% step 4: 2nd level analysis

%% ============= UNPACK ARGS.
xph          = args.xphase;
subject_list = args.subject_list;

%*************** output basename
basename    = args.analysis_basename;
xsub_groups = args.g_sub;
        
%% ============= SETUP FILE NAMES
%*************** base filename
ph1.basename = sprintf('%s_%s_%s', args.phase_name{xph}, args.mask_name, args.epi_name); 
ph2.basename = sprintf('%s_%s_tr%s_blk_%s',...
    ph1.basename, args.regress_type, ...
    num2str(args.shift_TRs), args.rest);
ph3.basename = sprintf('%s_featsel_thresh%s', ph2.basename, num2str(args.featSelThresh));
ph4.basename = sprintf('%s_decoding_setup', ph3.basename);

%*************** basename phase 5.
class_basename = sprintf('classified_%s_%s', ph4.basename, args.classifier);
args.grp_results_name = class_basename;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%============= 2ND LEVEL CLASSIFIED RESULTS STRUCTURE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for xsub = xsub_groups
    %*************** setup subject & directories
    args.subject_id = subject_list(xsub).name;
    dirs            = setup_directory(dirs, args);
    
    fprintf('(+) 2nd level mvpaout: s%s_%s\n', num2str(xsub), args.subject_id);
    
    %% ============= DEFINE PENALTY
    if ~(args.fixed_penalty{xph})
        penalty_basename = sprintf('classified_%s_%s', ph3.basename, args.classifier);
        load(sprintf('%s/penalty_check_%s.mat', dirs.mvpa.output{xph}, penalty_basename));%#ok<*LOAD> %'penalty_check'
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