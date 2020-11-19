function[] = clearmem_localizer_classification_03(args, dirs)
% step 3: parse the mvpa outcome 

%% ============= UNPACK ARGS.
xph = args.xphase;

%% ============= LOGGING
diary off;
fprintf('(+) using scratch dir: %s\n', dirs.mvpa.scratch{xph});
%*************** turn on diary to capture analysis output
diary(sprintf('%s/diary.txt', dirs.mvpa.scratch{xph}));
fprintf('running code: %s\n', mfilename)
fprintf('#####################################################################\n\n');
disp(args);
fprintf('#####################################################################\n');

%% ============= Load classification_01 file
args.regress_type = args.train_regress_type;

%*************** ph1. base filename
ph1.basename = sprintf('%s_%s_zscored_%s', args.phase_name{xph}, args.mask_name, args.epi_name); 
ph2.basename = sprintf('%s_%s_%s%dtr_blk_%s',...
    ph1.basename, args.level, args.regress_type, ...
    args.shift_TRs, args.rest);
ph3.basename = sprintf('%s_featsel_thresh%s', ph2.basename, num2str(args.featSelThresh));

%*************** ph4. base filenames
class_basename  = sprintf('classified_%s_%s', ph3.basename, args.classifier);

%% ============= load penalty_check
if ~(args.fixed_penalty{xph})
    penalty_rest = args.rest;%'rest'; %
    
    %*************** base filename
    penalty_ph2.basename = sprintf('%s_%s_%s%dtr_blk_%s',...
        ph1.basename, args.level, args.regress_type, ...
        args.shift_TRs, penalty_rest);
    penalty_ph3.basename = sprintf('%s_featsel_thresh%s', penalty_ph2.basename, num2str(args.featSelThresh));
    penalty_basename = sprintf('classified_%s_%s', penalty_ph3.basename, args.classifier);
    
    %*************** load penalty / definde penalty + classification
    load(sprintf('%s/penalty_check_%s.mat', dirs.mvpa.output{xph}, penalty_basename));%'pen_check'
    
    [xacc, whichmax] = max(pen_check.performance); %#ok<*NODEF>
    max_penalty      = pen_check.penalty(whichmax);
    args.xpenalty    = max_penalty;
    
    fprintf('... max_penalty: %s: acc: %1.4f\n', num2str(max_penalty), xacc);
    
end

%*************** basename phase 4.
ph4.basename = sprintf('%s_penalty%s', class_basename, num2str(args.xpenalty));

%*************** load phase 4.
fname = sprintf('%s/%s.mat', dirs.mvpa.output{xph}, ph4.basename);
load(fname);%'ph4'

fprintf('\n... loaded classification results (penalty: %s) of %s: %s\n', ...
    num2str(args.xpenalty), args.subject_id, fname);

%% ============= 05: PARSE THE RESULTS
fprintf('\n(+) parse the results\n');

[mvpaout] = parse_mvpa_results_localizer(args, ph4, dirs); %#ok<*NASGU>

fname     = sprintf('%s/mvpaout_%s.mat', dirs.mvpa.parse{xph}, ph4.basename);
save(fname, 'mvpaout', '-v7.3');

end