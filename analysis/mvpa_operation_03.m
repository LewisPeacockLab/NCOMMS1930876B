function[] = mvpa_operation_03(args, dirs)
% step 3: parse the mvpa outcome/ AUC 1st level

%% ============= UNPACK ARGS.
xph = args.xphase;

%% ============= LOGGING
fprintf('(+) using scratch dir: %s\n', dirs.mvpa.scratch{xph});
%*************** turn on diary to capture analysis output
diary off;
diary(sprintf('%s/diary.txt', dirs.mvpa.scratch{xph}));
fprintf('running code: %s\n', mfilename)
fprintf('#####################################################################\n\n');
disp(args);
fprintf('#####################################################################\n');

%% ============= SETUP FILE NAMES
%*************** base filename
ph1.basename   = sprintf('%s_%s_%s', args.phase_name{xph}, args.mask_name, args.epi_name); 
ph2.basename   = sprintf('%s_%s_tr%s_blk_%s',...
    ph1.basename, args.regress_type, ...
    num2str(args.shift_TRs), args.rest);
ph3.basename   = sprintf('%s_featsel_thresh%s', ph2.basename, num2str(args.featSelThresh));
ph4.basename   = sprintf('%s_decoding_setup', ph3.basename);
class_basename = sprintf('classified_%s_%s', ph4.basename, args.classifier);

%% ============= load penalty_check
if ~(args.fixed_penalty{xph})
    penalty_basename = sprintf('classified_%s_%s', ph3.basename, args.classifier);
    load(sprintf('%s/penalty_check_%s.mat', dirs.mvpa.output{xph}, penalty_basename));%'penalty_check'
    xpen = pen_check;
    
    [xacc, whichmax] = max(xpen.performance);
    max_penalty      = xpen.penalty(whichmax);
    fprintf('... max_penalty: %s: acc: %1.4f\n', num2str(max_penalty), xacc);
    
    args.xpenalty    = max_penalty;
end

%% ============= 06: PARSE THE RESULTS
%*************** basename phase 4.
ph5.basename = sprintf('%s_penalty%s', class_basename, num2str(args.xpenalty));
fname        = sprintf('%s/%s.mat', dirs.mvpa.output{xph}, ph5.basename);
load(fname);%'ph5'

fprintf('\n... loaded classification results (penalty: %s) of %s: %s\n', ...
    num2str(args.xpenalty), args.subject_id, fname);

fprintf('\n(+) parse the results\n');

if args.four_oper_regress
    [mvpaout] = parse_mvpa_results_operation_four(args, ph5, dirs, 1);
else
    [mvpaout] = parse_mvpa_results_operation(args, ph5, dirs, 1); %#ok<*NASGU>
end

mvpaout_name = sprintf('%s/mvpaout_%s.mat', dirs.mvpa.parse{xph}, ph5.basename);
save(mvpaout_name, 'mvpaout', '-v7.3');

end