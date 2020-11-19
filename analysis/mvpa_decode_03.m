function[] = mvpa_decode_03(args, dirs)
% step 3: parse the mvpa outcome + 1st level analysis

%% ============= UNPACK ARGS.
args.xphase       = 2;
xph               = args.xphase;
mask_name         = args.mask_name;

%% ============= SETUP FILE NAMES
%*************** base file name
ph4.basename   = sprintf('%s_sh%d_%s_fselected_%s_%s_%s_%s_zepi', ...
    args.phase_name{xph}, args.shift_TRs, args.rest, mask_name, ...
    args.featSelThresh, args.level, args.epi_name);
ph5.basename   = sprintf('%s_%s', ph4.basename, args.regress_type);
class_basename = sprintf('decoding_%s_%s', ph5.basename, args.classifier);

%% ============= LOGGING
fprintf('(+) using scratch dir: %s\n', dirs.mvpa.scratch{xph});
%*************** turn on diary to capture analysis output
diary off;
diary(sprintf('%s/diary.txt', dirs.mvpa.scratch{xph}));
fprintf('running code: %s\n', mfilename)
fprintf('#####################################################################\n\n');
disp(args);
fprintf('#####################################################################\n');

%% ============= 08: PARSE THE RESULTS
fprintf('* STEP_08: parse the results\n');
fprintf('\n(+) subject:%s\n', args.subject_id);

%% ============= LOAD LOCALIZER PENALTY
fprintf(sprintf('... loading penalty from localizer: s_%s_%s\n', ...
    num2str(args.subject_num), args.subject_id));

%*************** loc base filename
ph1_loc.basename = sprintf('%s_%s_zscored_%s', args.phase_name{1}, mask_name, args.epi_name);
ph2_loc.basename = sprintf('%s_%s_%s%dtr_blk_%s',...
    ph1_loc.basename, args.level, args.train_regress_type, ...
    args.shift_TRs, args.rest);
ph3_loc.basename = sprintf('%s_featsel_thresh%s', ph2_loc.basename, num2str(args.featSelThresh));
loc_class_basename = sprintf('classified_%s_%s', ph3_loc.basename, args.classifier);

%*************** load penalty_check
load(sprintf('%s/penalty_check_%s.mat', dirs.mvpa.output{1}, loc_class_basename));%'pen_check'
[xacc, whichmax] = max(pen_check.performance); %#ok<*NODEF>
max_penalty      = pen_check.penalty(whichmax);
args.xpenalty    = max_penalty;

fprintf('... max_penalty: %s: acc: %1.4f\n', num2str(max_penalty), xacc);

%% ============= load results (phase 7)
%*************** ph7. basename
xbasename        = sprintf('%s_penalty%s', class_basename, num2str(args.xpenalty));
mvpa_parse_name  = sprintf('%s/mvpaout_%s.mat', dirs.mvpa.parse{xph}, xbasename);

%*************** load phase 7.
mvpa_result_name = sprintf('%s/%s.mat', dirs.mvpa.output{xph}, xbasename);

%% ============= 08: PARSE THE RESULTS
fprintf('\n... loading classification results (penalty: %s) of %s: %s\n', ...
    num2str(args.xpenalty), args.subject_id, mvpa_result_name);

load(mvpa_result_name);%'ph7'

fprintf('\n(+) parse the results\n');

%*************** parse the results
[mvpaout] = parse_mvpa_results_study(args, ph7, dirs); %#ok<*NASGU>

fprintf('\n... saving parsed mvpa results of %s: %s\n', ...
    args.subject_id, mvpa_parse_name);

save(mvpa_parse_name, 'mvpaout', '-v7.3');

end%function