function[] = clearmem_mvpa_decode_03(args, dirs)
% parse the results

%% ============= UNPACK ARGS.
args.xphase       = 2;
xph               = args.xphase;
mask_name         = args.mask_name;
args.regress_type = args.test_regress_type;

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

%*************** reset ph4. filename
if args.class_selecting
    ph4.basename = sprintf('cate%s_%s', sprintf('%d',args.selected_category), ph4.basename);
end

%*************** ph5. base filename
ph5.basename = sprintf('%s_%s', ph4.basename, args.regress_type);

%*************** ph6. base filenames
ph6.basename     = sprintf('%s_decoding_setup_train+test', ph5.basename);

%*************** classifier name
class_basename  = sprintf('decoding_%s_%s', ph5.basename, args.classifier);

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

%*************** ph1. base filename
ph1_loc.basename = sprintf('%s_%s_zscored_%s', args.phase_name{1}, mask_name, args.epi_name); 

%*************** ph2. base filename
if strcmp(args.train_regress_type, 'shift')
    ph2_loc.basename = sprintf('%s_%s_%s%dtr_blk_%s',...
        ph1_loc.basename, args.level, args.train_regress_type, ...
        args.shift_TRs, args.rest);
elseif strcmp(args.train_regress_type, 'beta')
    ph2_loc.basename = sprintf('%s_%s_%s', ...
        ph1_loc.basename, args.level, args.regress_type);
end

%*************** reset ph2. filename
if args.class_selecting
    ph2_loc.basename = sprintf('cate%s_%s', sprintf('%d',args.selected_category), ph2_loc.basename);
end

%*************** ph3. base filenames
if args.featVox
    ph3_loc.basename = sprintf('%s_featsel_%svox', ph2_loc.basename, num2str(args.fsVosNum));
else
    ph3_loc.basename = sprintf('%s_featsel_thresh%s', ph2_loc.basename, num2str(args.featSelThresh));
end

%*************** ph4. base filenames
loc_class_basename  = sprintf('classified_%s_%s', ph3_loc.basename, args.classifier);

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

%*************** additional: parse the results
[mvpaout_add] = parse_mvpa_results_study_add(args, ph7); %#ok<*NASGU>

mvpa_parse_name_add = sprintf('%s/mvpaout_add_%s.mat', dirs.mvpa.parse{xph}, xbasename);

fprintf('\n... saving parsed mvpa results of %s: %s\n', ...
    args.subject_id, mvpa_parse_name_add);

save(mvpa_parse_name_add, 'mvpaout_add', '-v7.3');

%% ============= TIMECOURSE OF EVIDENCE

% analysis_timecourse_workingmemory(args, mvpaout, dirs)


end%function