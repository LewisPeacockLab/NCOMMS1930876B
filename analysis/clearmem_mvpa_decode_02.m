function[] = clearmem_mvpa_decode_02(args, dirs)
% train localizer + test study classification
% classification based on localizer penalty

%% ============= UNPACK ARGS.
args.xphase       = 2;
xph               = args.xphase;
mask_name         = args.mask_name;
args.regress_type = args.test_regress_type;

%*************** SETUP ARGS.
if args.featVox
    mask_selec = sprintf('%s_%svox', sprintf('%s_patterns_z_stat', ...
        args.train_phase), num2str(args.fsVosNum));%for ANOVA-ed mask
else
    mask_selec = sprintf('%s_thresh%s', sprintf('%s_patterns_z', ...
        args.train_phase), num2str(args.featSelThresh));%for ANOVA-ed mask
end

%*************** reset mask name
if args.class_selecting
    mask_selec = sprintf('cate%s_%s', sprintf('%d',args.selected_category), mask_selec);
end

%*************** args.mask_name
if strcmp(args.train_regress_type,'shift')
    args.mask_name = sprintf('%s_%s_%s_sh%d_blk_%s', mask_name, ...
        mask_selec, args.level, args.shift_TRs, args.rest);
elseif strcmp(args.train_regress_type,'beta')
    args.mask_name = sprintf('%s_%s_%s', mask_name, ...
        mask_selec, args.level);
end

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

%% ============= LOAD SUBJ. STRUCTURES (from classification_decode_01)
%*************** load phase 6.
fprintf('\n#####################################################################\n');
fprintf(sprintf('* loading ph6: s_%s_%s\n', num2str(args.subject_num), args.subject_id));

fname = sprintf('%s/ph6_%s.mat', dirs.mvpa.scratch{xph}, ph6.basename);
load(fname)%ph6

%*************** merge subj. structure
subj = ph6.subj; nm = ph6.nm; 
summarize(subj)

%% ============= LOAD PENALTY FROM LOCALIZER
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

%% ============= 07: MVPA CLASSIFICATION
%*************** train a classifier + do cross-validation.
fprintf('* STEP_07: decoding_classification\n');

fprintf('\n(+) subject:%s classifier: max penalty %s\n', ...
    args.subject_id, args.classifier, num2str(args.xpenalty));
fprintf('... decoding %s phase using %s phase\n', args.test_phase, args.train_phase);

[ph7.args, ph7.subj, ph7.results] = mvpa_ph04_classification(args, subj, nm);

%*************** basename phase 7.
ph7.basename  = sprintf('%s_penalty%s', class_basename, num2str(args.xpenalty));

%*************** save phase 7.
fname         = sprintf('%s/%s.mat', dirs.mvpa.output{xph}, ph7.basename);
save(fname,'ph7','-v7.3');

fprintf('\n... saved results_penalty%s of %s: %s\n', ...
    num2str(args.xpenalty), args.subject_id, fname);

end%function