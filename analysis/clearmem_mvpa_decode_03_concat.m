function[] = clearmem_mvpa_decode_03_concat(args, dirs)
% parse the results

%% ============= UNPACK ARGS.
args.xphase       = 2;
xph               = args.xphase;
mask_name         = args.mask_name;
args.regress_type = args.test_regress_type;

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

for xcate = 1:3
    %*************** reset ph2. filename
    ph2_loc_basename = sprintf('cate%s_%s', sprintf('%d', xcate), ph2_loc.basename);
    
    %*************** ph3. base filenames
    if args.featVox
        ph3_loc.basename = sprintf('%s_featsel_%svox', ph2_loc_basename, num2str(args.fsVosNum));
    else
        ph3_loc.basename = sprintf('%s_featsel_thresh%s', ph2_loc_basename, num2str(args.featSelThresh));
    end
    
    %*************** ph4. base filenames
    loc_class_basename  = sprintf('classified_%s_%s', ph3_loc.basename, args.classifier);
    
    %*************** load penalty_check
    load(sprintf('%s/penalty_check_%s.mat', dirs.mvpa.output{1}, loc_class_basename));%'penalty_check'
    
    max_penalty     = penalty_check.penalty(penalty_check.performance==max(penalty_check.performance)); %#ok<*NODEF>
    xpenalty(xcate) = min(max_penalty);
    
    %*************** reset xpenalty
    fprintf('...... subcategory %d: max penalty: %s\n', xcate, num2str(xpenalty(xcate)))
end

%% ============= LOAD THE RESULTS
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

%*************** concatenated decoding_results
cat_result_fname = sprintf('%s/decoding_concat_%s_%s_%s.mat', dirs.mvpa.output{xph}, ...
    ph4.basename, args.regress_type, args.classifier);

if args.load_mvpa_cat{xph}
    for xcate = 1:3
        %*************** result filename
        result_fname = sprintf('decoding_cate%s_%s_%s_%s_penalty%s.mat', ...
             sprintf('%d', xcate), ph4.basename, ...
             args.regress_type, args.classifier, num2str(xpenalty(xcate)));
         
        %*************** load result
        fprintf('\n... loading classification results (penalty: %s) of %s: %s\n', ...
            num2str(xpenalty(xcate)), args.subject_id, result_fname);
        t_result = load(sprintf('%s/%s', dirs.mvpa.output{xph}, result_fname));%#ok<*AGROW> %'ph7'
        
        xresult{xcate}.results                  = t_result.ph7.results;
        xresult{xcate}.masks.n_total_voxels     = t_result.ph7.subj.masks{1}.nvox;
        xresult{xcate}.masks.n_selected_voxels  = t_result.ph7.subj.masks{3}.nvox;
        xresult{xcate}.masks.n_selected_percent = (t_result.ph7.subj.masks{3}.nvox/t_result.ph7.subj.masks{1}.nvox)*100;
        
    end
    
    fprintf('\n(+) save the concatenated results\n');
    save(cat_result_fname, 'xresult', '-v7.3');
else
    fprintf('\n(+) load the concatenated results\n');
    load(cat_result_fname)% xresult
end

%% ============= 08: PARSE THE RESULTS
%*************** concatenated parsed results
cat_mvpaout_name = sprintf('%s/mvpaout_decoding_concat_%s_%s_%s.mat', dirs.mvpa.parse{xph}, ...
    ph4.basename, args.regress_type, args.classifier);

if args.parse_mvpa{xph}    
    %*************** parse the results
    [mvpaout]  = parse_mvpa_results_study_concat(args, xresult, dirs);
    save(cat_mvpaout_name, 'mvpaout', '-v7.3');
else
    fprintf('\n(+) load the parsed results\n');
    load(cat_mvpaout_name)%'mvpaout'
end

%% ============= TIMECOURSE OF EVIDENCE

analysis_timecourse_workingmemory_concat(args, mvpaout, dirs)

end%function