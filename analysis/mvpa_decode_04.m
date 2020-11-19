function[] = mvpa_decode_04(args, dirs)
% step 4: 2nd level analysis + timecourse

%% ============= UNPACK ARGS.
xph               = args.xphase;
mask_name         = args.mask_name;
subject_list      = args.subject_list;
xsub_groups       = args.filtered_subs;
param             = args.index{xph}.param;
n_category        = param.n_category;
n_subcategory     = param.n_subcategory;

if strcmp(args.level, 'category')
    %*************** create tables
    xcate_name = param.category_name;
elseif strcmp(args.level, 'subcategory')
    %*************** create tables
    for xcate = 1:n_category
        for xsubcate = 1:n_subcategory
            xcate_name{xsubcate + n_subcategory*(xcate-1)} = ...
                param.subcategory_name{xcate}{xsubcate};
        end
    end    
end

%*************** output basename
basename          = args.analysis_basename;
        
%% ============= SETUP FILE NAMES
%*************** base file name
ph4.basename   = sprintf('%s_sh%d_%s_fselected_%s_%s_%s_%s_zepi', ...
    args.phase_name{xph}, args.shift_TRs, args.rest, mask_name, ...
    args.featSelThresh, args.level, args.epi_name);
ph5.basename   = sprintf('%s_%s', ph4.basename, args.regress_type);
class_basename = sprintf('decoding_%s_%s', ph5.basename, args.classifier);

args.grp_results_name = class_basename;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%============= save MVPA parsed output in a group structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ============= 1ST LEVEL SUBJECT DATA
%% *************** collect 1st level results
for xsub = args.filtered_subs
    %*************** setup subject & directories
    args.subject_id = subject_list(xsub).name;
    dirs            = setup_directory(dirs, args);
    
    %% ============= LOAD LOCALIZER PENALTY
    %*************** load phase 6.
    fprintf(sprintf('... loading penalty from localizer: s_%s_%s\n', num2str(xsub), args.subject_id));
    
    %*************** loc base filename
    loc_ph1.basename   = sprintf('%s_%s_zscored_%s', args.phase_name{1}, args.mask_name, args.epi_name);
    loc_ph2.basename   = sprintf('%s_%s_%s%dtr_blk_%s',...
        loc_ph1.basename, args.level, args.regress_type, ...
        args.shift_TRs, args.rest);
    loc_ph3.basename   = sprintf('%s_featsel_thresh%s', loc_ph2.basename, num2str(args.featSelThresh));
    loc_class_basename = sprintf('classified_%s_%s', loc_ph3.basename, args.classifier);
    
    %*************** load penalty check of localizer
    loc_penalty      = load(sprintf('%s/penalty_check_%s.mat', dirs.mvpa.output{1}, loc_class_basename));%'penalty_check'
    [xacc, whichmax] = max(loc_penalty.pen_check.performance); %#ok<*NODEF>
    max_penalty      = loc_penalty.pen_check.penalty(whichmax);
    args.penalty     = max_penalty;
    
    fprintf('... max_penalty: %s: acc: %1.4f\n', num2str(max_penalty), xacc);
    
    %*************** reset xpenalty
    args.xpenalty = args.penalty;
    
    %*************** ph7. basename
    args.results_name = sprintf('%s_penalty%s', class_basename, num2str(args.xpenalty));
    
    %% ============= MVPAOUT
    %%*************** evidence
    mvpaout_name = sprintf('%s/mvpaout_%s.mat', dirs.mvpa.parse{xph}, args.results_name);
    xparse       = load(mvpaout_name);%'mvpaout'
    grp_mvpaout{xsub} = xparse.mvpaout; %#ok<*AGROW,*NASGU>
    
end

%%*************** 2ND LEVEL SAVE
fprintf('(+) save 2nd level mvpaout\n');

fname = fullfile(dirs.mvpa.group.out{xph}, sprintf('grp_mvpaout_%s.mat', basename));
save(fname, 'grp_mvpaout','-v7.3');

g_fsize = dir(fname);
fprintf('file size: %s GB\n', num2str(g_fsize.bytes/(10^9)));

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%============= 2ND LEVEL TIMECOURSE OF EVIDENCE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
grp_analysis_timecourse_workingmemory_sorted(args, grp_mvpaout, dirs)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%============= 2ND LEVEL: PERCEPTION ACCURACY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%*************** random effect
% mvpa_ROC{xph}.roc_out.rand_auc{xcate}(xsub,auc)
n_condition = length(xcate_name);
n_subs      = length(xsub_groups);
xmeas_names = {'classifier_accuracy','classifier_evidence'};
xoutput_dir = dirs.mvpa.group.out{xph};

%*************** filtered subject id num
for it = 1:length(xsub_groups)
    xsub = args.filtered_subs(it);
    sub_id_filtered(it) = str2double(subject_list(xsub).name(end-2:end));
end

xout_txt    = fopen(sprintf('%s/classifier_stats_%s_n%s.txt', xoutput_dir, basename, num2str(n_subs)), 'w+');

fprintf(xout_txt, '####################################################################\n');
fprintf(xout_txt, '#### one-sample T test with normal distribution (mean: chance level)\n');
fprintf(xout_txt, '####################################################################\n\n');

for xmeas = 1:length(xmeas_names)%1_auc, 2_auc_baseline, 3_4_classifier accuracy, 4_classifier evidence
    clear xmatrix xnames xp xstats
    
    fprintf(xout_txt, '=========================================\n');
    fprintf(xout_txt, '======== %s\n\n', xmeas_names{xmeas});
    
    %*************** variable names
    xnames{1}    = 'subject';
    xnames{2}    = 'id';
    
    for xcate = 1:n_condition
        xnames{xcate + 2} = xcate_name{xcate};
    end
    
    %*************** data matrix
    xmatrix(:,1) = xsub_groups';
    xmatrix(:,2) = sub_id_filtered';
    
    if xmeas == 1, xcsv = fullfile(dirs.mvpa.group.out{xph}, sprintf('grp_acc_matrix_%s_n%s.csv', basename, num2str(n_subs)));
    else,          xcsv = fullfile(dirs.mvpa.group.out{xph}, sprintf('grp_evi_matrix_%s_n%s.csv', basename, num2str(n_subs))); end
    
    xtable  = readtable(xcsv);
    tmatrix = table2array(xtable);
    xmatrix = horzcat(xmatrix, tmatrix);
    
    xchance = 1/n_condition;
    
    %*************** total performance
    xmatrix = horzcat(xmatrix, mean(xmatrix(:,3:end),2));
    xnames{end + 1} = 'total';
    
    %%============== one-sample ttest
    %*************** AUC baseline: 0.5
    for it_col = 3:(n_condition+3), [~, xp(it_col), ~, xstats{it_col}] = ttest(xmatrix(:,it_col), xchance); end
    
    for it_col = 3:(n_condition+3), fprintf(xout_txt, '%s\t', xnames{it_col}); end
    fprintf(xout_txt, '\n');
    for it_col = 3:(n_condition+3), fprintf(xout_txt, 'T(%s)=%4.4f ', num2str(xstats{it_col}.df), xstats{it_col}.tstat); end
    fprintf(xout_txt, '\n');
    for it_col = 3:(n_condition+3), fprintf(xout_txt, '%4.4f ', xp(it_col)); end
    fprintf(xout_txt, '\n');
    
    %*************** write tables to csv files
    xmatrix = vertcat(xmatrix, xp);
    xtable   = array2table(xmatrix, 'VariableNames', xnames);
    csv_name = sprintf('%s/%s_%s_n%s.csv', xoutput_dir, xmeas_names{xmeas}, basename, num2str(n_subs));
    writetable(xtable, csv_name,'WriteRowNames',true)
end

end