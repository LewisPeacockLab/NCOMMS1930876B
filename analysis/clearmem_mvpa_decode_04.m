function[] = clearmem_mvpa_decode_04(args, dirs)
% second level analysis
% Matlab R2014a

%% ============= UNPACK ARGS.
xph               = args.xphase;
mask_name         = args.mask_name;
args.regress_type = args.test_regress_type;
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

args.grp_results_name = class_basename;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%============= save MVPA parsed output in a group structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ============= 1ST LEVEL SUBJECT DATA
if args.group_mvpa{xph}
    %% *************** collect 1st level results
    for xsub = args.filtered_subs
        %*************** setup subject & directories
        args.subject_id = subject_list(xsub).name;
        dirs            = setup_directory(dirs, args);
        
        %% ============= LOAD LOCALIZER PENALTY
        %*************** load phase 6.
        fprintf(sprintf('... loading penalty from localizer: s_%s_%s\n', num2str(xsub), args.subject_id));
        
        %*************** loc_ph1. base filename
        loc_ph1.basename = sprintf('%s_%s_zscored_%s', args.phase_name{1}, args.mask_name, args.epi_name);
        
        %*************** loc_ph2. base filename
        if strcmp(args.regress_type, 'shift')
            loc_ph2.basename = sprintf('%s_%s_%s%dtr_blk_%s',...
                loc_ph1.basename, args.level, args.regress_type, ...
                args.shift_TRs, args.rest);
        elseif strcmp(args.regress_type, 'beta')
            loc_ph2.basename = sprintf('%s_%s_%s', ...
                loc_ph1.basename, args.level, args.regress_type);
        end
        
        %*************** loc_ph3. base filenames
        if args.featVox
            loc_ph3.basename = sprintf('%s_featsel_%svox', loc_ph2.basename, num2str(args.fsVosNum));
        else
            loc_ph3.basename = sprintf('%s_featsel_thresh%s', loc_ph2.basename, num2str(args.featSelThresh));
        end
        
        %*************** loc_ph4. base filenames
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
        
        mvpaout_name = sprintf('%s/mvpaout_add_%s.mat', dirs.mvpa.parse{xph}, args.results_name);
        xparse       = load(mvpaout_name);%'mvpaout_add'
        grp_mvpaout_add{xsub} = xparse.mvpaout_add;

    end
    
    %%*************** 2ND LEVEL SAVE
    fprintf('(+) save 2nd level mvpaout\n');
    
    fname = fullfile(dirs.mvpa.group.out{xph}, sprintf('grp_mvpaout_%s.mat', basename));
    save(fname, 'grp_mvpaout','-v7.3');
    
    fname = fullfile(dirs.mvpa.group.out{xph}, sprintf('grp_mvpaout_add_%s.mat', basename));
    save(fname, 'grp_mvpaout_add','-v7.3');
    
else
    %% ============= LOAD GROUP MVPAOUT
    %% *************** setup subject & directories
    args.subject_id = subject_list(1).name;
    dirs            = setup_directory(dirs, args);
    
    %%*************** 2ND LEVEL LOAD
    fprintf('(+) load 2nd level mvpaout\n');
    
    fname = fullfile(dirs.mvpa.group.out{xph}, sprintf('grp_mvpaout_%s.mat', basename));
    load(fname);%'grp_mvpaout'
    
    fprintf('(+) load 2nd level mvpaout add\n');
    
    fname = fullfile(dirs.mvpa.group.out{xph}, sprintf('grp_mvpaout_add_%s.mat', basename));
    load(fname);%'grp_mvpaout_add'
    
end

g_fsize = dir(fname);
fprintf('file size: %s GB\n', num2str(g_fsize.bytes/(10^9)));

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%============= 2ND LEVEL TIMECOURSE OF EVIDENCE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% use MATLAB 2017a
if strcmp(args.cluster,'local')
    grp_analysis_timecourse_workingmemory(args, grp_mvpaout, dirs)
    grp_analysis_timecourse_workingmemory_sorted(args, grp_mvpaout, dirs)
    grp_analysis_timecourse_workingmemory_sorted_cate(args, grp_mvpaout, dirs)
    grp_analysis_timecourse_workingmemory_add(args, grp_mvpaout_add, dirs)
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%============= CATEGORY BASED 2ND LEVEL TIMECOURSE OF EVIDENCE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% use 2017a
if strcmp(args.cluster,'local') && strcmp(args.level,'category') && ~(args.add_mvpa{xph})
    grp_analysis_timecourse_workingmemory_category(args, grp_mvpaout, dirs)
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%============= 2ND LEVEL: PERCEPTION ACCURACY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(args.cluster,'local')
    
    %% ============= 2nd LEVEL SUBJECT DATA: AUC / CLASSIFIER
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

end