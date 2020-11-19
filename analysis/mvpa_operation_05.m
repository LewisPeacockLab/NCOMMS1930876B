function[] = mvpa_operation_05(args, dirs)
% step 5: AUC/ROC 2nd

%% ============= UNPACK ARGS.
xph               = args.xphase;
args.regress_type = args.test_regress_type;
subject_list      = args.subject_list;

%*************** filtered subject id num
for it = 1:length(args.filtered_subs)
    xsub = args.filtered_subs(it);
    sub_id_filtered(it) = str2double(subject_list(xsub).name(end-2:end));
end

%*************** output basename
basename          = args.analysis_basename;
xoutput_dir       = dirs.mvpa.group.auc{xph};

if args.four_oper_regress
    it_conds = [1 2 4 5];
    cond_names = {'Maintain','RepCategory','Suppress','Clear'};
else
    it_conds = 1:5;
    cond_names = {'Maintain','RepCategory','RepSubcate','Suppress','Clear'};
end

n_target    = length(cond_names);
xsub_groups = args.g_sub;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%============= ROC for single subject + save mvpa results in group structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ============= load group mvpa results from clearmem_mvpa_operation_04
fname = sprintf('%s/mvpa_results_%s.mat', dirs.mvpa.group.mvpa_result{xph}, basename);
load(fname);%#ok<*LOAD> %, 'mvpa_results'

%% ============= SETUP
roc_matrix = []; perf_matrix = [];
roc_matrix_ex_rest = []; 

%% ============= 1ST LEVEL
for xsub = xsub_groups
    fprintf('%s roc: %s\n', args.phase_name{xph}, subject_list(xsub).name);
    
    %*************** reset desired
    args.subject_id = subject_list(xsub).name;
    dirs = setup_directory(dirs, args);
    
    clear index xmatrix
    load(fullfile(dirs.param, sprintf('%s_parameters_%s.mat', args.phase_name{2}, args.subject_id))); %index
    xmatrix  = index.matrix; xheader = index.header;
    xresults = mvpa_results{xsub}.results; %#ok<USENS>
    
    n_target = size(xresults.iterations(1).acts, 1);
    n_iter   = size(xresults.iterations, 2);
    for xiter = 1:n_iter
        n_vol(xiter) = size(xresults.iterations(xiter).acts, 2);
    end
    
    %***************** reset regs: shifting TR
    xregs.runs         = args.regs{xph}.selectors;
    xregs.operation_sh = zeros(1, sum(n_vol));% only operation window
    
    for it = 1:n_target
        xcond = it_conds(it);
        xunit = getDATA(xmatrix', xheader, {'condition'}, {xcond});
        
        xregs.timecourse(xunit) = xcond;
        
        xunit = find(getDATA(xmatrix', xheader, ...
            {'condition','manipulation'}, {xcond,1}));
        
        xregs.operation_sh(xunit + args.shift_TRs) = it;
    end
    
    %***************** remove spikes
    xspike = getDATA(xmatrix', xheader, {'spike'}, {1});
    xregs.operation_s(xspike) = 0;
    
    %*************** extract MVPA performance
    acts_array = [];
    for xiter = 1:n_iter
        acts_array = horzcat(acts_array, xresults.iterations(xiter).acts); %#ok<*AGROW>
    end
    
    xunit     = (xregs.operation_sh ~=0);
    xdesireds = xregs.operation_sh(xunit);
    xacts     = acts_array(:,xunit);
    
    %*************** ROC inputs: perfcurve(labels, scores, posclass)
    % scores: classifier predictions: (evidence)
    % labels: given true class labels: xdesireds
    % posclass: positive class label: (xcate)
    % T: (mean, lower bound, upper bound)
    for xcate = 1:n_target
        
        xlabels   = xdesireds;
        xscores   = xacts(xcate, :);
        xposclass = xcate;
        
        %*************** performance
        clear t_matrix
        
        t_matrix(:, 3:4) = [xlabels', xscores'];
        t_matrix(:, 1)   = xcate;
        t_matrix(:, 2)   = xsub;
        
        perf_matrix = vertcat(perf_matrix, t_matrix);
        
        %*************** perfcurve
        clear t_matrix
        
        [x_fa, y_hit, T, AUC] = perfcurve(xlabels, xscores, xposclass);
        
        t_matrix(:, 3:5) = [x_fa, y_hit, T];
        t_matrix(:, 1)   = xcate;
        t_matrix(:, 2)   = xsub;
        
        roc_matrix = vertcat(roc_matrix, t_matrix);
        
        mvpa_ROC{xph}.roc_out.rand_auc{xcate}(xsub, :) = [xsub, AUC];
        mvpa_ROC{xph}.baseline{xsub}{xcate} = 1/n_target;
        
    end
end

mvpa_ROC{xph}.roc_out.header = {'category','subject','x_fa','y_hit','threshold'};
mvpa_ROC{xph}.roc_out.matrix = roc_matrix;

if strcmp(args.rest, 'rest')
    mvpa_ROC{xph}.roc_out.matrix_ex_rest = roc_matrix_ex_rest;
end

%*************** save mvpa_ROC
fname = sprintf('%s/mvpa_roc_%s.mat', xoutput_dir, basename);
save(fname, 'mvpa_ROC', '-v7.3')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%============= WRITING STATS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xsub_groups = args.filtered_subs;

%%============= 2nd LEVEL SUBJECT DATA: AUC / CLASSIFIER
%*************** random effect
% mvpa_ROC{xph}.roc_out.rand_auc{xcate}(xsub,auc)
if args.four_oper_regress
    cond_names = {'Maintain','RepCategory','Suppress','Clear'};
else
    cond_names = {'Maintain','RepCategory','RepSubcate','Suppress','Clear'};
end

n_condition = length(cond_names);
n_subs      = length(xsub_groups);
xmeas_names = {'classifier_auc','classifier_auc_baseline','classifier_accuracy','classifier_evidence'};

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
        xnames{xcate + 2} = cond_names{xcate};
    end
    
    %*************** data matrix
    xmatrix(:,1) = xsub_groups';
    xmatrix(:,2) = sub_id_filtered';
    
    if (xmeas == 1) || (xmeas == 2)
        for xcate = 1:n_condition
            if xmeas == 1, xdat = mvpa_ROC{xph}.roc_out.rand_auc{xcate}(:,2);
            else,          xdat = ones(args.n_sub, 1) * 0.5;%  mvpa_ROC{xph}.baseline_criterion{xcate}(:,1);
            end
            
            xmatrix = horzcat(xmatrix, xdat(xsub_groups));
        end
        
        xchance = 0.5;%for auc
        
    elseif (xmeas == 3) || (xmeas == 4)
        if xmeas == 3, xcsv = fullfile(dirs.mvpa.group.out{xph}, sprintf('grp_acc_matrix_%s_n%s.csv', basename, num2str(n_subs)));
        else,          xcsv = fullfile(dirs.mvpa.group.out{xph}, sprintf('grp_evi_matrix_%s_n%s.csv', basename, num2str(n_subs))); end
        
        xtable  = readtable(xcsv);
        tmatrix = table2array(xtable);
        xmatrix = horzcat(xmatrix, tmatrix);
        
        xchance = 1/n_condition;
    end
    
    %*************** total performance
    xmatrix = horzcat(xmatrix, mean(xmatrix(:,3:end),2));
    xnames{end + 1} = 'total';
    
    %*************** rest classification excluding rest in AUC
    if (xmeas==1) && (strcmp(args.rest, 'rest'))
        for xcate = 1:(n_target-1)
            xmatrix = horzcat(xmatrix, mvpa_ROC{xph}.roc_out.rand_auc_ex_rest{xcate}(xsub_groups,2));
            xnames{xcate + 2 + n_condition} = sprintf('%s_ex_rest', cond_names{xcate});
        end
    end
    
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

%% ============= COMPARING: ACCURACY/EVIDENCE
clear xnames xconfusion
xmeas_names = {'accuracy','evidence'};
targ_array  = 1:n_condition;

fprintf(xout_txt, '\n\n####################################################################\n');
fprintf(xout_txt, '#### paired T test\n');
fprintf(xout_txt, '####################################################################\n\n');

for xmeas = 1:length(xmeas_names)%1_auc, 2_auc_baseline, 3_4_classifier accuracy, 4_classifier evidence
    clear xtable xheader xtarget xguess
    
    fprintf(xout_txt, '\n=========================================\n');
    fprintf(xout_txt, '======== %s\n\n', xmeas_names{xmeas});
    
    %*************** data matrix
    xcsv = fullfile(dirs.mvpa.group.out{xph}, ...
        sprintf('grp_%s_%s_n%s.csv', xmeas_names{xmeas}, basename, num2str(n_subs)));
    
    xtable  = readtable(xcsv);
    xheader = xtable.Properties.VariableNames;
    
    xtarget = unique(xtable.Row)';
    xguess  = xtable.Properties.VariableNames(~ismember(xtable.Properties.VariableNames, 'Row'));
    
    for xtarg = 1:length(xtarget)
        xunit = ismember(xtable.Row, xtarget{xtarg});
        xconfusion{xtarg} = xtable(xunit, :);
    end
    
    for xtarg = 1:length(xtarget)
        
        fprintf(xout_txt, '#### target: %s\n', cond_names{xtarg});
        
        for gg = 1:length(xtarget)
            xx = xconfusion{xtarg}(:, findCol(xheader, xguess{gg}));
            fprintf(xout_txt, 'guess %s: M=%4.4f, SE=%4.4f\n', cond_names{gg}, ...
                mean(xx.Variables), std(xx.Variables)/sqrt(length(xsub_groups)));
        end
        
        fprintf(xout_txt,'\n');
        
        clear xpvalue xstats
        it_targs = targ_array(~ismember(targ_array, xtarg));
        for it_y = 1:length(it_targs)
            ytarg = it_targs(it_y);
            
            %*************** ttest
            xx = xconfusion{xtarg}(:, findCol(xheader, xguess{xtarg}));
            yy = xconfusion{xtarg}(:, findCol(xheader, xguess{ytarg}));
            
            [~, xpvalue(it_y), ~, xstats] = ttest(xx.Variables, yy.Variables, 'Alpha', 0.05);
            
            fprintf(xout_txt, '%s vs. %s: T(%s)=%4.4f, P=%4.4f\n', ...
                cond_names{xtarg}, cond_names{ytarg}, num2str(xstats.df), xstats.tstat, xpvalue(it_y));
            
        end
        
        fprintf(xout_txt, '\n\n********* FDR p-value\n');
        
        [~, ~, ~, xadj_p] = fdr_bh(xpvalue, 0.05, 'pdep', 'no');
        
        for it_y = 1:length(it_targs)
            ytarg = it_targs(it_y);
            fprintf(xout_txt, '%s vs. %s: P=%4.4f\n', ...
                cond_names{xtarg}, cond_names{ytarg}, xadj_p(it_y));
        end
        
        fprintf(xout_txt, '\n');
        
    end
end

%% ============= COMPARING: REPLACES: ACCURACY/EVIDENCE
targ_names = {'target','nontarget'};
if ~(args.four_oper_regress)
    clear xconfusion
    xmeas_names = {'accuracy','evidence'};
    
    fprintf(xout_txt, '\n\n####################################################################\n');
    fprintf(xout_txt, '#### Replaces: paired T test\n');
    fprintf(xout_txt, '####################################################################\n\n');
    
    for xmeas = 1:length(xmeas_names)%1_auc, 2_auc_baseline, 3_4_classifier accuracy, 4_classifier evidence
        clear xtable xheader xtarget xguess
        
        fprintf(xout_txt, '\n=========================================\n');
        fprintf(xout_txt, '======== %s\n\n', xmeas_names{xmeas});
        
        %*************** data matrix
        xcsv = fullfile(dirs.mvpa.group.out{xph}, ...
            sprintf('grp_%s_%s_n%s.csv', xmeas_names{xmeas}, basename, num2str(n_subs)));
        
        xtable  = readtable(xcsv);
        xheader = xtable.Properties.VariableNames;
        
        xtarget = unique(xtable.Row)';
        xguess  = xtable.Properties.VariableNames(~ismember(xtable.Properties.VariableNames, 'Row'));
        
        for xtarg = 1:length(xtarget)
            xunit = ismember(xtable.Row, xtarget{xtarg});
            xconfusion{xtarg} = xtable(xunit, :);
        end
        
        %*************** 1_target (actual), 2_guessed (predicted)
        xx1 = xconfusion{2}(:, findCol(xheader, {xguess{2}, xguess{3}}));
        xx2 = xconfusion{3}(:, findCol(xheader, {xguess{3}, xguess{2}}));
        
        fprintf(xout_txt, '#### replaces: target vs. nontarg\n');
        
        xx = vertcat(xx1, xx2);
        
        for xtarg = 1:2
            yy = xx.Variables;
            
            fprintf(xout_txt, '%s: M=%4.4f, SE=%4.4f\n', targ_names{xtarg}, ...
                mean(yy(:,xtarg)), std(yy(:,xtarg))/sqrt(length(xsub_groups)));
        end
        
        fprintf(xout_txt,'\n');
        
        clear xpvalue xstats
        
        [~, xpvalue, ~, xstats] = ttest(yy(:,1), yy(:,2), 'Alpha', 0.05);
        
        fprintf(xout_txt, 'targ vs. nontarg: T(%s)=%4.4f, P=%4.4f\n', ...
            num2str(xstats.df), xstats.tstat, xpvalue);
        
    end
else
    %% ============= COMPARING CONDITIONS: suppress(3) vs. clear(4)
    clear xconfusion
    xmeas_names = {'accuracy','evidence'};
    
    fprintf(xout_txt, '\n\n####################################################################\n');
    fprintf(xout_txt, '#### Suppress vs. Clear: paired T test\n');
    fprintf(xout_txt, '####################################################################\n\n');
    
    for xmeas = 1:length(xmeas_names)%1_auc, 2_auc_baseline, 3_4_classifier accuracy, 4_classifier evidence
        clear xtable xheader xtarget xguess
        
        fprintf(xout_txt, '\n=========================================\n');
        fprintf(xout_txt, '======== %s\n\n', xmeas_names{xmeas});
        
        %*************** data matrix
        xcsv = fullfile(dirs.mvpa.group.out{xph}, ...
            sprintf('grp_%s_%s_n%s.csv', xmeas_names{xmeas}, basename, num2str(n_subs)));
        
        xtable  = readtable(xcsv);
        xheader = xtable.Properties.VariableNames;
        
        xtarget = unique(xtable.Row)';
        xguess  = xtable.Properties.VariableNames(~ismember(xtable.Properties.VariableNames, 'Row'));
        
        for xtarg = 1:length(xtarget)
            xunit = ismember(xtable.Row, xtarget{xtarg});
            xconfusion{xtarg} = xtable(xunit, :);
        end
        
        %*************** 1_target (actual), 2_guessed (predicted)
        % suppress(3) vs. clear(4)
        t_xx1 = xconfusion{3}(:, findCol(xheader, {xguess{3}, xguess{4}}));
        t_xx2 = xconfusion{4}(:, findCol(xheader, {xguess{4}, xguess{3}}));
        
        xx1 = t_xx1.Variables;
        xx2 = t_xx2.Variables;
        
        fprintf(xout_txt, '#### suppress vs. clear: target vs. nontarg\n');
        
        xx = vertcat(xx1, xx2);
        
        for xtarg = 1:2
            fprintf(xout_txt, '%s: M=%4.4f, SE=%4.4f\n', targ_names{xtarg}, ...
                mean(xx(:,xtarg)), std(xx(:,xtarg))/sqrt(length(xsub_groups)));
        end
        
        fprintf(xout_txt,'\n');
        
        clear xpvalue xstats
        
        [~, xpvalue, ~, xstats] = ttest(xx(:,1), xx(:,2), 'Alpha', 0.05);
        
        fprintf(xout_txt, 'targ vs. nontarg: T(%s)=%4.4f, P=%4.4f\n', ...
            num2str(xstats.df), xstats.tstat, xpvalue);
        
    end
    
    %% ============= COMPARING CONDITIONS: replace(2) vs. two_removals (3,4)
    clear xconfusion
    xmeas_names = {'accuracy','evidence'};
    
    fprintf(xout_txt, '\n\n####################################################################\n');
    fprintf(xout_txt, '#### Replace vs. Removals: paired T test\n');
    fprintf(xout_txt, '####################################################################\n\n');
    
    for xmeas = 1:length(xmeas_names)%1_auc, 2_auc_baseline, 3_4_classifier accuracy, 4_classifier evidence
        clear xtable xheader xtarget xguess xx1 xx2
        
        fprintf(xout_txt, '\n=========================================\n');
        fprintf(xout_txt, '======== %s\n\n', xmeas_names{xmeas});
        
        %*************** data matrix
        xcsv = fullfile(dirs.mvpa.group.out{xph}, ...
            sprintf('grp_%s_%s_n%s.csv', xmeas_names{xmeas}, basename, num2str(n_subs)));
        
        xtable  = readtable(xcsv);
        xheader = xtable.Properties.VariableNames;
        
        xtarget = unique(xtable.Row)';
        xguess  = xtable.Properties.VariableNames(~ismember(xtable.Properties.VariableNames, 'Row'));
        
        for xtarg = 1:length(xtarget)
            xunit = ismember(xtable.Row, xtarget{xtarg});
            xconfusion{xtarg} = xtable(xunit, :);
        end
        
        %*************** 1_target (actual), 2_guessed (predicted)
        % replace(2) vs. two_removals (3,4)
        t_xx1 = xconfusion{2}(:, findCol(xheader, {xguess{2}, xguess{3}, xguess{4}}));
        t_xx2 = xconfusion{3}(:, findCol(xheader, {xguess{3}, xguess{2}}));
        t_xx3 = xconfusion{4}(:, findCol(xheader, {xguess{4}, xguess{2}}));
        
        xx1 = t_xx1.Variables; xx1 = [xx1(:,1) mean(xx1(:,2:3),2)];
        xx2 = t_xx2.Variables;
        xx3 = t_xx3.Variables;
        
        xx = vertcat(xx1, xx2, xx3);
        
        fprintf(xout_txt, '#### suppress vs. clear: target vs. nontarg\n');
        
        for xtarg = 1:2
            
            fprintf(xout_txt, '%s: M=%4.4f, SE=%4.4f\n', targ_names{xtarg}, ...
                mean(xx(:,xtarg)), std(xx(:,xtarg))/sqrt(length(xsub_groups)));
        end
        
        fprintf(xout_txt,'\n');
        
        clear xpvalue xstats
        
        [~, xpvalue, ~, xstats] = ttest(xx(:,1), xx(:,2), 'Alpha', 0.05);
        
        fprintf(xout_txt, 'targ vs. nontarg: T(%s)=%4.4f, P=%4.4f\n', ...
            num2str(xstats.df), xstats.tstat, xpvalue);
        
    end
end

end