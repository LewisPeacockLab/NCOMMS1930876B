function[mvpaout] = parse_mvpa_results_study(args, ph, dirs)

%---------------------------------------------------------------------------
%*************** parse the results
%---------------------------------------------------------------------------
% timecourse: baseline = other categories except "rest" class

xph = args.xphase;
fprintf('\n(+) parse the results: %s\n', args.phase_name{xph});

xindex = args.index{xph};% param index from study

%*************** the study volume was not shifted but trimmed: 
%*************** args.index{xph}.matrix was not shifted/ only trimmed
%*************** args.regs{xph}.regressor was shifted for running classification
%*************** non-shifting for timecourse

decode.matrix    = xindex.matrix; %trimmed/non-shifted regressor
decode.header    = xindex.header;
n_header         = length(decode.header);

%% ============= GRAB ACTIVATIONS & PERFORMANCE
% xguess   = ph.results.iterations.perfmet.guesses; % only for perception
% xdesired = ph.results.iterations.perfmet.desireds;
% xcorrect = ph.results.iterations.perfmet.corrects;
% xacts    = ph.results.iterations.acts;

mvpaout.total_perf = ph.results.total_perf;
out_guess          = ph.results.iterations.perfmet.guesses;
out_acts           = ph.results.iterations.acts;

%% ============= UNPACK PARAMETERS
xparam             = xindex.param;
xcond_name         = xparam.conds_names;
n_condition        = length(xcond_name);
n_runs             = xparam.n_runs;
n_trials           = xparam.n_trials;
n_subcategory      = xparam.n_subcategory;
n_classes          = max(ph.results.iterations.perfmet.desireds);
dur_stim           = xparam.dur_stim;
dur_oper           = xparam.dur_manipulation;
dur_fix            = xparam.dur_iti;
n_trs              = args.tc_tr;

%*************** setup category index
subcate_index      = 1:n_subcategory;

if strcmp(args.level, 'category')
    class_name = xparam.category_name(args.selected_category);
elseif strcmp(args.level, 'subcategory')
    for it_cate = 1:length(args.selected_category)
        xcate = args.selected_category(it_cate);
        for xsubcate = 1:n_subcategory
            xunit = xsubcate + (xparam.n_subcategory * (it_cate-1));
            
            class_name{xunit} = xparam.subcategory_name{xcate}{xsubcate};
        end
    end
end

class_index = 1:length(class_name);

%*************** base name
base_name = args.analysis_basename;

%*************** xselected_class (out of 12)
if args.class_selecting
    
    if strcmp(args.level, 'category')
        xselected_class = args.selected_category;
        rest_class      = length(xselected_class) + 1;
        
    elseif strcmp(args.level, 'subcategory')

        xselected_class = [];
        
        for it_cate = args.selected_category
            xselected_class = horzcat(xselected_class, (1:n_subcategory) + n_subcategory * (it_cate-1));
        end
        
        rest_class = (length(args.selected_category) * n_subcategory) + 1;
    end
else
    xselected_class = 1:n_classes;
end

%% ============= RESET REGRESSOR
%*************** desired for class numbers (1:n_class)
decode.header{n_header + 1} = 'desired_class'; % xselected_class + rest_class
decode.header{n_header + 2} = 'desired_sh_class'; 
decode.matrix(findCol(decode.header, {'desired_class'}), :) = 0;
decode.matrix(findCol(decode.header, {'desired_sh_class'}), :) = 0;

out_desireds    = [];
sh_out_desireds = [];

for xrun = 1:n_runs
    tunit    = find(getDATA(decode.matrix', decode.header, {'run'}, {xrun}))';
    tcate    = decode.matrix(findCol(decode.header, {'category'}), tunit);
    tsubcate = decode.matrix(findCol(decode.header, {'subcategory'}), tunit);
    
    if strcmp(args.level, 'category')
        tclass = tcate;
    elseif strcmp(args.level, 'subcategory')
        tclass = tsubcate + ((tcate-1) * n_subcategory);
        tclass(tclass < 0) = 0;
    end
    
    out_desireds    = horzcat(out_desireds, tclass);
    sh_out_desireds = horzcat(sh_out_desireds, [zeros(1,args.shift_TRs) tclass(1:end-args.shift_TRs)]);
end

decode.matrix(findCol(decode.header, {'desired_class'}), :)    = out_desireds;
decode.matrix(findCol(decode.header, {'desired_sh_class'}), :) = sh_out_desireds;

xtable = array2table(decode.matrix', 'VariableNames', decode.header);
writetable(xtable, fullfile(dirs.param, sprintf('study_%s_volume_matrix_reset_%s.csv', args.level, args.subject_id)));

%% ============= SETUP STRUCTURE
% xdesired: 1_9 subcategory
% decode.classifier.guessed{xdesired}(xguess)
% decode.classifier.evidence{xdesired}(xguess)
% decode.classifier.accuracy(xdesired)

for xdesired = 1:n_classes
    decode.classifier.guessed{xdesired}  = [];
    decode.classifier.evidence{xdesired} = [];
    decode.classifier.accuracy(xdesired) = 0;
end

%% ============= CLASSIFICATION RESULTS
% xdesired: 1_9 subcategory
% decode.classifier.guessed{xdesired}(xguess)
% decode.classifier.evidence{xdesired}(xguess)
% decode.classifier.accuracy(xdesired)
% out_guess, out_acts

for xdesired = 1:n_classes
    
    tunit = find(getDATA(decode.matrix', decode.header,...
        {'desired_sh_class'}, {xdesired}));
    xspike = ~(decode.matrix(findCol(decode.header,{'spike'}), tunit));
    
    if sum(xspike~=0)
        xunit      = tunit(xspike);
        xguessed   = out_guess(xunit);
        act_matrix = out_acts(:,xunit);
        xevidence  = mean(act_matrix, 2);
        
        for xguess = 1:n_classes
            %*************** guessed category from classifier
            decode.classifier.guessed{xdesired}(xguess)  = numel(find(xguessed == xguess));
            
            %*************** evidence of the classifier
            decode.classifier.evidence{xdesired}(xguess) = xevidence(xguess);
        end
        
        %*************** accuracy of the classifier
        decode.classifier.accuracy(xdesired) = ...
            (decode.classifier.guessed{xdesired}(xdesired))/...
            (sum(decode.classifier.guessed{xdesired}));
    end
end

%% ============= ANALYSIS PARAMETERS
%*************** selected volxels
mvpaout.param.n_total_voxels     = ph.subj.masks{1}.nvox;
mvpaout.param.n_selected_voxels  = ph.subj.masks{3}.nvox;
mvpaout.param.n_selected_percent = (ph.subj.masks{3}.nvox/ph.subj.masks{1}.nvox)*100;

%% ============= VERIFYING CLASSIFICATION OF STUDY DATA
%*************** rest
if strcmp(args.rest, 'rest'), class_name{n_classes} = 'rest'; end

clear xaccuracy xevidence xtarg_accuracy target_names guess_names

for xcate = 1:n_classes
    target_names(xcate) = strcat(sprintf('t%s_', num2str(xcate)), class_name(xcate)); %#ok<*AGROW>
    guess_names(xcate)  = strcat(sprintf('g%s_', num2str(xcate)), class_name(xcate));
end

%*************** create tables
xaccuracy      = zeros(n_classes, n_classes);
xevidence      = zeros(n_classes, n_classes);
xtarg_accuracy = zeros(1, n_classes);
xtarg_evidence = zeros(1, n_classes);

for xcate = 1:n_classes
    if ~isempty(decode.classifier.guessed{xcate})
        xaccuracy(xcate,:)    = (decode.classifier.guessed{xcate})/(sum(decode.classifier.guessed{xcate}));
        xevidence(xcate,:)    = decode.classifier.evidence{xcate};
        xtarg_accuracy(xcate) = xaccuracy(xcate, xcate);
        xtarg_evidence(xcate) = xevidence(xcate, xcate);%target accuracy
    end
end

%*************** total accuracy
xtarg_accuracy(end + 1) = mean(xtarg_accuracy);
xtarg_evidence(end + 1) = mean(xtarg_evidence);

%*************** accuracy
table_accuracy = array2table(xaccuracy,...
    'RowNames', target_names, 'VariableNames', guess_names);

csv_name = sprintf('%s/table_accuracy_%s.csv', dirs.mvpa.parse{xph}, base_name);
writetable(table_accuracy, csv_name,'WriteRowNames',true)

%*************** evidence
table_evidence = array2table(xevidence,...
    'RowNames', target_names, 'VariableNames', guess_names);

csv_name = sprintf('%s/table_evidence_%s.csv', dirs.mvpa.parse{xph}, base_name);
writetable(table_evidence, csv_name,'WriteRowNames',true)

%*************** target accuracy
table_targ_accuracy = array2table(xtarg_accuracy,...
    'VariableNames', [target_names 'total']);

csv_name = sprintf('%s/table_target_accuracy_%s.csv', dirs.mvpa.parse{xph}, base_name);
writetable(table_targ_accuracy, csv_name,'WriteRowNames',true)

%*************** target evidence
table_targ_evidence = array2table(xtarg_evidence,...
    'VariableNames', [target_names 'total']);

csv_name = sprintf('%s/table_target_evidence_%s.csv', dirs.mvpa.parse{xph}, base_name);
writetable(table_targ_evidence, csv_name,'WriteRowNames',true)

%*************** selected volxels
table_selected_voxels = table(mvpaout.param.n_total_voxels, ...
    mvpaout.param.n_selected_voxels, mvpaout.param.n_selected_percent, ...
    'VariableNames', {'total_voxels', 'selected_voxels', 'percentage'}');

csv_name = sprintf('%s/table_n_voxels_%s.csv', dirs.mvpa.parse{xph}, base_name);
writetable(table_selected_voxels, csv_name,'WriteRowNames',true)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%*************** TIMECOURSE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ============= SETUP STRUCTURE
%*************** working memory contents
%--------------- subcategory level
% 1. maintain            : {1} subtarget, {2} nonsubtarget(2), {3} mean(nontarget(6))
% 2. replace category    : {1} subtarget, {2} nonsubtarget(2), {3} new_subtarget, {4} new_nonsubtarget(2), {5} mean(nontarget(3))
% 3. replace subcategory : {1} subtarget, {2} nonsubtarget(1), {3} new_subtarget(1), {4} mean(nontarget(6))
% 4. target suppress     : {1} subtarget, {2} nonsubtarget(2), {3} mean(nontargets(6))
% 5. global clear        : {1} subtarget, {2} nonsubtarget(2), {3} mean(nontargets(6))
%--------------- category level
% 1. maintain            : {1} target, {2} mean(nontargets)
% 2. replace category    : {1} target, {2} new category, {3} nontargets
% 3. replace subcategory : {1} target, {2} mean(nontargets)
% 4. target suppress     : {1} target, {2} mean(nontargets)
% 5. global clear        : {1} target, {2} mean(nontargets)

for xcond = 1:n_condition
    if strcmp(args.level,'subcategory')
        if xcond==2, n_targ = 5; elseif xcond==3, n_targ = 4; else n_targ = 3; end
    elseif strcmp(args.level,'category')
        if xcond==2, n_targ = 3; else n_targ = 2; end
    end
    
    %*************** bined trs
    for xbin = 1:2
        if xbin==1, it_ntrs = dur_stim + dur_oper + max(dur_fix); 
        else        it_ntrs = n_trs - (dur_stim + dur_oper + min(dur_fix)); end
        
        for xtr = 1:it_ntrs
            for xtarg = 1:n_targ
                decode.bins{xbin}.condition{xcond}.target{xtarg}.evidence{xtr} = [];
            end
        end
    end
    
    for xtr = 1:n_trs
        for xtarg = 1:n_targ
            decode.timecourse.condition{xcond}.target{xtarg}.evidence{xtr} = [];
        end
        
        for xcate = 1:n_classes
            for xx = 1:n_classes
                
                decode.cate_timecourse.condition{xcond}.category{xcate}.target{xx}.evidence{xtr} = [];
                
                for xnewcate = 1:n_classes
                    decode.cate_timecourse.condition{xcond}.category{xcate}.newcategory{xnewcate}.target{xx}.evidence{xtr} = [];\
                end
            end
        end
    end
end

%% ============= COLLECTING EVIDENCE: PER CONDITION
%*************** working memory contents
% timecourse: un-shifted tr: -args.shift_TRs
% condition: 1_maintain,2_replace_category,3_replace_item,4_target_suppress,5_global_clear
% xtarg: 1_targ, 2_non/newtarg(switch), 3_residuals(for cond_2, 3)
% xtr: 12 presentation + 5:9 fixation = 21 max
% decode.timecourse.condition{xcond}.target{xtarg}.evidence{xtr}

cond_1=[]; cond_2=[]; cond_3=[]; cond_4=[]; cond_5=[];

if strcmp(args.level, 'subcategory')
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %--------------- subcategory level
    % 1. maintain            : {1} subtarget, {2} nonsubtarget(2), {3} mean(nontarget(6))
    % 2. replace category    : {1} subtarget, {2} nonsubtarget(2), {3} new_subtarget(1), {4} new_nonsubtarget(2), {5} mean(nontarget(3))
    % 3. replace subcategory : {1} subtarget, {2} nonsubtarget(1), {3} new_subtarget(1), {4} mean(nontarget(6))
    % 4. target suppress     : {1} subtarget, {2} nonsubtarget(2), {3} mean(nontargets(6))
    % 5. global clear        : {1} subtarget, {2} nonsubtarget(2), {3} mean(nontargets(6))
    
    for xrun = 1:n_runs
        for xtrial = 1:n_trials(xrun)
            clear xsubtarg xnonsubtarg xnew_subtarg xnew_nonsubtarg xnontarg
            xunit_trial = find(getDATA(decode.matrix', decode.header,...
                {'run','trial'}, {xrun, xtrial}));
            xunit_stim  = xunit_trial(1:dur_stim);
            xunit_manip = xunit_trial(~ismember(xunit_trial, xunit_stim));
            
            xcond       = unique(decode.matrix(findCol(decode.header,{'condition'}), xunit_trial));
            xcate       = unique(decode.matrix(findCol(decode.header,{'category'}), xunit_stim));
            xsubcate    = unique(decode.matrix(findCol(decode.header,{'subcategory'}), xunit_stim));
            
            %*************** number of cells
            if xcond==2, n_targ = 5; elseif xcond==3, n_targ = 4; else n_targ = 3; end
            
            %*************** target: subtarget(1)
            xsubtarg    = xsubcate + (n_subcategory * (xcate-1));
            
            %*************** new target: new_subtarget(1)
            if (xcond==2 || xcond==3)
                xnew_cate    = unique(decode.matrix(findCol(decode.header,{'new_category'}), xunit_manip));
                xnew_subcate = unique(decode.matrix(findCol(decode.header,{'new_subcategory'}), xunit_manip));
                
                xnew_subtarg = xnew_subcate + (n_subcategory * (xnew_cate-1));
            end
            
            %*************** nonsubtarget
            if xcond==3 % replace subcategory: nonsubtarget(1)
                xnonsubcate = subcate_index(~ismember(subcate_index, [xsubcate, xnew_subcate]));
            else % nonsubtarget(2)
                xnonsubcate = subcate_index(~ismember(subcate_index, xsubcate));
            end
            
            xnonsubtarg = xnonsubcate + (n_subcategory * (xcate-1));
            
            %*************** new nonsubtarget
            if xcond==2 % replace category: new_nonsubtarget(2)
                xnew_nonsubcate = subcate_index(~ismember(subcate_index, xnew_subcate));
                xnew_nonsubtarg = xnew_nonsubcate + (n_subcategory * (xnew_cate-1));
            end
            
            %*************** mean(nontarg)
            if xcond==2
                xnontarg = class_index(~ismember(class_index, ...
                    [xsubtarg, xnonsubtarg, xnew_subtarg, xnew_nonsubtarg]));
            elseif xcond==3
                xnontarg = class_index(~ismember(class_index, ...
                    [xsubtarg, xnonsubtarg, xnew_subtarg]));
            else
                xnontarg = class_index(~ismember(class_index, ...
                    [xsubtarg, xnonsubtarg]));
            end
            
            %*************** verification
            if xcond==1
                cond_1 = horzcat(cond_1, [xsubtarg, xnonsubtarg, xnontarg]);
            elseif xcond==2
                cond_2 = horzcat(cond_2, [xsubtarg, xnonsubtarg, xnew_subtarg, xnew_nonsubtarg, xnontarg]);
            elseif xcond==3
                cond_3 = horzcat(cond_3, [xsubtarg, xnonsubtarg, xnew_subtarg, xnontarg]);
            elseif xcond==4
                cond_4 = horzcat(cond_4, [xsubtarg, xnonsubtarg, xnontarg]);
            elseif xcond==5
                cond_5 = horzcat(cond_5, [xsubtarg, xnonsubtarg, xnontarg]);
            end
            
            %********************************
            %*************** collect evidence
            on_unit   = xunit_stim(1);% - args.shift_TRs;
            off_unit  = on_unit + n_trs - 1;%xunit_manip(1) + dur_manip + hrf_trs - 1;% - args.shift_TRs;
            
            xtc_unit  = on_unit:off_unit;
            xtc_spike = ~(decode.matrix(findCol(decode.header,{'spike'}), xtc_unit));
            
            for xtr = 1:length(xtc_unit)
                if xtc_spike(xtr)~=0
                    act_matrix = out_acts(:, xtc_unit(xtr));
                    
                    %*************** subtarget: presented {1}
                    decode.timecourse.condition{xcond}.target{1}.evidence{xtr} = ...
                        horzcat(decode.timecourse.condition{xcond}.target{1}.evidence{xtr},...
                        act_matrix(xsubtarg));
                    
                    %*************** nonsubtarget {2}
                    decode.timecourse.condition{xcond}.target{2}.evidence{xtr} = ...
                        horzcat(decode.timecourse.condition{xcond}.target{2}.evidence{xtr},...
                        mean(act_matrix(xnonsubtarg)));
                    
                    %*************** new subtarget {3}
                    if (xcond==2) || (xcond==3)
                        if sum(ismember(xselected_class, xnew_subtarg))
                            decode.timecourse.condition{xcond}.target{3}.evidence{xtr} = ...
                                horzcat(decode.timecourse.condition{xcond}.target{3}.evidence{xtr},...
                                act_matrix(xnew_subtarg));
                        end
                    end
                    
                    %*************** new nonsubtarget {4}
                    if xcond==2
                        decode.timecourse.condition{xcond}.target{4}.evidence{xtr} = ...
                            horzcat(decode.timecourse.condition{xcond}.target{4}.evidence{xtr},...
                            mean(act_matrix(xnew_nonsubtarg)));
                    end
                    
                    %*************** nontarg {n_targ}
                    decode.timecourse.condition{xcond}.target{n_targ}.evidence{xtr} = ...
                        horzcat(decode.timecourse.condition{xcond}.target{n_targ}.evidence{xtr},...
                        mean(act_matrix(xnontarg)));
                    
                    %*************** bined trs
                    if xtr <= length(xunit_trial)
                        xbin = 1; it_tr = xtr;
                    else
                        xbin = 2; it_tr = xtr - length(xunit_trial);
                    end
                    
                    %*************** subtarget: presented {1}
                    decode.bins{xbin}.condition{xcond}.target{1}.evidence{it_tr} = ...
                        horzcat(decode.bins{xbin}.condition{xcond}.target{1}.evidence{it_tr},...
                        act_matrix(xsubtarg));
                    
                    %*************** nonsubtarget {2}
                    decode.bins{xbin}.condition{xcond}.target{2}.evidence{it_tr} = ...
                        horzcat(decode.bins{xbin}.condition{xcond}.target{2}.evidence{it_tr},...
                        mean(act_matrix(xnonsubtarg)));
                    
                    %*************** new subtarget {3}
                    if (xcond==2) || (xcond==3)
                        if sum(ismember(xselected_class, xnew_subtarg))
                            decode.bins{xbin}.condition{xcond}.target{3}.evidence{it_tr} = ...
                                horzcat(decode.bins{xbin}.condition{xcond}.target{3}.evidence{it_tr},...
                                act_matrix(xnew_subtarg));
                        end
                    end
                    
                    %*************** new nonsubtarget {4}
                    if xcond==2
                        decode.bins{xbin}.condition{xcond}.target{4}.evidence{it_tr} = ...
                            horzcat(decode.bins{xbin}.condition{xcond}.target{4}.evidence{it_tr},...
                            mean(act_matrix(xnew_nonsubtarg)));
                    end
                    
                    %*************** nontarg {n_targ}
                    decode.bins{xbin}.condition{xcond}.target{n_targ}.evidence{it_tr} = ...
                        horzcat(decode.bins{xbin}.condition{xcond}.target{n_targ}.evidence{it_tr},...
                        mean(act_matrix(xnontarg)));
                end
            end
        end
    end
elseif strcmp(args.level, 'category')
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %--------------- category level
    % 1. maintain            : {1} target, {2} mean(nontargets)
    % 2. replace category    : {1} target, {2} new category, {3} nontargets
    % 3. replace subcategory : {1} target, {2} mean(nontargets)
    % 4. target suppress     : {1} target, {2} mean(nontargets)
    % 5. global clear        : {1} target, {2} mean(nontargets)
    
    for xrun = 1:n_runs
        for xtrial = 1:n_trials(xrun)
            xunit_trial = find(getDATA(decode.matrix', decode.header,...
                {'run','trial'}, {xrun, xtrial}));
            xunit_stim  = xunit_trial(1:dur_stim);
            xunit_manip = xunit_trial(~ismember(xunit_trial, xunit_stim));
            
            xcond       = unique(decode.matrix(findCol(decode.header,{'condition'}), xunit_trial));
            xcate       = unique(decode.matrix(findCol(decode.header,{'category'}), xunit_stim));
            
            %*************** number of cells
            if xcond==2, n_targ = 3; else n_targ = 2; end
            
            %*************** target
            xtarg       = xcate;
            
            %*************** new target
            if xcond==2    
                xnew_cate = unique(decode.matrix(findCol(decode.header,{'new_category'}), xunit_manip));
                xnew_targ = xnew_cate;
            end
            
            %*************** mean(nontarg)
            if xcond==2
                xnontarg = class_index(~ismember(class_index, [xtarg, xnew_targ]));
            else
                xnontarg = class_index(~ismember(class_index, xtarg));
            end
            
            %********************************
            %*************** collect evidence
            on_unit   = xunit_stim(1);% - args.shift_TRs;
            off_unit  = on_unit + n_trs - 1;%xunit_manip(1) + dur_manip + hrf_trs - 1;% - args.shift_TRs;
            
            xtc_unit  = on_unit:off_unit;
            xtc_spike = ~(decode.matrix(findCol(decode.header,{'spike'}), xtc_unit));
            
            for xtr = 1:length(xtc_unit)
                if xtc_spike(xtr)~=0
                    act_matrix = out_acts(:, xtc_unit(xtr));
                    
                    %*************** target: presented {1}
                    decode.timecourse.condition{xcond}.target{1}.evidence{xtr} = ...
                        horzcat(decode.timecourse.condition{xcond}.target{1}.evidence{xtr},...
                        act_matrix(xtarg));
                    
                    %*************** new target {2}
                    if xcond==2
                        decode.timecourse.condition{xcond}.target{2}.evidence{xtr} = ...
                            horzcat(decode.timecourse.condition{xcond}.target{2}.evidence{xtr},...
                            act_matrix(xnew_targ));
                    end
                    
                    %*************** nontarg {n_targ}
                    decode.timecourse.condition{xcond}.target{n_targ}.evidence{xtr} = ...
                        horzcat(decode.timecourse.condition{xcond}.target{n_targ}.evidence{xtr},...
                        mean(act_matrix(xnontarg)));
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%
                    %*************** bined trs
                    if xtr <= length(xunit_trial)
                        xbin = 1; it_tr = xtr;
                    else
                        xbin = 2; it_tr = xtr - length(xunit_trial);
                    end
                    
                    %*************** target: presented {1}
                    decode.bins{xbin}.condition{xcond}.target{1}.evidence{it_tr} = ...
                        horzcat(decode.bins{xbin}.condition{xcond}.target{1}.evidence{it_tr},...
                        act_matrix(xtarg));
                    
                    %*************** new target {2}
                    if xcond==2
                        decode.bins{xbin}.condition{xcond}.target{2}.evidence{it_tr} = ...
                            horzcat(decode.bins{xbin}.condition{xcond}.target{2}.evidence{it_tr},...
                            act_matrix(xnew_targ));
                    end
                    
                    %*************** nontarg {n_targ}
                    decode.bins{xbin}.condition{xcond}.target{n_targ}.evidence{it_tr} = ...
                        horzcat(decode.bins{xbin}.condition{xcond}.target{n_targ}.evidence{it_tr},...
                        mean(act_matrix(xnontarg)));
                end
            end
        end
    end
end

%% ============= SAVE MVPAOUT
mvpaout.decode = decode;

end
