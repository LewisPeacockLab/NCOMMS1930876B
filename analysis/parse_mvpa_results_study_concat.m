function[mvpaout] = parse_mvpa_results_study_concat(args, xresult, dirs)

%---------------------------------------------------------------------------
%*************** parse the results
%---------------------------------------------------------------------------
% subcategory-level / norest
% args.class_selecting = 1
% timecourse: baseline = other categories except "rest" class
% regressor was shifted

xph = args.xphase;
fprintf('\n(+) parse the results: %s\n', args.phase_name{xph});

xindex = args.index{xph};% param index from study

%% ============= UNPACK PARAMETERS
xparam             = xindex.param;
xcond_name         = xparam.conds_names;
n_condition        = length(xcond_name);
n_runs             = xparam.n_runs;
n_trials           = xparam.n_trials;
n_category         = xparam.n_category;
n_subcategory      = xparam.n_subcategory;
n_classes          = n_category * n_subcategory;
dur_stim           = xparam.dur_stim;
dur_manip          = xparam.dur_manipulation;
hrf_trs            = xparam.hrf_trs;
n_tc_trs           = xparam.n_tc_trs;

%*************** setup category index
subcate_index      = 1:n_subcategory;
class_index        = 1:n_classes;

for xcate = 1:n_category
    for xsubcate = 1:n_subcategory
        class_name{xsubcate + (n_subcategory * (xcate-1))} = ...
            xparam.subcategory_name{xcate}{xsubcate};
    end
end

%*************** base name
base_name = args.analysis_basename;

%% ============= GRAB ACTIVATIONS & PERFORMANCE
%*************** the study volume was not shifted but trimmed: 
%*************** args.index{xph}.matrix was not shifted/ only trimmed
%*************** args.regs{xph}.regressor was shifted for running classification
%*************** non-shifting for timecourse

decode.matrix    = xindex.matrix; %trimmed/non-shifted regressor
decode.header    = xindex.header;
n_header         = length(decode.header);

% xguess   = ph.results.iterations.perfmet.guesses; % only for perception
% xdesired = ph.results.iterations.perfmet.desireds;
% xcorrect = ph.results.iterations.perfmet.corrects;
% xacts    = ph.results.iterations.acts;
out_guess = []; out_acts = []; out_regs = [];
out_guess_each = []; out_regs_each = [];%for each category separately

for xcate = 1:3
    mvpaout.total_perf(xcate) = xresult{xcate}.results.total_perf;
    out_guess(xcate, :) = ...
        xresult{xcate}.results.iterations.perfmet.guesses + (n_subcategory * (xcate-1));
    out_acts = vertcat(out_acts, xresult{xcate}.results.iterations.acts);
    
    out_regs = vertcat(out_regs, ...
        (xresult{xcate}.results.iterations.perfmet.desireds) + (n_subcategory * (xcate-1)));
    
    %*************** each category
    out_guess_each(xcate, :) = ...
        xresult{xcate}.results.iterations.perfmet.guesses;
    out_regs_each = vertcat(out_regs_each, ...
        xresult{xcate}.results.iterations.perfmet.desireds);
end

%% ============= RESET REGRESSOR
%*************** desired for class numbers (1:n_class)
decode.header{n_header + 1} = 'desired_class'; % xselected_class + rest_class
decode.header{n_header + 2} = 'desired_sh_class'; 
decode.header{n_header + 3} = 'desired_sh_cate'; 
decode.header{n_header + 4} = 'desired_sh_subcate'; 

out_desireds    = [];
sh_out_desireds = [];
sh_out_cate     = []; 
sh_out_subcate  = []; 

for xrun = 1:n_runs
    xunit         = find(getDATA(decode.matrix', decode.header, {'run'}, {xrun}))';
    cate_array    = decode.matrix(findCol(decode.header, {'category'}), xunit);
    subcate_array = decode.matrix(findCol(decode.header, {'subcategory'}), xunit);
    
    xclass_array = subcate_array + (n_subcategory * (cate_array-1));
    xclass_array(xclass_array < 0) = 0;
    
    out_desireds    = horzcat(out_desireds, xclass_array);
    sh_out_desireds = horzcat(sh_out_desireds, ...
        [zeros(1,args.shift_TRs) xclass_array(1:end-args.shift_TRs)]);
    sh_out_cate     = horzcat(sh_out_cate, ...
        [zeros(1,args.shift_TRs) cate_array(1:end-args.shift_TRs)]);
    sh_out_subcate  = horzcat(sh_out_subcate, ...
        [zeros(1,args.shift_TRs) subcate_array(1:end-args.shift_TRs)]);
end

decode.matrix(findCol(decode.header, {'desired_class'}), :)      = out_desireds;
decode.matrix(findCol(decode.header, {'desired_sh_class'}), :)   = sh_out_desireds;
decode.matrix(findCol(decode.header, {'desired_sh_cate'}), :)    = sh_out_cate;
decode.matrix(findCol(decode.header, {'desired_sh_subcate'}), :) = sh_out_subcate;

%% ============= SETUP STRUCTURE
% xdesired: 1_9 subcategory, 10: rest
% decode.classifier.guessed{xdesired}(xguess)
% decode.classifier.evidence{xdesired}(xguess)
% decode.classifier.accuracy(xdesired)

for xdesired = 1:n_classes
    decode.classifier.guessed{xdesired}  = zeros(1,n_subcategory);
    decode.classifier.evidence{xdesired} = zeros(1,n_subcategory);
    decode.classifier.accuracy(xdesired) = 0;
end

%% ============= CLASSIFICATION RESULTS
% xdesired: 1_9 subcategory, 10: rest
% decode.classifier.guessed{xdesired}(xguess)
% decode.classifier.evidence{xdesired}(xguess)
% decode.classifier.accuracy(xdesired)
% out_guess, out_acts

for xcate = 1:n_category
    selected_evidence = out_acts((1:3) + (n_subcategory * (xcate-1)), :);
    
    for xsubcate = 1:n_subcategory
        xdesired = xsubcate + (n_subcategory * (xcate-1));
        tunit = find(getDATA(decode.matrix', decode.header,...
            {'desired_sh_class'}, {xdesired}));
        
        xspike = ~(decode.matrix(findCol(decode.header,{'spike'}), tunit));
        
        if sum(xspike~=0)
            xunit     = tunit(xspike);
            xguessed  = out_guess_each(xcate, xunit);
            xevidence = mean(selected_evidence(:, xunit), 2);
            
            for xguess = 1:n_subcategory
                %*************** guessed category from classifier
                decode.classifier.guessed{xdesired}(xguess) = numel(find(xguessed == xguess));
                
                %*************** evidence of the classifier
                decode.classifier.evidence{xdesired}(xguess) = xevidence(xguess);
            end
            
            %*************** accuracy of the classifier
            decode.classifier.accuracy(xdesired) = ...
                (decode.classifier.guessed{xdesired}(xsubcate))/...
                (sum(decode.classifier.guessed{xdesired}));
        end
    end
end

%% ============= ANALYSIS PARAMETERS
%*************** selected volxels
for xcate = 1:n_category
    mvpaout.param{xcate}.n_total_voxels     = xresult{xcate}.masks.n_total_voxels;
    mvpaout.param{xcate}.n_selected_voxels  = xresult{xcate}.masks.n_selected_voxels;
    mvpaout.param{xcate}.n_selected_percent = xresult{xcate}.masks.n_selected_percent;
end

%% ============= VERIFYING CLASSIFICATION OF STUDY DATA
clear xaccuracy xevidence xtarg_accuracy target_names guess_names

for xclass = 1:n_classes
    target_names(xclass) = strcat(sprintf('t%s_', num2str(xclass)), class_name(xclass)); %#ok<*AGROW>
    guess_names(xclass)  = strcat(sprintf('g%s_', num2str(xclass)), class_name(xclass));
end

%*************** create tables
xaccuracy = zeros(n_classes, n_classes);
xevidence = zeros(n_classes, n_classes);
xtarg_accuracy = zeros(1, n_classes);

for xcate = 1:n_category
    for xsubcate = 1:n_subcategory
        xclass = xsubcate + (n_subcategory * (xcate-1));
        xunit  = (1:n_subcategory) + (n_subcategory * (xcate-1));
        xaccuracy(xclass,xunit) = (decode.classifier.guessed{xclass})/(sum(decode.classifier.guessed{xclass}));
        xevidence(xclass,xunit) = decode.classifier.evidence{xclass};
        xtarg_accuracy(xclass)  = xaccuracy(xclass, xclass);
    end
end

%*************** total accuracy
xtarg_accuracy(end + 1) = mean(xtarg_accuracy);

%*************** accuracy
table_accuracy = array2table(xaccuracy,...
    'RowNames', target_names, 'VariableNames', guess_names);

csv_name = sprintf('%s/table_percept_acc_%s.csv', dirs.mvpa.parse{xph}, base_name);
writetable(table_accuracy, csv_name,'WriteRowNames',true)

%*************** evidence
table_evidence = array2table(xevidence,...
    'RowNames', target_names, 'VariableNames', guess_names);

csv_name = sprintf('%s/table_percept_evi_%s.csv', dirs.mvpa.parse{xph}, base_name);
writetable(table_evidence, csv_name,'WriteRowNames',true)

%*************** target accuracy
table_targ_accuracy = array2table(xtarg_accuracy,...
    'VariableNames', [class_name 'total']);

csv_name = sprintf('%s/table_percept_targ_acc_%s.csv', dirs.mvpa.parse{xph}, base_name);
writetable(table_targ_accuracy, csv_name,'WriteRowNames',true)

%*************** selected volxels
total_vox = []; sel_vox = []; perce_vox = [];
for xcate = 1:n_category
    total_vox = vertcat(total_vox, mvpaout.param{xcate}.n_total_voxels);
    sel_vox   = vertcat(sel_vox, mvpaout.param{xcate}.n_selected_voxels);
    perce_vox = vertcat(perce_vox, mvpaout.param{xcate}.n_selected_voxels);
end

table_selected_voxels = table(total_vox, sel_vox, perce_vox, ...
    'VariableNames', {'total_voxels', 'selected_voxels', 'percentage'}',...
    'RowNames', xparam.category_name);

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

for xcond = 1:n_condition
    if strcmp(args.level,'subcategory')
        if xcond==2, n_targ = 5; elseif xcond==3, n_targ = 4; else n_targ = 3; end
    elseif strcmp(args.level,'category')
        if xcond==2, n_targ = 3; else n_targ = 2; end
    end
    
    for xtarg = 1:n_targ
        for xtr = 1:n_tc_trs
            decode.timecourse.condition{xcond}.target{xtarg}.evidence{xtr} = [];
            decode.timecourse_each.condition{xcond}.target{xtarg}.evidence{xtr} = [];%for each category separately
        end
    end
end

%% ============= COLLECTING EVIDENCE: CATEGORY
%*************** working memory contents
% timecourse: un-shifted tr: -args.shift_TRs
% condition: 1_maintain,2_replace_category,3_replace_item,4_target_suppress,5_global_clear
% xtarg: 1_targ, 2_non/newtarg(switch), 3_residuals(for cond_2, 3)
% xtr: 12 presentation + 5:9 fixation = 21 max
% decode.timecourse.condition{xcond}.target{xtarg}.evidence{xtr}

cond_1=[]; cond_2=[]; cond_3=[]; cond_4=[]; cond_5=[];

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------- subcategory level
% 1. maintain            : {1} subtarget, {2} nonsubtarget(2), {3} mean(nontarget(6))
% 2. replace category    : {1} subtarget, {2} nonsubtarget(2), {3} new_subtarget(1), {4} new_nonsubtarget(2), {5} mean(nontarget(3))
% 3. replace subcategory : {1} subtarget, {2} nonsubtarget(1), {3} new_subtarget(1), {4} mean(nontarget(6))
% 4. target suppress     : {1} subtarget, {2} nonsubtarget(2), {3} mean(nontargets(6))
% 5. global clear        : {1} subtarget, {2} nonsubtarget(2), {3} mean(nontargets(6))

for xrun = 1:n_runs
    for xtrial = 1:n_trials
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
            cond_1 = vertcat(cond_1, [xsubtarg, xnonsubtarg, xnontarg]);
        elseif xcond==2
            cond_2 = vertcat(cond_2, [xsubtarg, xnonsubtarg, xnew_subtarg, xnew_nonsubtarg, xnontarg]);
        elseif xcond==3
            cond_3 = vertcat(cond_3, [xsubtarg, xnonsubtarg, xnew_subtarg, xnontarg]);
        elseif xcond==4
            cond_4 = vertcat(cond_4, [xsubtarg, xnonsubtarg, xnontarg]);
        elseif xcond==5
            cond_5 = vertcat(cond_5, [xsubtarg, xnonsubtarg, xnontarg]);
        end
        
        %********************************
        %*************** collect evidence
        on_unit   = xunit_stim(1);% - args.shift_TRs;
        off_unit  = xunit_manip(1) + dur_manip + hrf_trs - 1;% - args.shift_TRs;
        
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
                    decode.timecourse.condition{xcond}.target{3}.evidence{xtr} = ...
                        horzcat(decode.timecourse.condition{xcond}.target{3}.evidence{xtr},...
                        act_matrix(xnew_subtarg));
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
            end
        end
    end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SEPARATED TIMECOURSE PER CATEGORY
%--------------- subcategory level
% 1. maintain            : {1} subtarget, {2} nonsubtarget(2)
% 2. replace category    : {1} subtarget, {2} nonsubtarget(2), {3} new_subtarget(1), {4} new_nonsubtarget(2)
% 3. replace subcategory : {1} subtarget, {2} nonsubtarget(1), {3} new_subtarget(1)
% 4. target suppress     : {1} subtarget, {2} nonsubtarget(2) 
% 5. global clear        : {1} subtarget, {2} nonsubtarget(2)

for xrun = 1:n_runs
    for xtrial = 1:n_trials
        clear xsubtarg xnonsubtarg xnew_subtarg xnew_nonsubtarg xnontarg
        xunit_trial = find(getDATA(decode.matrix', decode.header,...
            {'run','trial'}, {xrun, xtrial}));
        xunit_stim  = xunit_trial(1:dur_stim);
        xunit_manip = xunit_trial(~ismember(xunit_trial, xunit_stim));
        
        xcond       = unique(decode.matrix(findCol(decode.header,{'condition'}), xunit_trial));
        xcate       = unique(decode.matrix(findCol(decode.header,{'category'}), xunit_stim));
        xsubcate    = unique(decode.matrix(findCol(decode.header,{'subcategory'}), xunit_stim));
                
        %*************** target: subtarget(1)
        xsubtarg    = xsubcate;% + (n_subcategory * (xcate-1));
        
        %*************** new target: new_subtarget(1)
        if (xcond==2 || xcond==3)
            xnew_cate    = unique(decode.matrix(findCol(decode.header,{'new_category'}), xunit_manip));
            xnew_subcate = unique(decode.matrix(findCol(decode.header,{'new_subcategory'}), xunit_manip));
            
            xnew_subtarg = xnew_subcate;% + (n_subcategory * (xnew_cate-1));
        end
        
        %*************** nonsubtarget
        if xcond==3 % replace subcategory: nonsubtarget(1)
            xnonsubcate = subcate_index(~ismember(subcate_index, [xsubcate, xnew_subcate]));
        else % nonsubtarget(2)
            xnonsubcate = subcate_index(~ismember(subcate_index, xsubcate));
        end
        
        xnonsubtarg = xnonsubcate;% + (n_subcategory * (xcate-1));
        
        %*************** new nonsubtarget
        if xcond==2 % replace category: new_nonsubtarget(2)
            xnew_nonsubcate = subcate_index(~ismember(subcate_index, xnew_subcate));
            xnew_nonsubtarg = xnew_nonsubcate;% + (n_subcategory * (xnew_cate-1));
        end
        
        %********************************
        %*************** collect evidence
        on_unit   = xunit_stim(1);% - args.shift_TRs;
        off_unit  = xunit_manip(1) + dur_manip + hrf_trs - 1;% - args.shift_TRs;
        
        xtc_unit  = on_unit:off_unit;
        xtc_spike = ~(decode.matrix(findCol(decode.header,{'spike'}), xtc_unit));
        
        for xtr = 1:length(xtc_unit)
            if xtc_spike(xtr)~=0
                selected_evidence = out_acts((1:3) + (n_subcategory * (xcate-1)), xtc_unit(xtr));
                
                if xcond==2 %repCate
                    new_selected_evidence = out_acts((1:3) + (n_subcategory * (xnew_cate-1)), xtc_unit(xtr));
                end
                                
                %*************** subtarget: presented {1}
                decode.timecourse_each.condition{xcond}.target{1}.evidence{xtr} = ...
                    horzcat(decode.timecourse_each.condition{xcond}.target{1}.evidence{xtr},...
                    selected_evidence(xsubtarg));
                
                %*************** nonsubtarget {2}
                decode.timecourse_each.condition{xcond}.target{2}.evidence{xtr} = ...
                    horzcat(decode.timecourse_each.condition{xcond}.target{2}.evidence{xtr},...
                    mean(selected_evidence(xnonsubtarg)));
                
                %*************** new subtarget {3}
                if (xcond==2) || (xcond==3)
                    decode.timecourse_each.condition{xcond}.target{3}.evidence{xtr} = ...
                        horzcat(decode.timecourse_each.condition{xcond}.target{3}.evidence{xtr},...
                        new_selected_evidence(xnew_subtarg));
                end
                
                %*************** new nonsubtarget {4}
                if xcond==2
                    decode.timecourse_each.condition{xcond}.target{4}.evidence{xtr} = ...
                        horzcat(decode.timecourse_each.condition{xcond}.target{4}.evidence{xtr},...
                        mean(new_selected_evidence(xnew_nonsubtarg)));
                end
            end
        end
    end
end

%% ============= SAVE MVPAOUT
mvpaout.decode = decode;

end
