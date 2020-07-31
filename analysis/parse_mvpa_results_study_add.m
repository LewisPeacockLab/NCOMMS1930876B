function[mvpaout_add] = parse_mvpa_results_study_add(args, ph)

%---------------------------------------------------------------------------
%*************** parse the results
%---------------------------------------------------------------------------
% timecourse: baseline = other categories except "rest" class

xph = args.xphase;
fprintf('\n(+) parse the results: %s\n', args.phase_name{xph});

%*************** the study volume was not shifted but trimmed: 
%*************** args.index{xph}.matrix was not shifted/ only trimmed
%*************** args.regs{xph}.regressor was shifted for running classification
%*************** non-shifting for timecourse
xindex  = args.index{xph};% param index from study
xmatrix = xindex.matrix; %trimmed/non-shifted regressor
xheader = xindex.header;

%% ============= GRAB ACTIVATIONS & PERFORMANCE

out_acts           = ph.results.iterations.acts;

%% ============= UNPACK PARAMETERS
xparam             = xindex.param;
xcond_name         = xparam.conds_names;
n_condition        = length(xcond_name);
n_runs             = xparam.n_runs;
n_trials           = xparam.n_trials;
n_category         = xparam.n_category;
n_subcategory      = xparam.n_subcategory;
n_trs              = args.tc_tr;
dur_sync           = args.dur_sync;%next trial 

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%*************** TIMECOURSE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ============= SETUP STRUCTURE
%*************** working memory contents
if strcmp(args.level, 'category')
    it_sames = 1:2;
elseif strcmp(args.level, 'subcategory')
    it_sames = 1:3;
end

for xsame = it_sames
    for xcond = 1:n_condition
        for xtr = 1:n_trs
            decode.timecourse.sameseq{xsame}.condition{xcond}.evidence{xtr}      = [];
        end
        
        %*************** sync trs
        for xtr = 1:dur_sync    
            decode.timecourse.sameseq{xsame}.sync.condition{xcond}.evidence{xtr} = [];
        end
    end
end

%% ============= N=N+1 COLLECTING EVIDENCE: PER CONDITION
%*************** working memory contents
% timecourse: un-shifted tr: -args.shift_TRs
% condition: 1_maintain,2_replace_category,3_replace_item,4_target_suppress,5_global_clear
% xtr: 12 presentation + 5:9 fixation = 21 max

% same category/subcategory for N+1 trial: collecting only target evidence
% decode.timecourse.sameseq{xsame}.condition{xcond}.evidence{xtr}
% decode.timecourse.sameseq{xsame}.sync.condition{xcond}.evidence{xtr}
% xsame: category: 1_same, 2_differ /subcategory: 1_same, 2_differ, 3_related

% 1_maintain/2_replace(category)/3_replace(subcategory)/4_suppress/5_clear

for xrun = 1:n_runs
    
    fprintf('run: %d | trial: ', xrun)
    
    for xtrial = 1:n_trials(xrun)-1
        
        fprintf('%s.', num2str(xtrial))
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        %*************** N trial
        tunit      = getDATA(xmatrix', xheader, ...
            {'run','trial','presentation'}, {xrun, xtrial, 1});
        xcate      = unique(xmatrix(findCol(xheader, {'category'}), tunit));
        xsubcate   = unique(xmatrix(findCol(xheader, {'subcategory'}), tunit));
        
        %*************** target: subtarget class
        xsubtarg    = xsubcate + (n_subcategory * (xcate-1));
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        %*************** N+1 trial
        ttunit     = getDATA(xmatrix', xheader, ...
            {'run','trial','presentation'}, {xrun, xtrial+1, 1});
        xxcate     = unique(xmatrix(findCol(xheader, {'category'}), ttunit));
        xxsubcate  = unique(xmatrix(findCol(xheader, {'subcategory'}), ttunit));
        
        %*************** from N trial
        xcond      = unique(xmatrix(findCol(xheader, {'condition'}), tunit));
        
        %*************** timecourse
        xunit     = find(getDATA(xmatrix', xheader, {'run', 'trial'}, {xrun, xtrial}));
        xunit_tc  = xunit(1):(xunit(1) + n_trs - 1);
        xspike    = ~(xmatrix(findCol(xheader, {'spike'}), xunit_tc));
        
        if sum(xspike)~=0
            for xtr = 1:length(xunit_tc)
                if xspike(xtr)
                    
                    act_matrix = out_acts(:, xunit_tc(xtr));
                    
                    %*************** same category for N+1 trial
                    % 1_category, 2_subcategory
                    % 1_same, 2_differ, 3_related
                    
                    if strcmp(args.level, 'category')
                        
                        xevidence  = act_matrix(xcate);
                        
                        if (xcate == xxcate)
                            xsame = 1;
                        else
                            xsame = 2;
                        end
                        
                    elseif strcmp(args.level, 'subcategory')
                        
                        xevidence  = act_matrix(xsubtarg);
                        
                        if (xcate == xxcate) && (xsubcate == xxsubcate)
                            xsame = 1;
                        elseif (xcate == xxcate) && (xsubcate ~= xxsubcate)
                            xsame = 3;
                        elseif (xcate ~= xxcate)
                            xsame = 2;
                        end
                    end
                    
                    %*************** target
                    decode.timecourse.sameseq{xsame}.condition{xcond}.evidence{xtr} = ...
                        horzcat(decode.timecourse.sameseq{xsame}.condition{xcond}.evidence{xtr}, ...
                        xevidence);
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%
                    %*************** sync trs
                    if xtr > length(xunit)
                        it_tr = xtr - length(xunit);
                        if it_tr <= dur_sync
                            
                            %*************** target
                            decode.timecourse.sameseq{xsame}.sync.condition{xcond}.evidence{it_tr} = ...
                                horzcat(decode.timecourse.sameseq{xsame}.sync.condition{xcond}.evidence{it_tr}, ...
                                xevidence);
                        end
                    end
                end
            end
        end
    end
    
    fprintf('\n')
end

%% ============= NEW ITEM FOR SUPPRESS/CLEAR
%*************** working memory contents
% cond: 4, 5
% decode.control.condition{xit}.new{xnew}.evidence{xtr}
% 1_maintain/2_replace(category)/3_replace(subcategory)/4_suppress/5_clear

if strcmp(args.level, 'category')
    %*************** number of condition trials
    clear new_cate_array
    it_conds = [4, 5];
    cate_array = 1:n_category;
    
    for xit = 1:length(it_conds)
        xcond = it_conds(xit);
        for xcate = 1:n_category
            clear cond_trials t_new_array
            
            cond_trials = unique(getDATA(xmatrix', xheader, ...
                {'condition', 'category'}, {xcond, xcate}, findCol(xheader, {'trial'})));
            it_news = cate_array(~ismember(cate_array, xcate));
            t_new_array = [];
            
            xhalf(1) = round(length(cond_trials)/2);
            xhalf(2) = length(cond_trials) - xhalf(1);
            
            for yit = 1:length(it_news)
                xnew = it_news(yit);
                t_new_array = horzcat(t_new_array, ones(1, xhalf(yit)) * xnew);
            end
            new_cate_array{xit, xcate}(1,:) = cond_trials;% trials
            new_cate_array{xit, xcate}(2,:) = shuffle(t_new_array);
        end
    end
    
    for xit = 1:length(it_conds)
        for xtr = 1:n_trs
            decode.control.condition{xit}.new_evidence{xtr} = [];
        end
    end
    
    %*************** collect evidence
    for xit = 1:length(it_conds)
        xcond = it_conds(xit);
        for xcate = 1:n_category
            for xtrial = 1:size(new_cate_array{xit, xcate}, 2)
                xnew_cate = new_cate_array{xit, xcate}(2, xtrial);
                
                %*************** timecourse
                xunit = find(getDATA(xmatrix', xheader, ...
                    {'condition','category','trial'}, ...
                    {xcond, xcate, new_cate_array{xit, xcate}(1, xtrial)}));
                xunit_tc = xunit(1):(xunit(1) + n_trs - 1);
                xspike   = ~(xmatrix(findCol(xheader, {'spike'}), xunit_tc));
                
                if sum(xspike)~=0
                    for xtr = 1:length(xunit_tc)
                        if xspike(xtr)
                            %*************** new category
                            act_matrix = out_acts(:, xunit_tc(xtr));
                            xnew_evidence  = act_matrix(xnew_cate);
                            
                            decode.control.condition{xit}.new_evidence{xtr} = ...
                                horzcat(decode.control.condition{xit}.new_evidence{xtr}, ...
                                xnew_evidence);
                        end
                    end
                end
            end
        end
    end
end

%% ============= SAVE MVPAOUT
mvpaout_add.decode = decode;

end
