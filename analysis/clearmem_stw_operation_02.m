function[] = clearmem_stw_operation_02(args, dirs)

xph = args.xphase;

%% ============= LOGGING
fprintf('(+) using scratch dir: %s\n', dirs.mvpa.scratch{xph});
%*************** turn on diary to capture analysis output
diary off;
diary(sprintf('%s/diary.txt', dirs.mvpa.scratch{xph}));
fprintf('running code: %s\n', mfilename)
fprintf('#####################################################################\n\n');
disp(args);
fprintf('#####################################################################\n');

%% ============= SETUP REGRESSORS
%*************** define regressors + selectors
%*************** regressors: shift or beta
xheader        = args.index{xph}.header;%from study
xmatrix        = args.index{xph}.matrix;
xparam         = args.index{xph}.param;

%*************** unpack parameters
condition_name = xparam.conds_names;
n_condition    = length(condition_name);
n_trials       = xparam.n_trials;
% n_runs         = xparam.n_runs;

%*************** sliding time window param
it_wins = 1:(args.tc_tr_disp+1);

%% ============= PARSE THE RESULTS

for xcond = 1:n_condition
    for xtarg = 1:n_condition
        for xtr = it_wins
            tc_stw.operation{xcond}.decoded_operation{xtarg}.tr{xtr} = [];
            tc_stw.filter.operation{xcond}.decoded_operation{xtarg}.tr{xtr} = [];
        end
    end
end

%% ************* TIMECOURSE: collecting evidence
for xwin = it_wins
    
    %% ============= LOAD MVPA RESULTS
    clear ph4
    fprintf('\n#####################################################################\n');
    fprintf('(+) window %s: ', num2str(xwin));
    fprintf('\n#####################################################################\n\n');
    
    %*************** ph1. base filename
    ph1.basename = sprintf('%s_%s_%s', args.phase_name{xph}, args.mask_name, args.epi_name);
    
    %*************** ph2. base filename
    ph2.basename = sprintf('win%s_%s_%s_tr%s_blk_%s',...
        num2str(xwin), ph1.basename, args.regress_type, num2str(args.shift_TRs), args.rest);
    
    %*************** ph3. base filename
    ph3.basename = sprintf('%s_featsel_thresh%s', ph2.basename, num2str(args.featSelThresh));
    
    %*************** basename phase 4.
    ph4.basename = sprintf('%s_penalty%s', ph3.basename, num2str(args.xpenalty));
    
    %*************** save phase 4.
%     fname = sprintf('%s/%s.mat', dirs.mvpa.output{xph}, ph4.basename);
    fname = sprintf('%s/%s_result.mat', dirs.mvpa.output{xph}, ph4.basename);
    
    fprintf('\n... loading results_penalty%s of %s: %s\n', ...
        num2str(args.xpenalty), args.subject_id, fname);
    load(fname);%'ph4': ph4.args, ph4.subj, ph4.results
    
    n_iters = length(xresults.iterations);
    
    %% ============= CORRECT EVIDENCE
    %*************** regressor
    clear test_idx out_act
    for xrun = 1:n_iters
        test_idx{xrun} = xresults.iterations(xrun).test_idx; %#ok<*AGROW>
        out_act{xrun}  = xresults.iterations(xrun).acts;
        
        for xtrial = 1:n_trials(xrun)
            tunit  = find(getDATA(xmatrix', xheader, ...
                {'run','trial'}, {xrun, xtrial}));
            
            xcond  = xmatrix(findCol(xheader, 'condition'), tunit(1));
            xunit  = tunit(1) + (xwin-1);
            xacts  = out_act{xrun}(:, (test_idx{xrun}==xunit));
            
            if ~(isempty(xacts))
                for xtarg = 1:n_condition
                    tc_stw.operation{xcond}.decoded_operation{xtarg}.tr{xwin} = ...
                        horzcat(tc_stw.operation{xcond}.decoded_operation{xtarg}.tr{xwin}, ...
                        xacts(xtarg, 1));
                end
            end
            
            if xtrial > 1
                %*************** filtered
                punit  = find(getDATA(xmatrix', xheader, ...
                    {'run','trial'}, {xrun, xtrial-1}));
                
                pcond  = xmatrix(findCol(xheader, 'condition'), punit(1));
                
                if xcond~=pcond
                    if ~(isempty(xacts))
                        for xtarg = 1:n_condition
                            tc_stw.filter.operation{xcond}.decoded_operation{xtarg}.tr{xwin} = ...
                                horzcat(tc_stw.filter.operation{xcond}.decoded_operation{xtarg}.tr{xwin}, ...
                                xacts(xtarg, 1));
                        end
                    end
                end
            end
        end
    end
end

%% ============= SAVE PAMVPA
tc_stw.info.test_idx = test_idx;
tc_stw.info.out_act  = out_act;
tc_stw.info.matrix   = xmatrix;
tc_stw.info.header   = xheader;

mvpaout_name = sprintf('%s/mvpaout_%s.mat', dirs.mvpa.parse{xph}, args.analysis_basename);
save(mvpaout_name, 'tc_stw', '-v7.3');

%%
% for xcond = 1:n_condition
%     clear xmean xstd
%     figure;
%     
%     for xtarg = 1:n_condition
%         for xwin = it_wins
%             xmean(xwin, xtarg) = mean(tc_stw.operation{xcond}.decoded_operation{xtarg}.tr{xwin});
%             xstd(xwin, xtarg)  = std(tc_stw.operation{xcond}.decoded_operation{xtarg}.tr{xwin});
%         end
%     end
%     
%     for xtarg = 1:n_condition
%         if xtarg==xcond
%             plot(xmean(:,xtarg), 'LineWidth', 2); 
%         else
%             plot(xmean(:,xtarg)); 
%         end
%         hold on;
%     end
%     legend(condition_name)
%     legend(condition_name,'AutoUpdate','off')
%     
%     for xtarg = 1:n_condition
%         errorbar(it_wins, xmean(:,xtarg), xstd(:,xtarg));
%     end
% end


% %%
% figure;
% 
% for xcond = 1:n_condition
%     for xtarg = 1:n_condition
%         
%         subplot(5,5,xtarg + n_condition * (xcond-1))
%         
%         hist(tc_stw.filter.operation{xcond}.decoded_operation{xtarg}.tr{1});
%     end
% end

end
