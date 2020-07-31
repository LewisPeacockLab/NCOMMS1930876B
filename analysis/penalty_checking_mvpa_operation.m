function[] = penalty_checking_mvpa_operation(args, dirs)

xph = args.xphase;

%% ============= Load classification_01 file
%*************** ph1. base filename
ph1.basename = sprintf('%s_%s_%s', args.phase_name{xph}, args.mask_name, args.epi_name); 

%*************** ph2. base filename
if strcmp(args.regress_type, 'shift')
    ph2.basename = sprintf('%s_%s_tr%s_blk_%s',...
                   ph1.basename, args.regress_type, ...
                   num2str(args.shift_TRs), args.rest);
elseif strcmp(args.regress_type, 'beta')
    ph2.basename = sprintf('%s_%s', ...
                   ph1.basename, args.regress_type);
end

%*************** ph3. base filename
if args.featVox
    ph3.basename = sprintf('%s_featsel_%svox', ph2.basename, num2str(args.fsVosNum));
else
    ph3.basename = sprintf('%s_featsel_thresh%s', ph2.basename, num2str(args.featSelThresh));
end

%*************** basename.
penalty_basename = sprintf('classified_%s_%s', ph3.basename, args.classifier);
    
    %% *************** classification
    for i = 1:numel(args.penalty)
        %*************** definde penalty + classification
        args.xpenalty = args.penalty(i);
        
        fprintf('\n(+) subject:%s classifier: %s penalty: %s (%s out of %s)\n', ...
            args.subject_id, args.classifier, num2str(args.xpenalty), ...
            num2str(i), num2str(numel(args.penalty)));
        fprintf('... classification %s phase \n', args.phase_name{xph});
        
        [ph5.args, ~, ph5.results] = mvpa_ph04_classification(args, subj, nm);
        
        [mvpaout] = parse_mvpa_results_operation(args, ph5, dirs, 1);
        
        %*************** write performance depending on penalty
        for xiter = 1:size(ph4.results.iterations, 2)
            penalty_check.xpenalty(i).iteration(xiter) = ...
                mvpaout.decode.classifier.iteration{xiter}.accuracy;
        end
        
        xtotal_perf = mean(penalty_check.xpenalty(i).iteration);
        
        fprintf(fpenalty_check,'%4.4f\t', xtotal_perf);
        
        penalty_check.performance = [penalty_check.performance xtotal_perf];
        
    end
    
    %% ============= save penalty_check
    save(sprintf('%s/penalty_check_%s.mat', dirs.mvpa.output{xph}, penalty_basename),...
        'penalty_check');
    
    %*************** write selected n_voxels
    if args.cross_valid % localizer
        if args.featSel
            for i=1:numel(unique(args.regs{xph}.selectors))
                n_selected_vox(i) = subj.masks{i+1}.nvox; %#ok<*AGROW>
                
                fprintf(fpenalty_check,'%s\t', num2str(n_selected_vox(i)));
            end
        else
            n_selected_vox(1) = subj.masks{1}.nvox; %#ok<*AGROW>
            
            fprintf(fpenalty_check,'%s\t', num2str(n_selected_vox(1)));
        end
    else% decoding
        t_mask = get_object(subj,'mask',nm.mask_selected{1}{1});
        n_selected_vox = t_mask.nvox;
        
        fprintf(fpenalty_check,'%s\t', num2str(n_selected_vox));
    end
    
    max_penalty = penalty_check.penalty(penalty_check.performance==max(penalty_check.performance));
    fprintf(fpenalty_check,'%s\n', num2str(min(max_penalty)));
    fprintf('... max_penalty: %s\n', num2str(min(max_penalty)));
    

    %% ============= load penalty performance check file
    load(sprintf('%s/penalty_check_%s.mat', dirs.mvpa.output{xph}, penalty_basename));%'penalty_check'
    
    max_penalty = penalty_check.penalty(penalty_check.performance==max(penalty_check.performance)); %#ok<*NODEF>
    fprintf('%s\n', num2str(min(max_penalty)));
    


%% ============= 05: MVPA CLASSIFICATION for MAX PENALTY
%*************** definde penalty + classification
args.xpenalty = min(max_penalty);

fprintf('\n(+) subject:%s classifier: %s, max penalty: %s\n', ...
    args.subject_id, args.classifier, num2str(args.xpenalty));
fprintf('... classification %s phase \n', args.phase_name{xph});

[ph5.args, ~, ph5.results] = mvpa_ph04_classification(args, subj, nm);

%*************** basename phase 4.
ph5.basename  = sprintf('%s_penalty%s', class_basename, num2str(args.xpenalty));

%*************** save phase 4.
fname = sprintf('%s/%s.mat', dirs.mvpa.output{xph}, ph5.basename);
save(fname,'ph5','-v7.3');

fprintf('\n... saved results_penalty%s of %s: %s\n', ...
    num2str(args.xpenalty), args.subject_id, fname);
    
end