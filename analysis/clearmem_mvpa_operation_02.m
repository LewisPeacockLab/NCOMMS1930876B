function[] = clearmem_mvpa_operation_02(args, dirs)

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

%*************** ph4. base filename
ph4.basename = sprintf('%s_decoding_setup', ph3.basename);

%*************** load phase 6.
fname = sprintf('%s/ph4_%s.mat', dirs.mvpa.scratch{xph}, ph4.basename);
load(fname);%,'ph4'

%*************** merge subj. structure
subj = ph4.subj; nm = ph4.nm;
summarize(subj)

%% ============= 05: MVPA CLASSIFICATION & PENALTY CHECK
%*************** train a classifier + do cross-validation.
fprintf('\n#####################################################################\n');
fprintf('* STEP_04: classification\n');

%*************** basename phase 5.
class_basename   = sprintf('classified_%s_%s', ph4.basename, args.classifier);

%% *************** penalty check for individualized penalty
if ~(args.fixed_penalty{xph})
    
    penalty_basename = sprintf('classified_%s_%s', ph3.basename, args.classifier);
    
    if args.penalty_check{xph}
        
        pen_check.penalty     = [];
        pen_check.performance = [];
        
        %% *************** classification
        if args.penalty_step{xph}, it_steps = 1:2;%1_broad/2_narrow
        else it_steps = 1; end
        
        for xstep = it_steps
            
            if xstep==1
                it_penalty = args.penalty;
                all_pen    = length(args.penalty);
            else
                it_penalty = args.narrow_penalty;
                all_pen    = length(args.penalty) + length(args.narrow_penalty);
            end
            
            for i = 1:numel(it_penalty)
                %*************** definde penalty + classification
                if xstep==1
                    args.xpenalty = it_penalty(i); it = i;
                else
                    args.xpenalty = it_penalty(i); it = i + numel(args.penalty);
                end
                
                fprintf('\n(+) subject:%s classifier: %s penalty: %s (%s out of %s)\n', ...
                    args.subject_id, args.classifier, num2str(args.xpenalty), ...
                    num2str(it), num2str(all_pen));
                fprintf('... classification %s phase \n', args.phase_name{xph});
                
                [ph5.args, ~, ph5.results] = mvpa_ph04_classification(args, subj, nm);
                
                if args.four_oper_regress
                    [mvpaout] = parse_mvpa_results_operation_four(args, ph5, dirs, 0);
                else
                    [mvpaout] = parse_mvpa_results_operation(args, ph5, dirs, 0);
                end
                
                %*************** write performance depending on penalty
                for xiter = 1:size(ph5.results.iterations, 2)
                    pen_check.xpenalty(it).iteration(xiter) = ...
                        mvpaout.decode.classifier.iteration{xiter}.accuracy;
                end
                
                xtotal_perf = mean(pen_check.xpenalty(it).iteration);
                
                pen_check.penalty     = [pen_check.penalty it_penalty(i)];
                pen_check.performance = [pen_check.performance xtotal_perf];
                
            end
            
            fprintf('\n');
            
            %% *************** define step 2 penalty
            if args.penalty_step{xph} && (xstep==1)
                [xacc, whichmax] = max(pen_check.performance);
                
                fprintf('... step1: max acc: %4.4f\n', xacc);
                
                if whichmax == 1
                    % If the best penalty is the first broad search penalty (zero), only
                    % search between that and the upper bound (there is no lower bound):
                    xmax    = pen_check.penalty(whichmax);
                    xmax_up = pen_check.penalty(whichmax+1);
                    args.narrow_penalty = xmax:(xmax_up/2)/(args.n_pen-1):xmax_up/2;
                    
                elseif whichmax == length(args.penalty)
                    % If the best penalty was the last broad search value, try
                    % NPENALTIES between the normal lower bound and an upper bound
                    % beyond the broad search penalties (i.e. for 10,000 the upper
                    % bound would be 50,000, from 10,000*10 / 2):
                    xmax    = pen_check.penalty(whichmax);
                    args.narrow_penalty = xmax/2:((xmax*10)/2-(xmax/2))/(args.n_pen-1):(xmax*10)/2;
                    
                else
                    % For when the best penalty is not the first or last:
                    xmax    = pen_check.penalty(whichmax);
                    xmax_up = pen_check.penalty(whichmax+1);
                    args.narrow_penalty = xmax/2:(xmax_up/2-xmax/2)/(args.n_pen-1):xmax_up/2;
                end
            end
        end
        
        %% ============= save penalty_check
        %*************** save mat
        save(sprintf('%s/penalty_check_%s.mat', dirs.mvpa.output{xph}, penalty_basename),...
            'pen_check','-v7.3');
        
        [xacc, whichmax] = max(pen_check.performance);
        max_penalty      = pen_check.penalty(whichmax);
        
        %*************** open penalty performance check file
        fpenalty_check = fopen(sprintf('%s/penalty_check_%s.txt', dirs.mvpa.output{xph}, penalty_basename),'w+');
        
        %*************** 1 row: names penalty
        fprintf(fpenalty_check,'penalty: ');
        
        for xx=1:length(pen_check.penalty)
            fprintf(fpenalty_check,'%s\t', num2str(pen_check.penalty(xx)));
        end
        
        for i=1:numel(unique(args.regs{xph}.selectors))
            fprintf(fpenalty_check,'n_voxels_%d\t', i);
        end
        
        fprintf(fpenalty_check,'max_penalty max_acc\n');
        
        %*************** 2 row: accuracy
        fprintf(fpenalty_check,'accuracy: ');
        
        for xx=1:length(pen_check.performance)
            fprintf(fpenalty_check,'%1.4f\t', pen_check.performance(xx));
        end
        
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
        
        fprintf(fpenalty_check,'%s %1.4f\t', num2str(max_penalty), xacc);
        fclose(fpenalty_check);
        
        fprintf('... max_penalty: %s: acc: %1.4f\n', num2str(max_penalty), xacc);
        
    else
        %% ============= load penalty performance check file
        load(sprintf('%s/penalty_check_%s.mat', dirs.mvpa.output{xph}, penalty_basename));%'pen_check'
        
        if (args.reset_regs{xph}) || (args.reset_6tr), xpen = pen_check; %#ok<*NODEF>
        else                                           xpen = penalty_check; end
        
        [xacc, whichmax] = max(xpen.performance);
        max_penalty      = xpen.penalty(whichmax);
        fprintf('... max_penalty: %s: acc: %1.4f\n', num2str(max_penalty), xacc);
        
    end
    
    args.xpenalty = max_penalty;
end

%% ============= 05: MVPA CLASSIFICATION for MAX PENALTY
%*************** definde penalty + classification
fprintf('\n(+) subject:%s classifier: %s, max penalty: %s\n', ...
    args.subject_id, args.classifier, num2str(args.xpenalty));
fprintf('... classification %s phase \n', args.phase_name{xph});

[ph5.args, ~, ph5.results] = mvpa_ph04_classification(args, subj, nm);

%*************** basename phase 5.
ph5.basename  = sprintf('%s_penalty%s', class_basename, num2str(args.xpenalty));

%*************** selected volxels
for xiter=1:length(unique(args.regs{xph}.selectors))
    ph5.param.n_total_voxels(xiter,1) = ph4.subj.masks{1}.nvox;
    ph5.param.n_voxels(xiter,1)       = ph4.subj.masks{xiter+1}.nvox;
    ph5.param.n_percent(xiter,1)      = (ph4.subj.masks{xiter+1}.nvox/ph4.subj.masks{1}.nvox)*100;
end

%*************** save phase 5.
fname = sprintf('%s/%s.mat', dirs.mvpa.output{xph}, ph5.basename);
save(fname,'ph5','-v7.3');

fprintf('\n... saved results_penalty%s of %s: %s\n', ...
    num2str(args.xpenalty), args.subject_id, fname);
    
end