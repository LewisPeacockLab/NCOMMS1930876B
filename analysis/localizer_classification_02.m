function[] = localizer_classification_02(args, dirs)
% step 2: classification (with all penalties if args.penalty_check = 1)
%         if args.penalty_check = 0, load maximum penalty

%% ============= UNPACK ARGS.
xph = args.xphase;

%% ============= LOGGING
diary off;
fprintf('(+) using scratch dir: %s\n', dirs.mvpa.scratch{xph});
%*************** turn on diary to capture analysis output
diary(sprintf('%s/diary.txt', dirs.mvpa.scratch{xph}));
fprintf('running code: %s\n', mfilename)
fprintf('#####################################################################\n\n');
disp(args);
fprintf('#####################################################################\n');

%% ============= Load classification_01 file
args.regress_type = args.train_regress_type;

%*************** base filename
ph1.basename = sprintf('%s_%s_zscored_%s', args.phase_name{xph}, args.mask_name, args.epi_name); 
ph2.basename = sprintf('%s_%s_%s%dtr_blk_%s',...
    ph1.basename, args.level, args.regress_type, ...
    args.shift_TRs, args.rest);
ph3.basename = sprintf('%s_featsel_thresh%s', ph2.basename, num2str(args.featSelThresh));

%*************** load phase 3.
fname = sprintf('%s/ph3_%s.mat', dirs.mvpa.scratch{xph}, ph3.basename);
load(fname);%,'ph3'

%*************** merge subj. structure
subj = ph3.subj; nm = ph3.nm;
summarize(subj)

%% ============= 04: MVPA CLASSIFICATION & PENALTY CHECK
fprintf('\n#####################################################################\n');
fprintf('* STEP_04: classification\n');

%*************** basename phase 4.
class_basename   = sprintf('classified_%s_%s', ph3.basename, args.classifier);

%% *************** penalty check for individualized penalty
if ~(args.fixed_penalty{xph})
    
    penalty_basename = sprintf('classified_%s_%s', ph3.basename, args.classifier);
    
    if args.penalty_check{xph}
        
        pen_check.penalty     = [];
        pen_check.performance = [];
        
        %% *************** classification
        for xstep = 1:2%1_broad/2_narrow
            
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
                
                [ph4.args, ph4.subj, ph4.results] = mvpa_ph04_classification(args, subj, nm);
                
                %*************** write performance depending on penalty
                for xiter = 1:size(ph4.results.iterations, 2)
                    pen_check.xpenalty(i).iteration(xiter) = ph4.results.iterations(xiter).perf;
                end
                
                pen_check.penalty     = [pen_check.penalty it_penalty(i)];
                pen_check.performance = [pen_check.performance ph4.results.total_perf];
                
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
        %% ============= loac penalty performance check file
        load(sprintf('%s/penalty_check_%s.mat', dirs.mvpa.output{xph}, penalty_basename));%'penalty_check'
        
        [xacc, whichmax] = max(pen_check.performance); %#ok<*NODEF>
        max_penalty      = pen_check.penalty(whichmax);
        fprintf('... max_penalty: %s: acc: %1.4f\n', num2str(max_penalty), xacc);
        
    end
    
    args.xpenalty = max_penalty;
end

%% ============= 04: MVPA CLASSIFICATION for MAX PENALTY
%*************** definde penalty + classification
fprintf('\n(+) subject:%s classifier: %s, max penalty: %s\n', ...
    args.subject_id, args.classifier, num2str(args.xpenalty));
fprintf('... classification %s phase \n', args.phase_name{xph});

[ph4.args, ph4.subj, ph4.results] = mvpa_ph04_classification(args, subj, nm);

%*************** basename phase 4.
ph4.basename  = sprintf('%s_penalty%s', class_basename, num2str(args.xpenalty));

%*************** save phase 4.
fname = sprintf('%s/%s.mat', dirs.mvpa.output{xph}, ph4.basename);
save(fname,'ph4','-v7.3');

fprintf('\n... saved results_penalty%s of %s: %s\n', ...
    num2str(args.xpenalty), args.subject_id, fname);

end