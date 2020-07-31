function[mvpaout] = parse_mvpa_results_localizer(args, ph, dirs)

    %---------------------------------------------------------------------------
    %*************** parse the results
    %---------------------------------------------------------------------------
    
    xph = args.xphase;
    fprintf('\n(+) parse the results: %s\n', args.phase_name{xph});
    
    %% ============= UNPACK ARGS & PARAMS
    mvpaout.total_perf = ph.results.total_perf;
    n_validation       = numel(ph.results.iterations);
    param              = args.index{xph}.param;
    
    n_subcategory      = param.n_subcategory;
    n_category         = param.n_category;
    it_categories      = 1:n_category;
    it_subcategories   = 1:n_subcategory;
    
    if strcmp(args.level, 'subcategory') && args.class_selecting
        it_categories = args.selected_category;
       
        if args.subclass_selecting%only when selected category is one
            it_subcategories = args.selected_subcategory;
        end
        
    elseif strcmp(args.level, 'category') && args.class_selecting
        it_categories = args.selected_category;
    end
    
    %*************** output basename
    base_name = args.analysis_basename;

    %% ============= PARSE THE RESULTS
    for xiter=1:n_validation
        %*************** pfmet: based on shifted regressor / (startTR:end)
        pfmet.guesses  = ph.results.iterations(xiter).perfmet.guesses;
        pfmet.acts     = ph.results.iterations(xiter).acts;
        pfmet.desireds = ph.results.iterations(xiter).perfmet.desireds;
        
        %*************** correctness
        xcorr          = (pfmet.desireds) - (ph.results.iterations(xiter).perfmet.guesses);
        pfmet.corr     = zeros(1, length(pfmet.desireds));
        pfmet.corr(xcorr==0) = 1;

        n_target       = max(pfmet.desireds);

        for xcate=1:n_target

            xunit     = pfmet.desireds == xcate;
            t_guess   = pfmet.guesses(xunit);
            t_acts    = pfmet.acts(:,xunit);
            t_correct = pfmet.corr(xunit);

            for xguess=1:n_target
                mvpaout.decoding{xcate}(xiter,xguess) = numel(find(t_guess == xguess));
                mvpaout.evidence{xcate}(xiter,xguess) = mean(t_acts(xguess,:));    
            end
            
            mvpaout.accuracy{xcate}(xiter) = mean(t_correct);
        end
    end

    %% ============= ANALYSIS PARAMETERS
    %*************** selected volxels
    for xiter=1:n_validation
        n_total_voxels(xiter,1) = ph.subj.masks{1}.nvox;
        n_voxels(xiter,1)       = ph.subj.masks{xiter+1}.nvox;
        n_percent(xiter,1)      = (ph.subj.masks{xiter+1}.nvox/ph.subj.masks{1}.nvox)*100;
    end

    %% ============= WRITE TABLE TO CSV FILES FOR PYTHON
    if strcmp(args.level, 'category')
        %*************** create tables
        for xcate = 1:length(it_categories)
            itcate = it_categories(xcate);
            xcate_name{xcate} = param.category_name{itcate};
        end
    elseif strcmp(args.level, 'subcategory')
        %*************** create tables
        for xcate = 1:length(it_categories)
            itcate = it_categories(xcate);
            for xsubcate = 1:length(it_subcategories)
                itsubcate = it_subcategories(xsubcate);
                xcate_name{xsubcate + length(it_subcategories)*(xcate-1)} = ...
                    param.subcategory_name{itcate}{itsubcate};
            end
        end
    end
    
    %*************** rest
    if strcmp(args.rest, 'rest')
        xcate_name{n_target} = 'rest';
    end
    
    %*************** create tables
    for xcate = 1:n_target
        clear n_total tacc tevi
        
        guess_names(xcate)  = strcat(sprintf('g%s_', num2str(xcate)), xcate_name(xcate));
        target_names(xcate) = strcat(sprintf('t%s_', num2str(xcate)), xcate_name(xcate));
    
        for xiter = 1:n_validation
            n_total            = sum(mvpaout.decoding{xcate}(xiter,:));
            
            tacc(xiter,:)      = mvpaout.decoding{xcate}(xiter,:)/n_total;
            tevi(xiter,:)      = mvpaout.evidence{xcate}(xiter,:);
        end
        
        xaccuracy(xcate,:)     = mean(tacc); %#ok<*AGROW>
        xevidence(xcate,:)     = mean(tevi);
        targ_accuracy{xcate+1} = xaccuracy(xcate, xcate);%target accuracy
        targ_evidence{xcate+1} = xevidence(xcate, xcate);%target accuracy
    end
    
    %*************** create tables
    table_accuracy = array2table(xaccuracy,...
        'RowNames', target_names, 'VariableNames', guess_names);
    
    table_evidence = array2table(xevidence,...
        'RowNames', target_names, 'VariableNames', guess_names);
    
    targ_accuracy{1}     = args.subject_id;
    targ_accuracy{end+1} = mean([targ_accuracy{2:length(xcate_name)+1}]);
    targ_accuracy{end+1} = ph.args.xpenalty;
    targ_accuracy{end+1} = n_total_voxels;
    targ_accuracy{end+1} = n_voxels;
    targ_accuracy{end+1} = n_percent;
    
    targ_evidence{1}     = args.subject_id;
    
    table_targ_accuracy  = cell2table(targ_accuracy, 'VariableNames', ...
        horzcat('subject_id', target_names, 'total_accuracy','penalty',...
        'total_voxels', 'selected_voxels', 'percentage'));
    
    table_targ_evidence  = cell2table(targ_evidence, 'VariableNames', ...
        horzcat('subject_id', target_names));
    
    %*************** write tables to csv files
    csv_name = sprintf('%s/table_accuracy_%s.csv', dirs.mvpa.parse{xph}, base_name);
    writetable(table_accuracy, csv_name,'WriteRowNames',true)
    
    csv_name = sprintf('%s/table_evidence_%s.csv', dirs.mvpa.parse{xph}, base_name);
    writetable(table_evidence, csv_name,'WriteRowNames',true)
    
    csv_name = sprintf('%s/table_target_accuracy_%s.csv', dirs.mvpa.parse{xph}, base_name);
    writetable(table_targ_accuracy, csv_name,'WriteRowNames',true)
    
    csv_name = sprintf('%s/table_target_evidence_%s.csv', dirs.mvpa.parse{xph}, base_name);
    writetable(table_targ_evidence, csv_name,'WriteRowNames',true)
    
    %*************** selected volxels  
    table_selected_voxels = table(n_total_voxels, n_voxels, n_percent, ...
        'VariableNames', {'total_voxels', 'selected_voxels', 'percentage'}');
    
    csv_name = sprintf('%s/table_n_voxels_%s.csv', dirs.mvpa.parse{xph}, base_name);
    writetable(table_selected_voxels, csv_name,'WriteRowNames',true)   
    
end
