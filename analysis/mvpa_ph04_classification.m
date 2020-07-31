function[args, subj, results] = mvpa_ph04_classification(args, subj, nm)
    
    %---------------------------------------------------------------------------
    %*************** train a classifier + do cross-validatoin.
    %---------------------------------------------------------------------------    
    
    xph = args.xphase;
    
    %% ============= DETERMINE CLASSIFIER
    switch(args.classifier)      
        case 'bp'
            class_args.train_funct_name = 'train_bp';
            class_args.test_funct_name  = 'test_bp';
            class_args.nHidden          = 0;
            nIterations                 = args.repetitions;
            
        case 'L2logreg'            
            class_args.train_funct_name = 'train_L2_RLR';
            class_args.test_funct_name  = 'test_L2_RLR';
            class_args.lambda           = args.xpenalty;
            nIterations                 = 1;
            
        case 'logreg'
            class_args.train_funct_name = 'train_logreg';
            class_args.test_funct_name  = 'test_logreg';
            class_args.penalty          = args.xpenalty;
            nIterations                 = 1;
            
        case 'ridge'
            class_args.train_funct_name = 'train_ridge';
            class_args.test_funct_name  = 'test_ridge';
            class_args.penalty          = args.xpenalty;
            nIterations                 = 1;
            
        otherwise
            disp('(-) unknown classifier type!');
            return            
    end
    
    %% ============= CROSS-VALIDATION CLASSIFICATION
    %*************** getting rid of 1-of-n warnings using propval
    cv_args.perfmet_args.ignore_1ofn = true; 
    
    %*************** classification iteration
    for i=1:nIterations
        if strcmp(args.regress_type, 'shift')
            nm_regs = nm.conds_sh{xph};
        elseif strcmp(args.regress_type, 'beta')
            nm_regs = nm.conds{xph};
        end
        
        if args.cross_valid
            mask_name = sprintf('%s_patterns_z_thresh0.05', args.phase_name{xph});
%             subj.masks{2}.group_name;
        else
            mask_name = subj.masks{3}.name;
        end
        
        if args.cross_valid
            if args.cross_decoding
                fprintf('\n... classification: cross-validation (leave 1-out): testing on all timepoints\n');
                
                [subj, results] = cross_validation(subj, nm.pats_z{xph}, nm_regs,...
                    nm.decoding_runs, mask_name, class_args, cv_args);
                
            else
                fprintf('\n... classification: cross-validation (leave 1-out)\n');
                
                [subj, results] = cross_validation(subj, nm.pats_z{xph}, nm_regs,...
                    nm.run_xvalid{xph}, mask_name, class_args, cv_args);
            end
        else% decoding 
            fprintf('\n... decoding %d of %d repetitions)\n', i, nIterations);
            
            [subj, results] = cross_validation(subj, nm.combined_pats, nm.combined_regs,...
                nm.combined_runs, mask_name, class_args, cv_args);
        end
    end
    
    args.date.end = sprintf('%s_%s', datetime, args.hostname);
    
    fprintf('... end of classification: %s\n', args.date.end);

end