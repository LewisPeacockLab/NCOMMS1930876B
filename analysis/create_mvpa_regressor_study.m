function[regs] = create_mvpa_regressor_study(args, dirs)
%% =============== LOAD DESIGN INDEX
%***************** trim the matrix -> shift regressor
%***************** non-shifted decoding: shift TRs in the regressor

xph          = 2;
regs.header  = args.index{xph}.header;
regs.matrix  = args.index{xph}.matrix;
param        = args.index{xph}.param;

%***************** unpack parameters
n_category    = param.n_category;
n_subcategory = param.n_subcategory;
n_volumes     = param.n_cat_trim_vols;%size(args.index{xph}.matrix,2);%

category_name    = param.category_name;
subcategory_name = param.subcategory_name;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% GENERATE REGRESSOR FOR MVPA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%================ REGRESSOR/SELECTOR FOR SHIFTED EPI
%***************** selectors: runs
regs.selectors = regs.matrix(findCol(regs.header, {'run'}), :);

if ~strcmp(args.level, 'item')
    %***************** defining regressor names
    regress_name = [];
    if strcmp(args.level, 'category')
        regress_name = category_name(args.selected_category);
        rest_class   = length(args.selected_category) + 1;
        
    elseif strcmp(args.level, 'subcategory')
        for xcate = args.selected_category%n_category
            regress_name = [regress_name subcategory_name{xcate}];
        end %#ok<*AGROW>
        rest_class = (length(args.selected_category) * param.n_subcategory) + 1;
    end
    
    %***************** regressors
    regs.regressor_name   = regress_name;
    regs.regressors       = zeros(length(regress_name), n_volumes);
    regs.regressors_index = zeros(1, n_volumes);
    
    for it_cate = 1:length(args.selected_category)
        xcate = args.selected_category(it_cate);
        for xsubcate = 1:n_subcategory
            xunit = find(getDATA(regs.matrix', regs.header, ...
                {'category','subcategory'}, {xcate, xsubcate})) + args.shift_TRs; % ------------ shifting TRs
            
            if strcmp(args.level, 'category'),
                xregs = it_cate;
            elseif strcmp(args.level, 'subcategory'),
                xregs = xsubcate + (n_subcategory * (it_cate-1));
            end
            
            regs.regressors(xregs, xunit) = 1;
            regs.regressors_index(xunit)  = xregs;
        end
    end
    
    %% =============== REST CATEGORY
    if strcmp(args.rest, 'rest')
        xregs = rest_class;
        regs.regressor_name{xregs} = 'rest';
        regs.regressors(xregs,:)   = zeros(1, n_volumes);
        
        sum_regressors = sum(regs.regressors);
        regs.regressors(xregs, sum_regressors==0) = 1;
        regs.regressors_index(sum_regressors==0)  = xregs;
        
        %% =============== PICK CLASSES=0
        %***************** set 0 for excluded category in rest regressor
        if args.class_selecting
            classes   = 1:n_category;
            xcate_cut = classes(~ismember(classes, args.selected_category));
            
            xregs = rest_class;
            
            for it_cate = xcate_cut
                xunit   = getDATA(regs.matrix', regs.header, {'category'}, {it_cate});
                
                regs.regressors(xregs, xunit) = 0;
                regs.regressors_index(xunit)  = 0;
            end
        end
    end
    
    %% =============== SAVE REGRESSOR
    if args.class_selecting
        save(fullfile(dirs.mvpa.regressor{xph}, ...
            sprintf('study_classification_regressor_%s_%s_%dtr_%s_%s.mat', ...
            sprintf('%d',args.selected_category), args.train_regress_type, ...
            args.shift_TRs, args.level, args.rest)),'regs');
    else
        save(fullfile(dirs.mvpa.regressor{xph}, ...
            sprintf('study_classification_regressor_%s_%dtr_%s_%s.mat', ...
            args.train_regress_type, args.shift_TRs, ...
            args.level, args.rest)),'regs');
    end
end
end