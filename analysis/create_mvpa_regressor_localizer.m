function[regs] = create_mvpa_regressor_localizer(args, dirs)

%% =============== LOAD DESIGN INDEX
%***************** shifted decoding: shift TRs during decoding

xph              = 1;
regs.header      = args.index{xph}.header;
regs.matrix      = args.index{xph}.matrix;
param            = args.index{xph}.param;

%***************** unpack parameters
n_category       = param.n_category;
n_subcategory    = param.n_subcategory;
n_runs           = param.n_runs;% selected run
n_cat_volumes    = param.n_cat_trim_vols;%size(args.index{xph}.matrix,2);

category_name    = param.category_name;
subcategory_name = param.subcategory_name;
item_name        = param.item_name;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% GENERATE REGRESSOR FOR MVPA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%================ REGRESSOR/SELECTOR FOR SHIFTED EPI
%***************** selectors: runs
regs.selectors = regs.matrix(findCol(regs.header, {'run'}), :);

%***************** defining regressor names
regress_name = [];
if strcmp(args.level, 'category')
    regress_name = category_name(args.selected_category);
    rest_class   = length(args.selected_category) + 1;
    
elseif strcmp(args.level, 'subcategory')
    if ~(args.subclass_selecting)
        for xcate = args.selected_category%n_category
            regress_name = [regress_name subcategory_name{xcate}];
        end %#ok<*AGROW>
    else
        for i = 1:length(args.selected_category)
            xcate = args.selected_category(i);%n_category
            for j = 1:length(args.selected_subcategory)
                xsubcate = args.selected_subcategory(i,j);
                regress_name = [regress_name {subcategory_name{xcate}{xsubcate}}];
            end
        end %#ok<*AGROW>
    end
    rest_class = length(regress_name) + 1;
    
elseif strcmp(args.level, 'item')
    for xcate = args.selected_category%n_category
        for xsubcate = 1:n_subcategory 
            regress_name = [regress_name item_name{xcate}{xsubcate}]; 
        end
    end
end

%% ***************** regressors
regs.regressor_name   = regress_name;

if ~strcmp(args.level, 'item')
    regs.regressors       = zeros(length(regress_name), n_cat_volumes);
    regs.regressors_index = zeros(1, n_cat_volumes);
    
    for xrun = 1:n_runs
        for it_cate = 1:length(args.selected_category)
            xcate = args.selected_category(it_cate);
            
            for it_subcate = 1:length(args.selected_subcategory)
                xsubcate = args.selected_subcategory(it_subcate);
                xunit = getDATA(regs.matrix', regs.header, ...
                    {'run','category','subcategory'}, ...
                    {xrun, xcate, xsubcate});
                
                if strcmp(args.level, 'category')
                    xregs = it_cate;
                elseif strcmp(args.level, 'subcategory')
                    xregs = it_subcate + (length(args.selected_subcategory) * (it_cate-1));
                end
                
                regs.regressors(xregs, xunit) = 1;
                regs.regressors_index(xunit)  = xregs;
            end
        end
    end
    
    %% =============== REST CATEGORY
    if strcmp(args.rest, 'rest')
        xregs = rest_class;
        regs.regressor_name{xregs} = 'rest';
        regs.regressors(xregs,:)   = zeros(1, n_cat_volumes);
        
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
    xname = sprintf('%s_%dtr_%s_%s.mat', ...
            args.train_regress_type, args.shift_TRs, args.level, args.rest);
        
    if args.class_selecting
        xname = sprintf('cate%s_%s', sprintf('%d',args.selected_category), xname);
    end
    
    if args.subclass_selecting
        xname = sprintf('subcate%s_%s', sprintf('%d',args.selected_subcategory), xname);
    end
    
    fname = fullfile(dirs.mvpa.regressor{xph}, ...
        sprintf('localizer_classification_regressor_%s', xname));
    
    save(fname,'regs','-v7.3');   
    
end
end