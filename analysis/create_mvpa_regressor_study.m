function[regs] = create_mvpa_regressor_study(args, dirs)
%% =============== LOAD DESIGN INDEX
%***************** trim the matrix -> shift regressor
%***************** non-shifted decoding: shift TRs in the regressor

xph          = args.xphase;
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
        regress_name = category_name;
    elseif strcmp(args.level, 'subcategory')
        for xcate = 1:n_category
            regress_name = [regress_name subcategory_name{xcate}];
        end %#ok<*AGROW>
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
            
            if strcmp(args.level, 'category')
                xregs = it_cate;
            elseif strcmp(args.level, 'subcategory')
                xregs = xsubcate + (n_subcategory * (it_cate-1));
            end
            
            regs.regressors(xregs, xunit) = 1;
            regs.regressors_index(xunit)  = xregs;
        end
    end
    
    %% =============== SAVE REGRESSOR
    save(fullfile(dirs.mvpa.regressor{xph}, ...
        sprintf('study_classification_regressor_%s_%dtr_%s_%s.mat', ...
        args.train_regress_type, args.shift_TRs, ...
        args.level, args.rest)),'regs');
end
end