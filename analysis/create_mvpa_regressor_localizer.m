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
    regress_name = category_name;    
elseif strcmp(args.level, 'subcategory')
    for xcate = 1:n_category
        regress_name = [regress_name subcategory_name{xcate}];
    end %#ok<*AGROW>    
elseif strcmp(args.level, 'item')
    for xcate = 1:n_category
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
    
    %% =============== SAVE REGRESSOR
    xname = sprintf('%s_%dtr_%s_%s.mat', ...
            args.train_regress_type, args.shift_TRs, args.level, args.rest);
    fname = fullfile(dirs.mvpa.regressor{xph}, ...
        sprintf('localizer_classification_regressor_%s', xname));
    
    save(fname,'regs','-v7.3');   
    
end
end