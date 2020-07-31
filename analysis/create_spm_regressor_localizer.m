function[regs] = create_spm_regressor_localizer(args, dirs)
%%%%%%%%%%%%%%%%%% creating SPM regressors
%%%%%%%%%%%%%%%%%% Multiple condition: names, onsets, durations
%%%%%%%%%%%%%%%%%% Multiple regressors: motion correction

%% =============== LOAD DESIGN INDEX
% cf.: the matrix is for trimmed TR, but the epi is not trimmed
% So, I need to use not trimmed regressor for GLM analyses

xph           = 1;
xheader       = args.index{xph}.header;
xmatrix       = args.index{xph}.matrix_notrim';
xparam        = args.index{xph}.param;

%% *************** unpack parameters
n_category    = xparam.n_category;
n_subcategory = xparam.n_subcategory;
n_runs        = xparam.n_runs;% selected run
n_blocks      = xparam.n_blocks;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% GENERATE REGRESSOR FOR SPM GLM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% names, onsets, durations

%%%%%%%%%%%%%%%%%% CATEGORY LEVEL
if strcmp(args.level, 'category')
    clear names onsets durations
    %***************** names
    names     = {'face','fruit','scene'};
    
    %***************** onsets & durations
    onsets    = cell(1,n_category);
    durations = cell(1,n_category);
    
    for xrun = 1:n_runs
        for xblk = 1:n_blocks(xrun)
            xunit  = getDATA(xmatrix, xheader, {'run','block'},{xrun, xblk});
            
            xcate  = unique(xmatrix(xunit, findCol(xheader, {'category'})));
            xonset = xmatrix(xunit, findCol(xheader, {'it_volume'}));
            
            onsets{xcate}    = horzcat(onsets{xcate}, xonset(1));
            durations{xcate} = horzcat(durations{xcate}, xonset(end) - xonset(1) + 1);
        end
    end
%%%%%%%%%%%%%%%%%% SUBCATEGORY LEVEL    
elseif strcmp(args.level, 'subcategory')
    %***************** names
    for xcate = 1:n_category
        for xsubcate = 1:n_subcategory
            names{xsubcate + (n_subcategory * (xcate - 1))} = ...
                xparam.subcategory_name{xcate}{xsubcate};
        end
    end
    
    %***************** onsets & durations
    onsets    = cell(1, n_category * n_subcategory);
    durations = cell(1, n_category * n_subcategory);
    
    for xrun = 1:n_runs
        for xblk = 1:n_blocks(xrun)
            xunit    = getDATA(xmatrix, xheader, {'run','block'},{xrun, xblk});
            
            xcate = unique(xmatrix(xunit, findCol(xheader, {'category'})));
            xsubcate = unique(xmatrix(xunit, findCol(xheader, {'subcategory'})));
            xonset   = xmatrix(xunit, findCol(xheader, {'it_volume'}));
            
            it_subcate = xsubcate + (n_subcategory * (xcate - 1));
            
            onsets{it_subcate}    = horzcat(onsets{it_subcate}, xonset(1));
            durations{it_subcate} = horzcat(durations{it_subcate}, xonset(end) - xonset(1) + 1);
        end
    end
%%%%%%%%%%%%%%%%%% ITEM LEVEL
elseif strcmp(args.level, 'item')
    clear names onsets durations
    n_item = [];
    
    for xcate = 1:n_category
        for xsubcate = 1:n_subcategory
            
            xunit = getDATA(xmatrix, xheader, ...
                {'category','subcategory'}, {xcate, xsubcate});
            
            it_items = unique(xmatrix(xunit, findCol(xheader, {'item'})));
            n_item(xsubcate + (n_subcategory * (xcate-1))) = length(it_items);
            
            for xitem = it_items'
                it_runs = unique(getDATA(xmatrix, xheader, ...
                    {'category','subcategory','item'}, ...
                    {xcate, xsubcate, xitem},...
                    findCol(xheader, {'run'})));
                
                xcell = xitem + sum(n_item(1:(xsubcate + n_subcategory * (xcate-1)) - 1));
                onsets{xcell} = []; durations{xcell} = [];
                
                %***************** names
                names{xcell}  = xparam.item_name{xcate}{xsubcate}{xitem}(1:end-4);
                
                %***************** onsets & durations
                for xrun = it_runs'
                    xunit  = getDATA(xmatrix, xheader, ...
                        {'run','category','subcategory','item'},...
                        {xrun, xcate, xsubcate, xitem});
                    
                    xonset = xmatrix(xunit, findCol(xheader, {'it_volume'}));
                    
                    onsets{xcell}    = horzcat(onsets{xcell}, xonset(1));
                    durations{xcell} = horzcat(durations{xcell}, 0);
                end
            end
        end
    end
end

%***************** save in the structure
regs.names     = names;
regs.onsets    = onsets;
regs.durations = durations;

%***************** SAVE REGRESSOR
save(fullfile(dirs.rsa.roi.regressor, ...
    sprintf('localizer_spm_regressor_%s.mat', args.level)),'names','onsets','durations','-v7.3');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% MOTION REGRESSOR FOR SPM GLM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(args.subject_id, 'clearmem_v1_sub032') 
    n_all_runs = 4;
elseif strcmp(args.subject_id, 'clearmem_v1_sub033')
    n_all_runs = 3;
else
    n_all_runs = 5;
end

t_motion      = load(fullfile(dirs.motion, 'motion_localizer.txt'));
selected_runs = unique(xmatrix(:, findCol(xheader, {'it_run'})));
motion_matrix = t_motion(:, 1:6);
n_vols        = args.index{xph}.param.n_volumes;
vol_index = ones(1, n_all_runs) * double(n_vols(1));
vol_index(selected_runs(~ismember(n_vols, n_vols(1)))) = ...
    n_vols(~ismember(n_vols, n_vols(1)));

xmotion = []; xrun_regs = [];

for xit = 1:length(selected_runs)
    xrun = selected_runs(xit);
    
    xunit   = (1:n_vols(xit)) + sum(n_vols(1:xit-1));
    if ~strcmp(args.subject_id, 'clearmem_v1_sub033')
        yunit = (1:vol_index(xrun)) + sum(vol_index(1:xrun-1));
    else
        yunit = (1:vol_index(xit)) + sum(vol_index(1:xit-1));
    end
    xmotion = vertcat(xmotion, motion_matrix(yunit, :)); %#ok<*AGROW>
    
    t_run_regs = zeros(length(find(xunit)), n_runs);
    t_run_regs(:, xit) = 1;
    
    xrun_regs  = vertcat(xrun_regs, t_run_regs);
end

%***************** save in the structure
regs.multi_regress = horzcat(xmotion, xrun_regs);

%***************** SAVE REGRESSOR
dlmwrite(fullfile(dirs.rsa.roi.regressor, 'localizer_multi_regressor.txt'),...
    regs.multi_regress);

end