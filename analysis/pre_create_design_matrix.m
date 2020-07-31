%***************** creating TR-based design matrix 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOCALIZER DESIGN MATRIX
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;

% addpath /Applications/spm12/
% addpath ~/github/lewpealab/mvpa/
% addpath(genpath('~/github/GroupICATv4.0b/icatb/'))
    
addpath ./homemade/
addpath ./spm_custom/

%% =============== SETUP DIRECTORY
curr_dir       = pwd; cd ~
dirs.home      = pwd; cd(curr_dir);
dirs.project   = fullfile(dirs.home,'lewpealab_dpbx','STUDY','Clearmem');
dirs.fmri_data = fullfile(dirs.project,'imaging_data','utaustin');%'utaustin'

dirs.protocols = fullfile(dirs.fmri_data, 'protocols');

%% =============== READ PROTOCOLS / PARSE THE PARAMETERS
%***************** run: #/5, block: 18/run, mini_trial: 3/block, trial: 54/run
%***************** category x subcategory x item: 3 x 3 x 6
%***************** 460ms tr, 805 volumes, 
%***************** stim_dur: 3tr, iti: jitter 5:10

fprintf('(+) creating design matrix for localizer\n')

xprotocol_xls            = fullfile(dirs.protocols, 'funcloc_design.xlsx');
[xnum_matrix, read_text] = xlsread(xprotocol_xls);
xtext_matrix             = read_text(2:end, :);
read_header              = read_text(1,:);

%% =============== EXTRACT PARAMETERS
%***************** event
event_name    = unique(xtext_matrix(:, ismember(read_header, 'category')))';
category_name = event_name(~ismember(event_name, {'fix_ITI','fix_block','instruct'}));%face,fruit,land

for xcate = 1:length(category_name)
    %************* category unit
    xunit_cate{xcate}       = ismember(xtext_matrix(:, ismember(read_header, 'category')), category_name{xcate}); %#ok<*SAGROW,*AGROW>

    %************* subcategory names
    subcategory_name{xcate} = unique(xtext_matrix(xunit_cate{xcate}, ismember(read_header, 'subcategory')))';
    
    for xsubcate = 1:length(subcategory_name{xcate})
        %********* subcategory unit
        xunit_subcate{xcate}{xsubcate} = ...
            ismember(xtext_matrix(:, ismember(read_header, 'subcategory')), subcategory_name{xcate}{xsubcate});
        %********* item names
        item_name{xcate}{xsubcate} = ...
            unique(xtext_matrix(xunit_subcate{xcate}{xsubcate}, ismember(read_header, 'image')))';
        
        for xitem = 1:length(item_name{xcate}{xsubcate})
            %***** subcategory unit
            xunit_item{xcate}{xsubcate}{xitem} = ...
                ismember(xtext_matrix(:, ismember(read_header, 'image')), item_name{xcate}{xsubcate}{xitem});
        end    
    end
end

%***************** nmber of parameters
n_category    = length(category_name);
n_subcategory = length(subcategory_name{1});
n_item        = length(item_name{1}{1});

n_runs        = max(xnum_matrix(:,ismember(read_header,'run')));
n_blocks      = max(xnum_matrix(:,ismember(read_header,'block')));
n_trials      = max(xnum_matrix(:,ismember(read_header,'trial')))/2;
n_minitrials  = 3;
n_volumes     = 805;

%***************** durations
tr_in_ms      = 460; % 0.46s 
dur_stim      = 3;
dur_ibi       = 13;
dur_run       = n_volumes * (tr_in_ms/1000);% in sec
dur_jitter    = unique(xnum_matrix(ismember(xtext_matrix(:,ismember(read_header, 'category')),'fix_ITI'), ...
    ismember(read_header,'TRs')))';

%% *************** trial-based matrix
t_header = {'run','block','trial','category','subcategory','item','image_id','onsets'};
t_matrix = zeros(size(xnum_matrix, 1), length(t_header));

xcols = {'run','block','onsets'};
t_matrix(:, findCol(t_header, xcols)) = xnum_matrix(:, findCol(read_header, xcols)); 

for xcate = 1:n_category
    t_matrix(xunit_cate{xcate}, ismember(t_header, 'category')) = xcate;
    
    for xsubcate = 1:n_subcategory
        t_matrix(xunit_subcate{xcate}{xsubcate}, ismember(t_header, 'subcategory')) = xsubcate;
        
        for xitem = 1:n_item
            xitem_id = xitem + (n_item * (xsubcate-1)) + ...
                ((n_item * n_subcategory) * (xcate-1));
            
            t_matrix(xunit_item{xcate}{xsubcate}{xitem}, ismember(t_header, 'item'))     = xitem;
            t_matrix(xunit_item{xcate}{xsubcate}{xitem}, ismember(t_header, 'image_id')) = xitem_id;
            
            image_names{xitem_id, 1} = xitem_id;
            image_names{xitem_id, 2} = item_name{xcate}{xsubcate}{xitem}(1:end-4);
        end
    end
end

%***************** trial numbers
for xrun = 1:n_runs
    xtrial = 0;
    for xrow = find(t_matrix(:, ismember(t_header, 'run')) == xrun)'
        
        if sum(t_matrix(xrow, findCol(t_header, {'category','subcategory'})))
            xtrial = xtrial + 1;
            t_matrix(xrow, ismember(t_header, 'trial')) = xtrial;
        end
    end
end

%***************** save trial index
trial_table = array2table(t_matrix, 'VariableNames', t_header);
writetable(trial_table, fullfile(dirs.protocols, 'localizer_trial_design_index.csv')); 

%***************** save image index
image_table = cell2table(image_names, 'VariableNames', {'image_id','image_name'});
writetable(image_table, fullfile(dirs.protocols, 'localizer_image_index.csv')); 

%% =============== TR BASED DATA MATRIX
%***************** volume in TR / it_volumes: volume id
%***************** block: 0_ibi/instruction 
%***************** trial: 0_ibi, including iti
%***************** stimulus: 1_presentation, 0_rest (fixation)
%***************** onset_time: onset in ms

xheader = {'volume','run','block','trial','category',...
    'subcategory','item','image_id','stimulus'};
xmatrix = zeros(length(xheader), n_volumes * n_runs);

xindex  = 0:(n_volumes * n_runs)-1;

xmatrix(findCol(xheader, {'volume'}), :) = mod(xindex,n_volumes) + 1;
xmatrix(findCol(xheader, {'run'}), :)    = fix(xindex/n_volumes) + 1;

for xrun = 1:n_runs
    for xblock = 1:n_blocks
        clear onset 
        %***************** extracting onset
        xunit  = find(getDATA(t_matrix, t_header, {'run','block'}, {xrun, xblock}));
        
        onset(1) = t_matrix(xunit(1), findCol(t_header, {'onsets'}));
        onset(4) = t_matrix(xunit(end) + 1, findCol(t_header, {'onsets'})) - 1;
        
        for xminitrial = 2:n_minitrials;
            xtrial = xminitrial + (n_minitrials * (xblock-1));
            
            onset(xminitrial) = ...
                getDATA(t_matrix, t_header, {'run','block','trial'}, ...
                {xrun, xblock, xtrial}, findCol(t_header, {'onsets'}));
        end
        
        onset = onset + (n_volumes * (xrun-1));
        
        %***************** block
        xmatrix(findCol(xheader, {'block'}), onset(1):onset(4)) = xblock;
        
        %***************** trial
        for xminitrial = 1:n_minitrials;
            xtrial = xminitrial + (n_minitrials * (xblock-1));
            
            xmatrix(findCol(xheader, {'trial'}), ...
                onset(xminitrial):onset(xminitrial+1)) = xtrial;
        end
    end
end

%***************** volume-based design matrix
for xrun = 1:n_runs
    for xtrial = 1:n_trials
        
        xunit     = getDATA(t_matrix, t_header, {'run','trial'}, {xrun, xtrial});
        xcate     = t_matrix(xunit, findCol(t_header,{'category'}));
        xsubcate  = t_matrix(xunit, findCol(t_header,{'subcategory'}));
        xitem     = t_matrix(xunit, findCol(t_header,{'item'}));
        ximage_id = t_matrix(xunit, findCol(t_header,{'image_id'}));
        
        xunit     = find(getDATA(xmatrix', xheader, {'run','trial'}, {xrun, xtrial}));
        xonset    = xunit(1);
        end_onset = xunit(end);
        
        xmatrix(findCol(xheader,{'category'}), xonset:end_onset)           = xcate;
        xmatrix(findCol(xheader,{'subcategory'}), xonset:end_onset)        = xsubcate;
        xmatrix(findCol(xheader,{'item'}), xonset:end_onset)               = xitem;
        xmatrix(findCol(xheader,{'image_id'}), xonset:end_onset)           = ximage_id;
        xmatrix(findCol(xheader,{'stimulus'}), xonset:(xonset+dur_stim-1)) = 1;
    end
end

%% =============== save regressor table
xtable = array2table(xmatrix', 'VariableNames', xheader);
writetable(xtable, fullfile(dirs.protocols, 'localizer_volume_design_matrix.csv'));

%% =============== SETUP PARAMETER STRUCTURE
clear param

param.category_name     = category_name;%'face', 'fruit', 'land'
param.subcategory_name  = subcategory_name;
param.item_name         = item_name;
param.n_category        = n_category;% 3
param.n_subcategory     = n_subcategory;% 3/category
param.n_item            = n_item;% 6/subcategory
param.n_volumes         = n_volumes;% 805/run
param.n_runs            = n_runs;% n_run out of all 5 runs
param.n_blocks          = n_blocks;% 18/run
param.n_trials          = n_trials;% 54/run
param.n_minitrials      = n_minitrials;%/block

param.tr                = tr_in_ms;% in ms.
param.dur_stim          = dur_stim;% stimulus presentation in tr
param.dur_ibi           = dur_ibi;% inter-block duration
param.dur_run           = dur_run;% in sec: 370.3s = 6.17min
param.dur_jitter        = dur_jitter;

param.image_index       = image_names; %#ok<*STRNU>

%% =============== save params
fname = fullfile(dirs.protocols, 'localizer_params.mat');
save(fname, 'param', '-v7.3'); 

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STUDY DESIGN MATRIX
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%=============== READ PROTOCOLS / PARSE THE PARAMETERS
%***************** run: #/6, trial: 60/run
%***************** category x subcategory x item: 3 x 3 x 6
%***************** category: faces, fruit, landscapes
%***************** 460ms tr, 1165(trimmed)/1175 volumes, 

fprintf('(+) creating design matrix for study\n')

xprotocol_xls            = fullfile(dirs.protocols, 'workingmem_design.xlsx');
[xnum_matrix, read_text] = xlsread(xprotocol_xls);
xtext_matrix             = read_text(2:end, :);
read_header              = read_text(1,:);

%% =============== EXTRACT PARAMETERS
%***************** event
event_name    = unique(xtext_matrix(:, ismember(read_header, 'cat')))';
category_name = event_name(~ismember(event_name, {''}));%'faces,fruit,landscapes

for xcate = 1:length(category_name)
    %************* category unit
    xunit_cate{xcate}       = ismember(xtext_matrix(:, ismember(read_header, 'cat')), category_name{xcate}); %#ok<*AGROW>
    xunit_new_cate{xcate}   = ismember(xtext_matrix(:, ismember(read_header, 'new_item_cat')), category_name{xcate});
    
    %************* subcategory names
    subcategory_name{xcate} = unique(xtext_matrix(xunit_cate{xcate}, ismember(read_header, 'subcat')))';
    
    for xsubcate = 1:length(subcategory_name{xcate})
        %********* subcategory unit
        xunit_subcate{xcate}{xsubcate} = ...
            ismember(xtext_matrix(:, ismember(read_header, 'subcat')), subcategory_name{xcate}{xsubcate});
        xunit_new_subcate{xcate}{xsubcate} = ...
            ismember(xtext_matrix(:, ismember(read_header, 'new_item_subcat')), subcategory_name{xcate}{xsubcate});
        
        %********* item names
        item_name{xcate}{xsubcate} = ...
            unique(xtext_matrix(xunit_subcate{xcate}{xsubcate}, ismember(read_header, 'image')))';
        
        for xitem = 1:length(item_name{xcate}{xsubcate})
            %***** subcategory unit
            xunit_item{xcate}{xsubcate}{xitem} = ...
                ismember(xtext_matrix(:, ismember(read_header, 'image')), item_name{xcate}{xsubcate}{xitem});
        end    
    end
end

%***************** nmber of parameters
n_category     = length(category_name);
n_subcategory  = length(subcategory_name{1});
n_item         = length(item_name{1}{1});

n_runs         = max(xnum_matrix(:,ismember(read_header,'run')));
n_trials       = max(xnum_matrix(:,ismember(read_header,'trial')));
n_volumes      = 1175;

%***************** conditions 
% cond_names = unique(text_body(:, ismember(header, 'trial_type')));
t_conds        = xtext_matrix(:, ismember(read_header, 'trial_type'));
cond_names     = {'maintain','repCat','repItem','target','global'};%Fix_0
adj_cond_names = {'maintain','repCat','repSubcate','target','global'};
cond_array     = zeros(length(t_conds), 1);

for xcond = 1:length(cond_names);
    cond_array(ismember(t_conds, cond_names{xcond})) = xcond;
end

%***************** operation
% oper_names = unique(text_body(:, ismember(header, 'operation')));
t_opers       = xtext_matrix(:, ismember(read_header, 'operation'));
oper_names    = {'Maintain','Switch','Suppress','Clear'};%Fix_0
oper_array    = zeros(length(t_opers), 1);

for xoper = 1:length(oper_names);
    oper_array(ismember(t_opers, oper_names{xoper})) = xoper;
end

%***************** durations
tr_in_ms      = 460; % 0.46s 
dur_stim      = 6; % in TR
dur_manip     = 6;
dur_iti       = 5:9; % jitter
dur_run       = n_volumes * (tr_in_ms/1000);% in sec

%% *************** trial-based matrix
t_header = {'run','trial','condition','operation',...
    'category','subcategory','item','image_id',...
    'new_category','new_subcategory',...
    'trial_trs','stim_trs','manip_trs','fix_trs','onsets'};
t_matrix = zeros(size(xnum_matrix, 1), length(t_header));

xcols = {'run','onsets'};
t_matrix(:, findCol(t_header, xcols)) = xnum_matrix(:, findCol(read_header, xcols)); 

t_matrix(:, ismember(t_header, 'trial'))     = xnum_matrix(:,ismember(read_header,'series'));
t_matrix(:, ismember(t_header, 'condition')) = cond_array;
t_matrix(:, ismember(t_header, 'operation')) = oper_array;
t_matrix(:, ismember(t_header, 'trial_trs')) = xnum_matrix(:,ismember(read_header,'trial_TRs'));
t_matrix(:, ismember(t_header, 'stim_trs'))  = dur_stim;
t_matrix(:, ismember(t_header, 'manip_trs')) = dur_manip;
t_matrix(:, ismember(t_header, 'fix_trs'))   = xnum_matrix(:,ismember(read_header,'fix_TRs'));
t_matrix(:, ismember(t_header, 'onsets'))    = xnum_matrix(:,ismember(read_header,'onsets'));

xunit = t_matrix(:, ismember(t_header, 'condition')) == 0;
t_matrix(xunit, ismember(t_header, 'stim_trs'))      = 0;

for xcate = 1:n_category
    t_matrix(xunit_cate{xcate}, ismember(t_header, 'category')) = xcate;
    t_matrix(xunit_new_cate{xcate}, ismember(t_header, 'new_category')) = xcate;
    
    for xsubcate = 1:n_subcategory
        t_matrix(xunit_subcate{xcate}{xsubcate}, ismember(t_header, 'subcategory')) = xsubcate;
        t_matrix(xunit_new_subcate{xcate}{xsubcate}, ismember(t_header, 'new_subcategory')) = xsubcate;
        
        for xitem = 1:n_item
            xitem_id = xitem + (n_item * (xsubcate-1)) + ...
                ((n_item * n_subcategory) * (xcate-1));
            
            t_matrix(xunit_item{xcate}{xsubcate}{xitem}, ismember(t_header, 'item'))     = xitem;
            t_matrix(xunit_item{xcate}{xsubcate}{xitem}, ismember(t_header, 'image_id')) = xitem_id;
            
            image_names{xitem_id, 1} = xitem_id;
            image_names{xitem_id, 2} = item_name{xcate}{xsubcate}{xitem}(1:end-4);
        end
    end
end

%***************** save trial index
trial_table = array2table(t_matrix, 'VariableNames', t_header);
writetable(trial_table, fullfile(dirs.protocols, 'study_trial_design_index.csv'));

%***************** save image index
image_table = cell2table(image_names, 'VariableNames', {'image_id','image_name'});
writetable(image_table, fullfile(dirs.protocols, 'study_image_index.csv'));

%% =============== TR BASED DATA MATRIX
%***************** volume in TR 
%***************** trial: 0_FIX, each trial: stim (12 tr) + iti (5:9 tr)
%***************** presentation: 1_stim, 2_manipulation, 0_fixation
%***************** manipulation: 0_stim, 1_instruction + fixation

xheader  = {'volume','run','trial','condition','operation',...
    'category','subcategory','item','image_id',...
    'new_category','new_subcategory','presentation','manipulation'};
xmatrix  = zeros(length(xheader), n_volumes * n_runs);

xindex  = 0:(n_volumes * n_runs)-1;

xmatrix(findCol(xheader, {'volume'}), :) = mod(xindex,n_volumes) + 1;
xmatrix(findCol(xheader, {'run'}), :)    = fix(xindex/n_volumes) + 1;
xmatrix(findCol(xheader, {'manipulation'}), :) = 1;

%***************** volume-based design matrix
for xseq = 1:size(t_matrix,1)
    if t_matrix(xseq, findCol(t_header, {'trial'}))~=0
        xrun       = t_matrix(xseq, findCol(t_header, {'run'}));
        xonset     = t_matrix(xseq, findCol(t_header, {'onsets'})) + (n_volumes * (xrun-1));
        trial_vols = xonset + (0:t_matrix(xseq, findCol(t_header, {'trial_trs'})) - 1);
        stim_vols  = trial_vols(1:dur_stim);
        manip_vols = (1:dur_manip)+ max(stim_vols);
        fix_vols   = trial_vols(~ismember(trial_vols, [stim_vols manip_vols]));
        
        xmatrix(findCol(xheader,{'trial'}), trial_vols)        = t_matrix(xseq, findCol(t_header, {'trial'}));
        xmatrix(findCol(xheader,{'condition'}), trial_vols)    = t_matrix(xseq, findCol(t_header, {'condition'}));
        xmatrix(findCol(xheader,{'operation'}), trial_vols)    = t_matrix(xseq, findCol(t_header, {'operation'}));
        
        xmatrix(findCol(xheader,{'category'}), stim_vols)      = t_matrix(xseq, findCol(t_header, {'category'}));
        xmatrix(findCol(xheader,{'subcategory'}), stim_vols)   = t_matrix(xseq, findCol(t_header, {'subcategory'}));
        xmatrix(findCol(xheader,{'item'}), stim_vols)          = t_matrix(xseq, findCol(t_header, {'item'}));
        xmatrix(findCol(xheader,{'image_id'}), stim_vols)      = t_matrix(xseq, findCol(t_header, {'image_id'}));
        xmatrix(findCol(xheader,{'presentation'}), stim_vols)  = 1;% 1_stim, 2_manipulation, 0_fixation
        xmatrix(findCol(xheader,{'presentation'}), manip_vols) = 2;% 1_stim, 2_manipulation, 0_fixation
        xmatrix(findCol(xheader,{'manipulation'}), stim_vols)  = 0;% 0_stim, 1_instruction + fixation
        
        xmatrix(findCol(xheader,{'new_category'}), [manip_vols fix_vols])    = t_matrix(xseq, findCol(t_header, {'new_category'}));
        xmatrix(findCol(xheader,{'new_subcategory'}), [manip_vols fix_vols]) = t_matrix(xseq, findCol(t_header, {'new_subcategory'}));
    end
end

%% =============== save regressor table
xtable = array2table(xmatrix', 'VariableNames', xheader);
writetable(xtable, fullfile(dirs.protocols, 'study_volume_design_matrix.csv'));

%% =============== SETUP PARAMETER STRUCTURE
clear param

param.category_name     = category_name;%'face', 'fruit', 'land'
param.subcategory_name  = subcategory_name;
param.item_name         = item_name;
param.conds_names       = adj_cond_names;
param.operation_names   = oper_names;

param.n_category        = n_category;% 3
param.n_subcategory     = n_subcategory;% 3/category
param.n_item            = n_item;% 6/subcategory
param.n_volumes         = n_volumes;% 1175/run
param.n_runs            = n_runs;
param.n_trials          = n_trials;% 60/run

param.tr                = tr_in_ms;% in ms.
param.dur_stim          = dur_stim;% stimulus presentation in tr
param.dur_manipulation  = dur_manip;% manipulation instruction presentation in tr
param.dur_iti           = dur_iti;% inter-trial duration
param.dur_run           = dur_run;% in sec: 370.3s = 6.17min
param.hrf_trs           = round(10000/param.tr);%10s
param.n_tc_trs          = dur_stim + dur_manip + param.hrf_trs;%timecourse trs

param.image_index       = image_names;

%% =============== save params
fname = fullfile(dirs.protocols, 'study_params.mat');
save(fname, 'param', '-v7.3'); 
