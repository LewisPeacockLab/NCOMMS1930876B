function[] = clearmem_mvpa_decode_05(args, dirs)
% second level analysis
% Matlab R2014a

%% ============= UNPACK ARGS.
xph               = args.xphase;
mask_name         = args.mask_name;
args.regress_type = args.train_regress_type;
param             = args.index{xph}.param;
subject_list      = args.subject_list;
xsub_grp          = args.filtered_subs;

%*************** subject id num
for xsub = 1:args.n_sub
    sub_id(xsub) = str2double(subject_list(xsub).name(end-2:end));
end

%*************** filtered subject id num
for it = 1:length(xsub_grp)
    xsub = args.filtered_subs(it);
    sub_id_filtered(it) = str2double(subject_list(xsub).name(end-2:end));
end

%*************** output basename
basename          = args.analysis_basename;
xoutput_dir       = dirs.mvpa.group.auc{xph};

n_subcategory     = param.n_subcategory;
n_category        = param.n_category;
it_categories     = 1:n_category;

if strcmp(args.level, 'subcategory') && args.class_selecting
    it_categories = args.selected_category;
end

if strcmp(args.level, 'category')
    %*************** create tables
    xcate_name = param.category_name;
    n_class    = n_category;
elseif strcmp(args.level, 'subcategory')
    %*************** create tables
    for xcate = 1:length(it_categories)
        itcate = it_categories(xcate);
        for xsubcate = 1:n_subcategory
            xcate_name{xsubcate + n_subcategory*(xcate-1)} = ...
                param.subcategory_name{itcate}{xsubcate};
        end
    end
    
    n_class = n_category * n_subcategory;
end

if strcmp(args.rest, 'rest'), xcate_name{end + 1} = 'rest'; end

n_target     = length(xcate_name);
xsub_groups  = args.filtered_subs;

xcond_color  = args.cond_color;
xcate_color  = args.cate_color;
xbase_color  = args.base_color;

%*************** output basename
basename          = args.analysis_basename;
        
%% ============= SETUP FILE NAMES
%*************** ph4. base file name
if strcmp(args.regress_type, 'shift')
    ph4.basename = sprintf('%s_sh%d_%s_fselected_%s_%s_%s_%s_zepi', ...
        args.phase_name{xph}, args.shift_TRs, args.rest, mask_name, ...
        args.featSelThresh, args.level, args.epi_name); 
elseif strcmp(args.regress_type, 'beta')
    ph4.basename = sprintf('%s_fselected_%s_%s_%s_%s_zepi', ...
        args.phase_name{xph}, mask_name, ...
        args.featSelThresh, args.level, args.beta_name); 
end

%*************** reset ph4. filename
if args.class_selecting
    ph4.basename = sprintf('cate%s_%s', sprintf('%d',args.selected_category), ph4.basename);
end

%*************** ph5. base filename
ph5.basename = sprintf('%s_%s', ph4.basename, args.regress_type);

%*************** ph6. base filenames
ph6.basename     = sprintf('%s_decoding_setup_train+test', ph5.basename);

%*************** classifier name
class_basename  = sprintf('decoding_%s_%s', ph5.basename, args.classifier);

args.grp_results_name = class_basename;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%============= save MVPA results in a group structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if args.mvpa_out{xph}
    %% ============= LOAD EXISTING FILE
    %*************** load mvpa_results
    fname = sprintf('%s/mvpa_results_%s.mat', dirs.mvpa.group.mvpa_result{xph}, basename);
    if exist(fname, 'file'), load(fname); end % 'mvpa_results'
    
    %% ============= SAVE FILE
    %% *************** collect 1st level results
    for xsub = args.filtered_subs
        %*************** setup subject & directories
        args.subject_id = subject_list(xsub).name;
        dirs            = setup_directory(dirs, args);
        
        %% ============= LOAD LOCALIZER PENALTY
        %*************** load phase 6.
        fprintf(sprintf('... loading penalty from localizer: s_%s_%s\n', num2str(xsub), args.subject_id));
        
        %*************** loc_ph1. base filename
        loc_ph1.basename = sprintf('%s_%s_zscored_%s', args.phase_name{1}, args.mask_name, args.epi_name);
        
        %*************** loc_ph2. base filename
        if strcmp(args.regress_type, 'shift')
            loc_ph2.basename = sprintf('%s_%s_%s%dtr_blk_%s',...
                loc_ph1.basename, args.level, args.regress_type, ...
                args.shift_TRs, args.rest);
        elseif strcmp(args.regress_type, 'beta')
            loc_ph2.basename = sprintf('%s_%s_%s', ...
                loc_ph1.basename, args.level, args.regress_type);
        end
        
        %*************** loc_ph3. base filenames
        if args.featVox
            loc_ph3.basename = sprintf('%s_featsel_%svox', loc_ph2.basename, num2str(args.fsVosNum));
        else
            loc_ph3.basename = sprintf('%s_featsel_thresh%s', loc_ph2.basename, num2str(args.featSelThresh));
        end
        
        %*************** loc_ph4. base filenames
        loc_class_basename = sprintf('classified_%s_%s', loc_ph3.basename, args.classifier);
        
        %*************** load penalty check of localizer
        loc_penalty      = load(sprintf('%s/penalty_check_%s.mat', dirs.mvpa.output{1}, loc_class_basename));%'penalty_check'
        [xacc, whichmax] = max(loc_penalty.pen_check.performance); %#ok<*NODEF>
        max_penalty      = loc_penalty.pen_check.penalty(whichmax);
        args.penalty     = max_penalty;
        
        fprintf('... max_penalty: %s: acc: %1.4f\n', num2str(max_penalty), xacc);
        
        %*************** reset xpenalty
        args.xpenalty = args.penalty;
        
        %% ============= MVPAOUT
        %*************** ph7. basename
        ph7.basename = sprintf('%s_penalty%s', class_basename, num2str(args.xpenalty));
        mvpa_result_name = sprintf('%s/%s.mat', dirs.mvpa.output{xph}, ph7.basename);
        
        %*************** result
        fprintf('\n... loading classification results (penalty: %s) of %s: %s\n', ...
            num2str(args.xpenalty), args.subject_id, mvpa_result_name);
        
        xresult = load(mvpa_result_name);%'ph7'
        mvpa_results{xsub} = xresult.ph7.results; %#ok<*AGROW,*NASGU>
        
    end
    
    %%*************** 2ND LEVEL SAVE
    fprintf('(+) save 2nd level mvpa results\n');
    
    fname = sprintf('%s/mvpa_results_%s.mat', dirs.mvpa.group.mvpa_result{xph}, basename);
    save(fname, 'mvpa_results','-v7.3');
    
else
    
    %% ============= LOAD GROUP MVPA RESULTS
    %% *************** setup subject & directories
    args.subject_id = subject_list(1).name;
    dirs            = setup_directory(dirs, args);
    
    %%*************** 2ND LEVEL LOAD
    fprintf('(+) load 2nd level mvpa results\n');
    
    fname = sprintf('%s/mvpa_results_%s.mat', dirs.mvpa.group.mvpa_result{xph}, basename);
    load(fname);%'mvpa_results'
    
    g_fsize = dir(fname);
    fprintf('file size: %s GB\n', num2str(g_fsize.bytes/(10^9)));
    
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%============= ROC for single subject + save mvpa results in group structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


end