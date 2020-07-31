function[] = clearmem_impmap_rsa_decode_01(args, dirs)
%*************** create template by averaging the extracted pattern
%*************** category level

%% ============= SETUP DIRECTORY
output_dir      = dirs.rsa.pattern;

%% ============= SETUP PARAMETERS
n_category      = args.index{1}.param.n_category;
n_subcategory   = args.index{1}.param.n_subcategory;
category_names  = {'face','fruit','scene'};
conds_names     = args.index{2}.param.conds_names;
n_condition     = length(conds_names);
n_trials        = args.index{2}.param.n_trials;%of study
cate_members    = 1:n_category;
n_trs           = args.index{2}.param.n_tc_trs;

if strcmp(args.level, 'category')
    n_class = n_category;
elseif strcmp(args.level, 'subcategory')
    n_class = n_category * n_subcategory;
end
    
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%============= LOCALIZER: 1ST LEVEL 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ============= UNPACK ARGS.
xph     = 1;
xmatrix = args.regs{xph}.matrix;%3180
xheader = args.regs{xph}.header;
n_runs  = args.index{xph}.param.n_runs;

%% ============= LOADING MVPA RESULTS
%*************** ph1. base filename
ph1.basename = sprintf('%s_%s_zscored_%s', args.phase_name{xph}, args.mask_name, args.epi_name); 

%*************** ph2. base filename
ph2.basename = sprintf('%s_%s_%s%dtr_blk_%s',...
    ph1.basename, args.level, args.train_regress_type, ...
    args.shift_TRs, args.rest);

%*************** ph3. base filenames
ph3.basename   = sprintf('%s_featsel_thresh%s', ph2.basename, num2str(args.featSelThresh));

%*************** ph4. base filenames
class_basename = sprintf('classified_%s_%s', ph3.basename, args.classifier);

%*************** load penalty
load(sprintf('%s/penalty_check_%s.mat', dirs.mvpa.output{xph}, class_basename));%'penalty_check'
max_penalty   = penalty_check.penalty(penalty_check.performance==max(penalty_check.performance)); %#ok<*NODEF>
args.xpenalty = min(max_penalty);

%*************** basename phase 4.
ph4.basename  = sprintf('%s_penalty%s', class_basename, num2str(args.xpenalty));

%*************** load phase 4.
fname = sprintf('%s/%s.mat', dirs.mvpa.output{xph}, ph4.basename);

fprintf('\n... loading classification results (penalty: %s) of %s: %s\n', ...
    num2str(args.xpenalty), args.subject_id, fname);

load(fname);%'ph4'

subj    = ph4.subj; 
results = ph4.results;

%% ============= CREATE MASK
%*************** find overlapped features selected across runs (iteration)
n_iters = length(results.iterations);

%*************** Get pats and weights for each crossval iteration

for xclass = 1:n_class
    %*************** create imp_map vol
    t_imp_map_conds{xclass} = zeros([size(xmask_iter) n_class]);
    
    for i = 1:n_iters
        %*************** iteration mask
        xmask_iter = subj.masks{i + 1}.mat;
        
        %*************** get weight
        iter    = results.iterations(i);
        weights = iter.scratchpad.weights;
        
        [n_voxels_iter, n_conds_iter] = size(weights);
        
        
        
        
        xmask_vol(xmask_vol==1) = 
    end
end
%% 
    
    
for i = 1:n_iters
    
    
    
    
    % Get data for this iteration
    
    maskname  = iter.created.maskname;
    
    % Grab network weights
    
    weights = iter.scratchpad.weights;
    
    [n_voxels_iter, n_conds_iter] = size(weights);
    
    assert(n_voxels_pats == n_voxels_iter,'[*] n_voxels dont match');
    assert(n_conds_regs == n_conds_iter,'[*] n_conds dont match');
    
    for j = 1:n_conds_regs
        
        cond_weights   = weights(:,j);
        
        pos_pats = cond_pats_mean > 0;
        neg_pats = cond_pats_mean < 0;
        
        pos_weights = cond_weights > 0;
        neg_weights = cond_weights < 0;
        
        pos_pos = pos_pats & pos_weights;
        neg_neg = neg_pats & neg_weights;
        
        cond_impmap = 0*cond_weights;
        cond_impmap(pos_pos) = 100 * cond_pats_mean(pos_pos) .* cond_weights(pos_pos);
        cond_impmap(neg_neg) = -100 * cond_pats_mean(neg_neg) .* cond_weights(neg_neg);
        cond_impmap_normed   = cond_impmap/n_iters;
        
        % Save a separate impmap for each condition & each iteration (diff #'s of voxels!)
         
        %*************** reset group name
        if args.class_selecting
            impmap_group = sprintf('impmap_%s_%s_%s_cate%s_cond%d', ...
                args.level, args.rest, args.imp_type, sprintf('%d',args.selected_category), j);
        else
            impmap_group = sprintf('impmap_%s_%s_%s_cond%d', args.level, args.rest, args.imp_type, j);
        end
        
        impmap_name  = sprintf('%s_%d',impmap_group,i);
        
        subj = initset_object(subj,'pattern', impmap_name, cond_impmap_normed, ...
            'masked_by', maskname, 'group_name', impmap_group);
        
        created.by = mfilename;
        subj = add_created(subj,'pattern',impmap_name,created);
        
        %*************** reset header info for converting to nifty file
        xheader = get_objfield(subj,'pattern', iter.created.patname,'header');
        subj    = set_objfield(subj,'pattern', impmap_name, 'header', xheader);
        
        %*************** create empty 3d bold with reference brain
        [~, refer_head, refer_vol] = icatb_read_gzip_nii(sprintf('%s.gz', refer_mask));

        % capture mask & create new volume
        mask_matrix = get_mat(subj,'mask', maskname);
        mask_dims   = size(mask_matrix);
        
        masked_pattern_matrix = zeros(mask_dims(1),mask_dims(2),mask_dims(3),1);
        relative_map_of_mask  = find(mask_matrix);
        masked_pattern_matrix(relative_map_of_mask) = cond_impmap_normed;  %#ok<*FNDSB>
        
        ximp_name = sprintf('%s.nii', impmap_name);
        refer_vol.fname = fullfile(xout_dir, ximp_name);
        
        spm_write_vol(refer_vol, masked_pattern_matrix);
        
        %*************** set mat_vol with importance map
        xvol.pat_mat = masked_pattern_matrix;
        xvol.V       = refer_vol;
        xvol.header  = refer_head;
        
        subj = set_objfield(subj,'pattern', impmap_name, 'mat_vol', xvol, 'ignore_absence', true);
        
        delete(fullfile(xout_dir, ximp_name));
    end
end


%% ============= extract patterns & weight (classification)
% for xcate = 1:n_category
%     %*************** z-scored patterns
%     xpats{xph}.allTRs = subj.patterns{2};
%     
% end





end




