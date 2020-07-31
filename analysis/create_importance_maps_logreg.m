function [subj] = create_importance_maps_logreg(subj, results, args, dirs)

% Calculate discrimination maps
%
% * subj    = subject structure (MVPA)
% * results = results structure from cross_validation analysis in MVPA toolbox
% * type    = 'polyn'|'mcduff' : polyn(neg*neg=pos), mcduff(neg*neg=neg)
% * xheader = reference headear

xph        = args.xphase;
type       = args.imp_type;
refer_mask = args.refer_mask;
n_ori_pats = length(subj.patterns);
imp_type   = {'pos','neg','comb'};
xout_dir   = dirs.mvpa.imp_map{xph};
n_imp_type = length(imp_type);

%% Setup stuff

assert(or(strcmp(type,'polyn'),strcmp(type,'mcduff')),'[*] bad impmap type');
n_iters = length(results.iterations);

%% Get pats and weights for each crossval iteration

xsubj = subj;

for i = 1:n_iters
    
    iter = results.iterations(i);
    
    % Get data for this iteration
    
    maskname  = iter.created.maskname;
    pats      = get_masked_pattern(xsubj, iter.created.patname, maskname);
    selectors = get_mat(xsubj,'selector', iter.created.selname);
    regs      = get_mat(xsubj,'regressors', iter.created.regsname);
    
    % Grab only selected TRs
    
    pats = pats(:,logical(selectors));
    regs = regs(:,logical(selectors));
    
    % Grab network weights
    
    weights = iter.scratchpad.weights;
    
    [n_conds_regs, n_trs_regs] = size(regs);
    [n_voxels_pats, n_trs_pats] = size(pats);
    [n_voxels_iter, n_conds_iter] = size(weights);
    
    assert(n_voxels_pats == n_voxels_iter,'[*] n_voxels dont match');
    assert(n_conds_regs == n_conds_iter,'[*] n_conds dont match');
    assert(n_trs_regs == n_trs_pats,'[*] n_trs dont match');
    
    for j = 1:n_conds_regs
        
        cond_trs       = logical(regs(j,:));
        
        cond_pats      = pats(:,cond_trs);
        cond_pats_mean = mean(cond_pats,2);
        
        cond_weights   = weights(:,j);
        
        switch type
            
            case 'polyn'
                
                % importance map = 100 * avg_input * weight
                
                cond_impmap = 100 * cond_pats_mean .* cond_weights;
                cond_impmap_normed = cond_impmap/n_iters;
                
            case 'mcduff'
                
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
        end
        
        % Save a separate impmap for each condition & each iteration (diff #'s of voxels!)
         
        %*************** reset group name
        if args.class_selecting
            impmap_group = sprintf('impmap_%s_%s_%s_cate%s_cond%d', ...
                args.level, args.rest, args.imp_type, sprintf('%d',args.selected_category), j);
        else
            impmap_group = sprintf('impmap_%s_%s_%s_cond%d', args.level, args.rest, args.imp_type, j);
        end
        
        if (xph==1) && args.reset_regs{xph}
            impmap_group = sprintf('nw_%s', impmap_group);
        end
        
        impmap_name  = sprintf('%s_%d',impmap_group,i);
        
        xsubj = initset_object(xsubj,'pattern', impmap_name, cond_impmap_normed, ...
            'masked_by', maskname, 'group_name', impmap_group);
        
        created.by = mfilename;
        xsubj = add_created(xsubj,'pattern',impmap_name,created);
        
        %*************** reset header info for converting to nifty file
        xheader = get_objfield(xsubj,'pattern', iter.created.patname, 'header');
        xsubj   = set_objfield(xsubj,'pattern', impmap_name, 'header', xheader);
        
        %*************** create empty 3d bold with reference brain
        [~, refer_head, refer_vol] = icatb_read_gzip_nii(sprintf('%s.gz', refer_mask));

        % capture mask & create new volume
        mask_matrix = get_mat(xsubj,'mask', maskname);
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
        
        xsubj = set_objfield(xsubj,'pattern', impmap_name, 'mat_vol', xvol, 'ignore_absence', true);
        
        delete(fullfile(xout_dir, ximp_name));
    end
end

%% ============= MEAN IMPORTANCE MAP PER CATEGORIES
n_patterns     = length(subj.patterns);
xepi_dim       = get_objfield(subj, 'mask', args.mask_name, 'matsize');
ximp_map_conds = cell(1, n_conds_regs);

%*************** mean of importance map per category
for xreg = 1:n_conds_regs
    clear t_imp_map_conds xvol impmap_name
    
    xunit = (xreg:n_conds_regs:(n_iters * n_conds_regs)) + n_ori_pats;
    t_imp_map_conds = zeros([xepi_dim length(xunit)]);
    
    for i = 1:length(xunit)
        xpat = xunit(i);
        it_imp_map = xsubj.patterns{xpat}.mat_vol.pat_mat;
        t_imp_map_conds(:,:,:,i) = it_imp_map;
    end
        
    %*********** mean importance value separately for 
    % positive / negative / combined
    for xmap = 1:n_imp_type, ximp_map_conds{xreg}{xmap} = zeros(size(it_imp_map)); end
    
    for x = 1:size(t_imp_map_conds, 1)
        fprintf('%s.', num2str(x));
        
        for y = 1:size(t_imp_map_conds, 2)
            for z = 1:size(t_imp_map_conds, 3)
                pos_imp = find(t_imp_map_conds(x,y,z,:) > 0);
                neg_imp = find(t_imp_map_conds(x,y,z,:) < 0);
                com_imp = find(t_imp_map_conds(x,y,z,:) ~= 0);
                
                %*********** positive
                if ~isempty(pos_imp)
                    ximp_map_conds{xreg}{1}(x,y,z) = mean(t_imp_map_conds(x,y,z,pos_imp));
                end
                %*********** negative (absolute value)
                if ~isempty(neg_imp)
                    ximp_map_conds{xreg}{2}(x,y,z) = mean(abs(t_imp_map_conds(x,y,z,neg_imp)));
                end
                %*********** combined (absolute value)
                if ~isempty(com_imp)
                    ximp_map_conds{xreg}{3}(x,y,z) = mean(abs(t_imp_map_conds(x,y,z,com_imp)));
                end
            end
        end
    end
    fprintf('\n');
    
    %*************** zscoring combined maps
    for xmap = 1:n_imp_type
        zimp_map_conds{xreg}{xmap} = zeros(size(ximp_map_conds{xreg}{xmap})); %#ok<*AGROW>
    end
    
    xcomb_map = 3;%comb
    xcord = find(ximp_map_conds{xreg}{xcomb_map});
    xx    = ximp_map_conds{xreg}{xcomb_map}(xcord);
    z_xx  = zscore(xx);
    zimp_map_conds{xreg}{xcomb_map}(xcord) = z_xx;
    
    %*************** positive/negative
    for xmap = 1:2
        xcord = find(ximp_map_conds{xreg}{xmap});
        z_xx  = zimp_map_conds{xreg}{xcomb_map}(xcord);
        
        zimp_map_conds{xreg}{xmap}(xcord) = z_xx;
    end
    
    %*************** set object field for condition mean
    cond_name = args.regs{args.xphase}.regressor_name{xreg};
    
    %*************** reset group name
    if args.class_selecting
        impmap_group = sprintf('impmap_%s_%s_cate%s_%s_%s_%s_mean_cond', ...
            args.phase_name{args.xphase}, args.level, ...
            sprintf('%d',args.selected_category), args.mask_name, args.rest, args.imp_type);
    else
        impmap_group = sprintf('impmap_%s_%s_%s_%s_%s_mean_cond', ...
            args.phase_name{args.xphase}, args.level, args.mask_name, args.rest, args.imp_type);
    end
        
    for i = 1:n_imp_type%1_pos, 2_neg, 3_combined
        mean_impmap_name  = sprintf('%s%d_%s_%s', impmap_group, xreg, cond_name, imp_type{i});
        mean_new_filename = fullfile(xout_dir, sprintf('%s.nii', mean_impmap_name));
        
        xunit = n_patterns + i + (n_imp_type * (xreg-1));
        subj.patterns{xunit}.name       = mean_impmap_name;
        subj.patterns{xunit}.header     = xheader;
        subj.patterns{xunit}.pat_mat    = zimp_map_conds{xreg}{i};
        subj.patterns{xunit}.matsize    = size(ximp_map_conds{xreg}{i});
        subj.patterns{xunit}.group_name = impmap_group;
        subj.patterns{xunit}.masked_by  = maskname;
        subj.patterns{xunit}.fname      = mean_new_filename;
        
        %*************** create nii.gz
        xcur_vol       = subj.patterns{xunit}.header.vol{1}(1);
        xcur_vol.fname = mean_new_filename;
        spm_write_vol(xcur_vol, zimp_map_conds{xreg}{i});
        
        gzip(mean_new_filename, xout_dir);
        delete(mean_new_filename);
    end    
end

%% ============= CONVERT SUBJECT TO STANDARD (MNI) 
% # $1: cluster 'blanca' | 'local'
% # $2: imp_map dir 
% # $3: subject's directory
% # $4: importance map basename
% # $5: mni_mask dir
fprintf('(+) normalize subject importance map\n');

system(sprintf('./imp_map_inversetrans_3mm.sh %s %s %s %s %s', ...
    args.cluster, dirs.mvpa.imp_map{xph}, dirs.subj_home, ...
    impmap_group, fullfile(dirs.fmri_data,'mni-masks')))

%% ============= SAVE NORMALIZED PATTERNS
%*************** mean of importance map per category
clear xpat
for xreg = 1:n_conds_regs
    
    cond_name = args.regs{args.xphase}.regressor_name{xreg};
    
    for i = 1:n_imp_type%1_pos, 2_neg
        xname  = sprintf('norm_%s%d_%s_%s', impmap_group, xreg, cond_name, imp_type{i});
        
        xnifty = fullfile(xout_dir, sprintf('%s.nii.gz', xname));
        
        %*************** load the norm nii
        [xpat{xreg}{i}, head, vol] = icatb_read_gzip_nii(xnifty);
        
        xunit = n_patterns + i + (n_imp_type * (xreg-1));
        subj.patterns{xunit}.norm.name    = xname;
        subj.patterns{xunit}.norm.header  = head;
        subj.patterns{xunit}.norm.vol     = vol;
        subj.patterns{xunit}.norm.pat_mat = xpat{xreg}{i};
        subj.patterns{xunit}.norm.matsize = size(xpat{xreg}{i});
        subj.patterns{xunit}.norm.fname   = xnifty;
    end    
end

end
