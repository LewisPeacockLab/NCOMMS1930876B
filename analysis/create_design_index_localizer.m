function[index] = create_design_index_localizer(args, dirs, reset)

if reset
    xph = 1;%localizer
    
    fprintf('(+) creating design index for localizer: %s\n', args.subject_id)
    
    %% =============== SELECTED RUNS
    n_runs        = length(dirs.runs{xph});
    selected_runs = zeros(1, n_runs);
    
    for xrun = 1:n_runs
        selected_runs(xrun) = str2double(dirs.runs{xph}{xrun}(end-1:end));
    end
    
    fprintf('...localizer runs: selected %d runs out of 5 runs: %s\n', n_runs, num2str(selected_runs))
    
    %% =============== SELECTED VOLUMES
    %***************** reset the matrix if # of volume/run is short
    n_fix = 0;%13;
    xcat_volume = []; xvol_idx = [];
    
    for xrun = 1:n_runs
        xepi = fullfile(dirs.runs{xph}{xrun}, sprintf('%s.%s.gz', args.epi_name, args.epiext));
        xvol = icatb_read_gzip_nii(xepi, 'read_hdr_only', 1);
        
        n_volumes(xrun) = xvol.dim(5); %#ok<*AGROW>
        
        if xrun == 1, n_vol_run = n_volumes(1) - n_fix; end
        
        tvolume = (1:n_volumes(xrun)) + sum(n_volumes(1:xrun-1));
        xcat_volume = horzcat(xcat_volume, tvolume);
        
        %***************** first 40TR trimmed volume index
        tvol_idx = zeros(1,n_volumes(xrun));
        tvol_idx((n_fix+1):n_volumes(xrun)) = ...
            (1:(n_volumes(xrun)-n_fix)) + (n_vol_run * (selected_runs(xrun) - 1));
        xvol_idx = horzcat(xvol_idx, tvol_idx);
    end
    
    fprintf('...localizer volumes: selected %s volumes out of 805 volumes: total %s\n', ...
        num2str(n_volumes), num2str(length(xcat_volume)))
    
    %% =============== READ PROTOCOLS / PARSE THE PARAMETERS
    %***************** run: #/5, block: 18/run, mini_trial: 3/block, trial: 54/run
    %***************** category x subcategory x item: 3 x 3 x 6
    %***************** 460ms tr, 805 volumes,
    %***************** stim_dur: 3tr, iti: jitter 5:10
    
    xprotocol_csv = fullfile(dirs.protocols, 'localizer_volume_design_matrix.csv');
    xtable        = readtable(xprotocol_csv);
    
    t_header      = xtable.Properties.VariableNames;
    t_matrix      = table2array(xtable);
    
    %***************** param
    load(fullfile(dirs.protocols, 'localizer_params.mat'));% param
    
    %% =============== DATA MATRIX
    %***************** volume in TR / it_volumes: volume id
    %***************** block: 0_ibi/instruction
    %***************** trial: 0_ibi, including iti
    %***************** stimulus: 1_presentation, 0_rest (fixation)
    %***************** onset_time: onset in ms
    
    xheader = {'volume','it_volume','vol_index','run','it_run','block','trial',...
        'category','subcategory','item','image_id','stimulus','spike'};
    xmatrix = zeros(length(xheader), sum(n_volumes));
    
    for xrun = 1:n_runs
        tunit    = find(getDATA(t_matrix, t_header, {'run'}, {selected_runs(xrun)}));
        src_unit = tunit(1:n_volumes(xrun));
        dst_unit = (1:n_volumes(xrun)) + sum(n_volumes(1:xrun-1));
        
        xmatrix(findCol(xheader, t_header), dst_unit) = ...
            t_matrix(src_unit, findCol(t_header, xheader))';
        
        %***************** run
        xmatrix(findCol(xheader, {'run'}), dst_unit) = xrun;
        
        %***************** it run
        xmatrix(findCol(xheader, {'it_run'}), dst_unit) = ...
            t_matrix(src_unit, findCol(t_header, {'run'}))';
    end
    
    xmatrix(findCol(xheader, {'it_volume'}), :) = 1:sum(n_volumes);
    xmatrix(findCol(xheader, {'vol_index'}), :) = xvol_idx;
    
    %***************** extract new parameters
    for xrun = 1:n_runs
        n_trials(xrun) = max(getDATA(xmatrix', xheader, ...
            {'run'}, {xrun}, findCol(xheader, {'trial'})));
        
        n_blocks(xrun) = max(getDATA(xmatrix', xheader, ...
            {'run'}, {xrun}, findCol(xheader, {'block'})));
    end
    
    %% =============== SETUP PARAMETER STRUCTURE
    param.n_volumes         = n_volumes;
    param.n_runs            = n_runs;% n_run out of all 5 runs
    param.n_trials          = n_trials;
    param.n_blocks          = n_blocks;
    
    param.n_trimmed_volumes = n_volumes - args.xtrim;% 795/run
    param.n_cat_trim_vols   = length(xcat_volume) - (args.xtrim * n_runs);
    
    %% =============== SPIKE ARRAY
    spike = zeros(1, sum(n_volumes));
    for xrun = 1:n_runs
        xspike = load(fullfile(dirs.runs{xph}{xrun},'motion_assess','confound.txt'));
        if ~isempty(xspike)
            xunit        = find(sum(xspike, 2)==1) + sum(n_volumes(1:xrun-1));
            spike(xunit) = 1;
        end
    end
    
    xspike_file = fullfile(dirs.param, sprintf('spike_%s_%s.mat', ...
        args.phase_name{xph}, args.subject_id));
    
    save(xspike_file, 'spike');
    
    xmatrix(findCol(xheader, {'spike'}), :) = spike;
    
    %% =============== TRIM MATRIX
    %***************** trim TRs: 10TRs are trimed
    xunit = [];
    for xrun = selected_runs
        tunit = find(getDATA(xmatrix', xheader, {'it_run'}, {xrun}));
        xunit = vertcat(xunit, tunit((args.trim_trs + 1):end));
    end
    
    %***************** reset volumes
    trimed_matrix = xmatrix(:, xunit);
    trimed_matrix(findCol(xheader, {'volume'}), :) = 1:(length(xcat_volume)-(args.xtrim * n_runs));
    
    %% =============== SAVE
    index.header        = xheader;
    index.matrix        = trimed_matrix;
    index.matrix_notrim = xmatrix;
    index.param         = param;
    
    save(fullfile(dirs.param, sprintf('localizer_parameters_%s.mat', args.subject_id)),'index');
    fprintf('...localizer_parameters_%s.mat was saved\n', args.subject_id);
    
else
    xx = load(fullfile(dirs.param, sprintf('localizer_parameters_%s.mat', args.subject_id)));%,'index');
    index = xx.index;
    fprintf('...localizer_parameters_%s.mat was loaded\n', args.subject_id);
end

end