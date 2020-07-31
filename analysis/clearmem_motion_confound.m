function[] = clearmem_motion_confound(args, dirs)
%% ============= UNPACK ARGS.
subject_list      = args.subject_list;
xsub_groups       = args.filtered_subs;

for it = 1:length(xsub_groups)
    xsub = xsub_groups(it);
    %*************** setup subject & directories
    args.subject_id = subject_list(xsub).name;
    dirs            = setup_directory(dirs, args);
    
    %% =============== SPIKE ARRAY
    for xph = 1:2
        clear spike
        xspike_file = fullfile(dirs.param, sprintf('spike_%s_%s.mat', ...
            args.phase_name{xph}, args.subject_id));
        
        load(xspike_file);%, 'spike');
        n_spike(it, xph) = size(find(spike),2);
    end
end

save(fullfile(dirs.protocols, 'n_spikes.mat'), 'n_spike');
fprintf('%s was saved\n', fullfile(dirs.protocols, 'n_spikes.mat'));

end