function[] = clearmem_params_creator(args, dirs)

%% ============= PARAMETER 
%*************** setup subject & directories

fprintf('################################\n');
fprintf('%s\n', args.subject_id);
fprintf('################################\n');

%*************** make directories
if ~isdir(dirs.param), mkdir(dirs.param); end

%*************** design parameters
create_design_index_localizer(args, dirs);
create_design_index_study(args, dirs);
    
end

