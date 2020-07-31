function[] = create_motion_plot(args, dirs)
% plot motion correction

%% ============= UNPACK ARGS.
phase_name = args.phase_name;

%% ============= LOAD MOTION CORRECTION PARAMETERS
for xph = args.xphase%1_loc,2_enco
    
    n_run{xph} = length(dirs.runs{xph});
    t_param    = [];
    
    for xrun=1:n_run{xph}
        clear t_par
        t_par   = load(fullfile(dirs.runs{xph}{xrun}, 'bold_mcf.par'));
        t_param = vertcat(t_param, t_par); %#ok<*AGROW>
        
        n_vol{xph}(xrun) = size(t_par,1);%per run
    end
    
    run_param  = zeros(sum(n_vol{xph}), n_run{xph}-1);
    
    if xph==1, n_par_run = n_run{xph}; else n_par_run = n_run{xph}-1; end
    
    for xrun=1:n_par_run
        run_param((1:n_vol{xph}(xrun)) + sum(n_vol{xph}(1:xrun-1)), xrun) = 1;
    end
 
    param{xph} = horzcat(t_param, run_param);
    
    %%========== SAVE MOTION REGRESSOR
    fname = fullfile(dirs.motion, ['motion_' phase_name{xph} '.txt']);
    dlmwrite(fname, param{xph}, 'delimiter', ' ')
   
end

%% ============= PLOT
cat_param = [param{1}(:,1:6); param{2}(:,1:6)];

fig = figure;
set(fig, 'Position', [0 0 1500 1000])

subplot(2,1,1)
h1 = plot(cat_param(:,1), 'c'); hold on
h2 = plot(cat_param(:,2), 'm'); hold on
h3 = plot(cat_param(:,3), 'b'); hold on
legend([h1, h2, h3], 'x', 'y', 'z');
plot([0 size(cat_param,1)], [0 0], '--k')

%*************** run boundary
for xph = args.xphase
    for xrun = 1:n_run{xph}
        xx = n_vol{xph}(xrun) + sum(n_vol{xph}(1:xrun-1)) + (sum(n_vol{1}) * (xph -1)); 
        plot([xx xx], [min(min(cat_param(:,1:3))) max(max(cat_param(:,1:3)))],'--k');
    end
end

title('Rotation')

subplot(2,1,2)
h4 = plot(cat_param(:,4), 'c'); hold on
h5 = plot(cat_param(:,5), 'm'); hold on
h6 = plot(cat_param(:,6), 'b'); hold on
legend([h4, h5, h6], 'x', 'y', 'z');
plot([0 size(cat_param,1)], [0 0], '--k')

%*************** run boundary
for xph = args.xphase
    for xrun = 1:n_run{xph}
        xx = n_vol{xph}(xrun) + sum(n_vol{xph}(1:xrun-1)) + (sum(n_vol{1}) * (xph -1)); 
        plot([xx xx], [min(min(cat_param(:,4:6))) max(max(cat_param(:,4:6)))],'--k');
    end
end

title('Translation')

%*************** save plot
fig_fname = fullfile(dirs.motion,'motion_correction');

savefig(fig, sprintf('%s.fig',fig_fname));
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 14 7])
print('-djpeg', sprintf('%s.jpg',fig_fname), '-r100')
saveas(fig, sprintf('%s.jpg',fig_fname), 'jpg')

close(fig);

end