function[] = clearmem_operation_importance_map_subj(args, dirs)

% individual-level (1st level)
%   1. Negative and Positive maps are calculated separately
%   2. combine negative absolute values into positive values
%   3. z-score the combined values
%   4. normalize to MNI-space

%% ============= UNPACK ARGS.
xph               = args.xphase;
args.regress_type = args.train_regress_type;

%*************** output basename
basename          = args.analysis_basename;

%*************** set FSL environment
setenv('FSLDIR', dirs.fsl);  % this to tell where FSL folder is
setenv('FSLOUTPUTTYPE', 'NIFTI_GZ'); % this to tell what the output type would be

%% ============= SETUP FILE NAMES
ph1.basename = sprintf('%s_%s_%s', args.phase_name{xph}, args.mask_name, args.epi_name);
ph2.basename = sprintf('%s_%s_tr%s_blk_%s',...
    ph1.basename, args.regress_type, ...
    num2str(args.shift_TRs), args.rest);
ph3.basename = sprintf('%s_featsel_thresh%s', ph2.basename, num2str(args.featSelThresh));
ph4.basename = sprintf('%s_decoding_setup', ph3.basename);

%*************** basename phase 5.
class_basename   = sprintf('classified_%s_%s', ph4.basename, args.classifier);
penalty_basename = sprintf('classified_%s_%s', ph3.basename, args.classifier);

%% ============= load penalty_check
load(sprintf('%s/penalty_check_%s.mat', dirs.mvpa.output{xph}, penalty_basename));%'pen_check'
[xacc, whichmax] = max(pen_check.performance); %#ok<*NODEF>
max_penalty      = pen_check.penalty(whichmax);
args.xpenalty    = max_penalty;

fprintf('... max_penalty: %s: acc: %1.4f\n', num2str(max_penalty), xacc);

%*************** load phase 4
fname        = sprintf('%s/ph4_%s.mat', dirs.mvpa.scratch{xph}, ph4.basename);
load(fname);%'ph4'

%*************** load phase 5:results
ph5.basename = sprintf('%s_penalty%s', class_basename, num2str(args.xpenalty));
fname        = sprintf('%s/%s.mat', dirs.mvpa.output{xph}, ph5.basename);
load(fname);%'ph5'

fprintf('\n... loaded classification results (penalty: %s) of %s: %s\n', ...
    num2str(args.xpenalty), args.subject_id, fname);

%% ============= 1ST-LEVEL IMPORTANCE MAP
xdir = fullfile(dirs.epi_mid, 'warp_param');
if ~isdir(xdir), mkdir(xdir); end

xname = sprintf('subj_imp_map_%s_%s', args.imp_type, basename);
fname = sprintf('%s/%s.mat', dirs.mvpa.imp_map{xph}, xname);

%*************** running 1st level importance map
fprintf('(+) 1st level importance map: %s\n', args.subject_id);

[subj] = create_importance_maps_logreg(ph4.subj, ph5.results, args, dirs); %#ok<*NASGU>

save(fname, 'subj','-v7.3');

end