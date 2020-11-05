function [] = make_temporal_residuals(subject_dir, subject_idx, visit, varargin)
%MAKE_TEMPORAL_RESIDUALS *Insert a one line summary here*
%   [] = make_temporal_residuals()
%
% Inputs:
%
% Outputs:
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 08-Jul-2020
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 
args = u_packargs(varargin, 0, ...
    'dyn_conc_dir', '',...
    'dynamic_dir', '',...
    'n_vols', 100,...
    'analysis_dirs', [],...
    'roi_mask_path',[],...
    'save_dir', [],...
    'use_exiting', 0,...
    'update_xtr_files', 0,...
    'plot', 0);

%Set up paths from input
save_dir = [subject_dir args.save_dir '\'];

if args.use_existing
    load([save_dir 'temporal_residuals.mat'], 'temporal_residuals');
else
    dyn_conc_dir = [subject_dir args.dyn_conc_dir '\'];

    %Set up ROI masks
    hdr = read_analyze_hdr([dyn_conc_dir 'Ct_sig1.hdr']);
    if ~isempty(args.roi_mask_path)
        roi_mask = load_img_volume([subject_dir...
            args.roi_mask_path]) > 0;
    else
        roi_mask = true(hdr.Dimensions([2 1 3]));
    end

    %Set up containers for residuals
    n_models = length(args.analysis_dirs);
    temporal_residuals.MAE = zeros(args.n_vols, n_models);
    temporal_residuals.pct95 = zeros(args.n_vols, n_models);
    temporal_residuals.MSE = zeros(args.n_vols, n_models);

    %Preload the concentration time-series
    C_t = get_dyn_vals([dyn_conc_dir 'Ct_sig'], args.n_vols, roi_mask);

    %Loop through each model, computing model residuals ACROSS all voxels for
    %each time
    fprintf('Making temporal residuals for %d, visit %d\n', ...
        subject_idx, visit);
    for i_m = 1:n_models

        analysis_dir = [subject_dir args.analysis_dirs{i_m}];

        for i_t = 1:args.n_vols
            fprintf('Processing volume %d\n', i_t);
            C_m = load_img_volume([analysis_dir '/Ct_mod' num2str(i_t) '.hdr']);

            %Compute residuals as difference between model fit and signal
            %concentration
            residuals = C_t(:,i_t) - C_m(roi_mask);
            residuals(~isfinite(residuals)) = [];

            %Take median absolute difference, 95th percentile and mean-squared
            %error
            temporal_residuals.MAE(i_t, i_m) = median(abs(residuals));
            temporal_residuals.pct95(i_t, i_m) = prctile(abs(residuals), 95);
            temporal_residuals.MSE(i_t, i_m) = ...
                mean(sum(residuals(...
                residuals < temporal_residuals.pct95(i_t, i_m)).^2));

        end
    end
    create_folder(save_dir);
    save([save_dir 'temporal_residuals.mat'], 'temporal_residuals');
end

%%
if args.update_xtr_files
    %Record the temporal residual for each volume in its xtr file so it can
    %be used in Madym fitting
    dynamic_root = [subject_dir args.dynamic_dir '\dyn_'];
    for i_t = 1:args.n_vols
        noise_sigma = min(temporal_residuals.MSE(i_t, :));
        xtr_path = [dynamic_root num2str(i_t) '.xtr'];
        write_xtr_file(xtr_path, 1, 'NoiseSigma', noise_sigma);
    end    
end
    
%%
if args.plot
    figure;
    subplot(3,1,1);
    plot(temporal_residuals.MAE);
    title(['Subject ' num2str(subject_idx) ': residuals per volume']);
    xlabel('# volume in DCE time-series');
    ylabel('MAE');

    subplot(3,1,2);
    plot(temporal_residuals.pct95);
    xlabel('Time-series # volume');
    ylabel('95^{th} percentile');

    subplot(3,1,3);
    plot(temporal_residuals.MSE);
    xlabel('Time-series # volume');
    ylabel('MSE');

end
%%


