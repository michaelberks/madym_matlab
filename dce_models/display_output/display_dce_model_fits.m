function [] = display_dce_model_fits(dyn_signals, model_signals, varargin)
%DISPLAY_DCE_MODEL_FITS visual display of time-series on a gridded sub-plot
%figure
%   [] = display_dce_model_fits(dyn_signals, model_signals, num_voxels, plot_rows, plot_cols)
%
% Inputs:
%      dyn_signals - N_vox x N_times
%
%      model_signals - N_vox x N_times x N_models
%
%      Optional inputs using name/value pairs:
%
%      voxel_idx ([]) - subset of voxels to display, if empty all voxels
%      displayed
%
%      plot_rows (3) - number of rows in figure grid
%
%      plot_cols (5) - number of columns in figure grid
%
%      model_colors ('rgmykc') - color to plot each model time-series
%
%      max_signal (0) - maximum scale on y-axis, if 0, uses max over all input signals
%
%      sse_thresh (1e7) - model fit SSE threshold at which voxels are highlighted
%
%      fig_dir ([]) - if not empty, path to folder where figures will be
%      saved, otherwise figures not saved
%
%
% Outputs:
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 12-Apr-2017
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 

args = u_packargs(varargin, 0, ...
    'voxel_idx', [],...
    'plot_rows', 3, ...
    'plot_cols', 5, ...
    'model_colors', [],...
    'max_signal', 0,...
    'sse_thresh', 1e7,...
    'fig_dir', []);
clear varargin;

%Useful values
[num_voxels, num_time_points, num_models] = size(model_signals);
plots_per_fig = args.plot_rows * args.plot_cols;

if isempty(args.voxel_idx) 
    args.voxel_idx = 1:num_voxels;
end
args.voxel_idx = args.voxel_idx(:);

%Check if we're savig figures
do_save = ~isempty(args.fig_dir);

%Compute max signal to plot
if ~args.max_signal
    args.max_signal = max(max(dyn_signals(:)), max(model_signals(:)));
end

if isempty(args.model_colors)
    args.model_colors = lines(num_models);
end

%Loop through voxels
for i_vox = 1:length(args.voxel_idx)
    vox_idx = args.voxel_idx(i_vox);
    
    if rem(i_vox, plots_per_fig) == 1
        if do_save && i_vox > 1
            exportfig([args.fig_dir 'fig' zerostr((i_vox-1)/plots_per_fig, 3) '.png']);
        end
        figure; plot_num = 1;
    else
        plot_num = plot_num + 1;
    end
    
    subplot(args.plot_rows, args.plot_cols, plot_num); 
    plot(dyn_signals(vox_idx,:)); hold all;
    
    sse_str = ', SSE = ';
    for i_model = 1:num_models
        model_signals_i = model_signals(vox_idx, :, i_model);
        plot(model_signals_i, [args.model_colors(i_model,:) '--'], 'linewidth', 1.5);

        sse = sum((dyn_signals(i_vox,:) - model_signals_i).^2);
        if args.sse_thresh && sse > args.sse_thresh
            plot(...
                [1 num_time_points-1  num_time_points-1 1 1],...
                [1 1 args.max_signal-1 args.max_signal-1 1], 'g', 'linewidth', 3.0);
        end
        sse_str = strcat(sse_str, num2str(sse,3));
        if i_model < num_models
           sse_str = strcat(sse_str, ', ');
        end
    end
    
    set(gca, 'ylim', [0 args.max_signal], 'xlim', [0 num_time_points]);
    if do_save
        title(['V ' num2str(i_vox)]);
    else
        title(['Vox ' num2str(i_vox)  sse_str]);
    end
end
