function [] = display_time_series(time_signals, varargin)
%DISPLAY_TIME_SERIES displays time-series in gridded sub-plot figures
%   [] = display_time_series(time_signals, varargin)
%
% Inputs:
%      time_signals - N_vox x N_times x N_sources
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
%      signal_colors ('rgmykc') - color to plot each model time-series
%
%      min/max_signal (0) - minimum/maximum scale on y-axis, if 0, uses 
%       min/max over all input signals
%
%      t_error ([]) - N_vox x N_signals of fit errors to display
%
%      error_label ([]) - Label to format display of errors
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
    'signal_colors', 'rgmykc',...
    'min_signal', [],...
    'max_signal', [],...
    't_error', [],...
    'error_label', [],...
    'fig_dir', []);
clear varargin;

%Useful values
[num_voxels, num_time_points, num_signals] = size(time_signals);
plots_per_fig = args.plot_rows * args.plot_cols;

if isempty(args.voxel_idx) 
    args.voxel_idx = 1:num_voxels;
end
args.voxel_idx = args.voxel_idx(:);

%Check if we're savig figures
do_save = ~isempty(args.fig_dir);

%Compute max signal to plot
if isempty(args.max_signal)
    args.max_signal = max(time_signals(:));
end
if isempty(args.min_signal)
    args.min_signal = min(time_signals(:));
end

display_errors = ~isempty(args.t_error);

%Loop through voxels
for i_vox = 1:length(args.voxel_idx)
    vox_idx = args.voxel_idx(i_vox);
    
    if rem(i_vox, plots_per_fig) == 1
        if do_save && i_vox > 1
            %exportfig([args.fig_dir 'fig' zerostr((i_vox-1)/plots_per_fig, 3) '.png']);
            saveas(gcf, ...
                [args.fig_dir 'fig' zerostr((i_vox-1)/plots_per_fig, 3) '.png']);
        end
        figure; plot_num = 1;
    else
        plot_num = plot_num + 1;
    end
    
    subplot(args.plot_rows, args.plot_cols, plot_num); 
    hold all;
    
    if display_errors
        err_str = [', ' args.error_label ' = '];
    else
        err_str = [];
    end
    
    for i_signal = 1:num_signals
        plot(time_signals(vox_idx, :,i_signal), ...
            args.signal_colors(i_signal));
        
        if display_errors
            err_str = strcat(err_str, num2str(args.t_error(i_vox, i_signal),3));
            if i_signal < num_signals
               err_str = strcat(err_str, ', ');
            end
        end
    end
    
    set(gca, 'ylim', [args.min_signal args.max_signal], 'xlim', [0 num_time_points]);
    title(['Vox ' num2str(i_vox) err_str]);
end
