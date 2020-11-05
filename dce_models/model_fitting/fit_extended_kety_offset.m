function [offset] = fit_extended_kety_offset(dyn_signals, dyn_times, aif, T1, varargin)
%FIT_EXTENDED_KETY_OFFSET fit just the offset parameter (to determine peak position)
%of the extended-Kety (Tofts) model
%   [offset] = fit_extended_kety_offset(dyn_signals, dyn_times, aif, T1, varargin)
%
% Inputs:
%      dyn_signals - N_vox x N_times array DCE time-series data
%
%      dyn_times - real world times associated with each time-point
%
%      aif - arterial input function
%
%      T1 - baseline T1 for each voxel
%
%      Optional inputs supplied as name/value pairs:
%
%     'FA', 20, FA of dynamic scans
%     'TR', 4, TR of dynamic scans
%     'relax_coeff', 14.0e-3, tissue relaxivity constant
%     'Ktrans_init', 0.397579, Initial values for parameters, either scalar or 1 per vox
%     'Vp_init', 0.015379,...
%     'Ve_init', 0.166493,...
%     'offset_init', 0,...
%     'num_baseline', 5, Number of times to compute baseline mean signal
%     'num_itr', 20, Max number of iterations
%     'opt_display', 'final', Display option during optimisation
%     'plot_rows', 3, Number of rows/cols in figure grid for displaying model fits
%     'plot_cols', 5,...
%     'fig_dir', [], If not empty, saves visual output to specfied folder
%     'debug', false); Switch on debugging visual output
%
% Outputs:
%
%      F_p, v_ecs, k_i, k_ef, f_a, aoffset, voffset - Model parameters, 1 per voxel of
%      the Gadoxetate model
%
%      model_signals - N_vox x N_times array of modelled concentration
%      time-series
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 29-Mar-2017
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 
% Unpack the arguments:
args = u_packargs(varargin, 0, ...
    'FA', 20,...
    'TR', 4,...
    'relax_coeff', 4.5e-3,...
    'Ktrans_init', 0.397579,...
    'Vp_init', 0.015379,...
    'Ve_init', 0.166493,...
    'offset_init', 0,...
    'num_baseline', 5,...
    'num_itr', 20,...
    'opt_display', 'final',...
    'plot_rows', 3,...
    'plot_cols', 5,...
    'fig_dir', [],...
    'debug', false);
clear varargin;

%Set up params structure from default/input args
params.FA = args.FA;
params.TR = args.TR;
params.relax_coeff = args.relax_coeff;
params.num_baseline = args.num_baseline;
params.debug = args.debug;

%Expand up the initial values if only single value shave been given
num_voxels = size(dyn_signals, 2);
if length(args.Ktrans_init) == 1
    args.Ktrans_init = ones(num_voxels,1)*args.Ktrans_init;
end
if length(args.Vp_init) == 1
    args.Vp_init = ones(num_voxels,1)*args.Vp_init;
end
if length(args.Ve_init) == 1
    args.Ve_init = ones(num_voxels,1)*args.Ve_init;
end
if length(args.offset_init) == 1
    args.offset_init = ones(num_voxels,1)*args.offset_init;
end

%Get mean signal at baseline (pre-bolus)
baseline_S0 = mean(dyn_signals(:,1:args.num_baseline),2);

%Set options for optimisation
% opt_options = optimoptions('fmincon',...
%     'Algorithm','interior-point',...
%     'Display', args.opt_display,...
%     'MaxIter', args.num_itr);

%Create containers for fitted model params
offset = zeros(num_voxels, 1);

for i_vox = 1:num_voxels
    %Set T0 for this voxel
    params.T1 = T1(i_vox);
    params.S0 = baseline_S0(i_vox);
    
    %Set initial values for extended Kety model params
    x_init = args.offset_init(i_vox);
    
    params.Ktrans = args.Ktrans_init(i_vox);
    params.Vp = args.Vp_init(i_vox);
    params.Ve = args.Ve_init(i_vox);
    
    %Make sure initial values for Ve are valid
    if ~params.Ve
        params.Ve = mean(args.Ve_init);
    end
    
    if args.debug
        if rem(i_vox, args.plot_rows*args.plot_cols) == 1
            if ~isempty(args.fig_dir) && i_vox > 1
                exportfig([args.fig_dir 'fig' zerostr((i_vox-1)/24, 3) '.png']);
            end
            figure; plot_num = 1;
        else
            plot_num = plot_num + 1;
        end
        subplot(args.plot_rows, args.plot_cols, plot_num);
        plot(dyn_signals(i_vox,:));
        axis([0 75 0 400]);
        hold all;
        [model_signals] = ...
            fit_model(dyn_times, aif, x_init(1), params);
        plot(model_signals, 'g', 'linewidth', 2);
    end
    
    %objective function, optimising
    obj_fun = @(x)compute_model_SSE(x,...
        dyn_signals(i_vox,:), dyn_times, aif,... %fixed variables in the objective function
        params); %auxilliary variables required to compute objective

    offset(i_vox) = fmincon(obj_fun, x_init, [], [], [], [], 0, 1);
    
    if args.debug
        [model_signals] = ...
            fit_model(dyn_times, aif, offset(i_vox), params);
        plot(model_signals, 'r', 'linewidth', 2);
    end
end

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
function [model_signals] = fit_model(dyn_times, aif, offset, params)
[model_concentrations] = extended_kety_model(dyn_times, aif, params.Ktrans, params.Vp, params.Ve, offset);
[model_signals_unscaled] = concentration_to_signal(model_concentrations, ...
    params.FA, params.TR, params.T1, 1, params.relax_coeff);

model_baseline = mean(model_signals_unscaled(:,1:params.num_baseline),2);
model_signals = params.S0*model_signals_unscaled / model_baseline;

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
function [SSE] = compute_model_SSE(x, dyn_signals, dyn_times, aif, params)

[model_signals] = fit_model(dyn_times, aif, x, params);

SSE = sum((dyn_signals - model_signals).^2,2);

if params.debug > 1
    plot(model_signals, '--');
end
    





