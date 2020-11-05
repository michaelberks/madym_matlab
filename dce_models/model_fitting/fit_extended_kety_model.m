function [Ktrans, Ve, Vp, offset] = fit_extended_kety_model(dyn_signals, dyn_times, aif, T1, varargin)
%FIT_EXTENDED_KETY_MODEL Fit the extended-Kety (Tofts) model to dynamic
%time-series
%   [Ktrans, Ve, Vp, offset] = fit_extended_kety_model(dyn_signals, dyn_times, aif, T1, varargin)
%
% Inputs:
%      dyn_signals - N_vox x N_times array DCE time-series data
%
%      dyn_times - real world times associated with each time-point
%
%      aif - arterial input function
%
%      pif - hepatic portal vein input function
%
%      T1 - baseline T1 for each voxel
%
%      Optional inputs supplied as name/value pairs:
%
%     'S0',[], baseline S0 if not using ratio method
%     'FA', 20, FA of dynamic scans
%     'TR', 4, TR of dynamic scans
%     'fit_to', 'signal', fit to signal or concentration?
%     'relax_coeff', 14.0e-3, tissue relaxivity constant
%     'Ktrans_init', 0.397579, Initial values for parameters, either scalar or 1 per vox
%     'Vp_init', 0.015379,...
%     'Ve_init', 0.166493,...
%     'offset_init', 0.0,...
%     'lower_bounds', [-10 -10 -10 0 0], Upper/lower bounds on parameters
%     'upper_bounds', [1 1 1 0.5 0.5],...
%     'num_baseline', 5, Number of times to compute baseline mean signal
%     'num_itr', 20, Max number of iterations
%     'opt_display', 'final', Display option during optimisation
%     'plot_rows', 3, Number of rows/cols in figure grid for displaying model fits
%     'plot_cols', 5,...
%     'fig_dir', [], If not empty, saves visual output to specfied folder
%     'save_path', [], If not empty, saves fitted parameters to specified file
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
    'S0',[],...
    'FA', 20,...
    'TR', 4,...
    'fit_to', 'signal',...
    'relax_coeff', 4.5e-3,...
    'Ktrans_init', 0.397579,...
    'Vp_init', 0.015379,...
    'Ve_init', 0.166493,...
    'offset_init', 0.0,...
    'lower_bounds', [-10 -10 -10 0],...
    'upper_bounds', [1 1 1 0.3],...
    'num_baseline', 5,...
    'num_itr', 20,...
    'opt_display', 'final',...
    'plot_rows', 3,...
    'plot_cols', 5,...
    'fig_dir', [],...
    'save_path', [],...
    'debug', false);
clear varargin;

%Check which space we're fitting in
if strcmpi(args.fit_to, 'signal')
    fit_to_signal = true;
    
elseif strcmpi(args.fit_to, 'concentration')
    fit_to_signal = false;
    
else
    warning([args.fit_to ' not valid for fit_to. Must be either signal or concentration. Using fit_to signal']);
    fit_to_signal = true;
end

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
if isempty(args.S0)
    baseline_S0 = mean(dyn_signals(:,1:args.num_baseline),2);
    params.num_baseline = args.num_baseline;
else
    baseline_S0 = args.S0;
    params.num_baseline = 0;
end

% %Set options for optimisation
opt_options = optimoptions('fmincon',...
    'Display', args.opt_display,...
    'MaxIter', args.num_itr);

%Create containers for fitted model params
Ktrans = zeros(num_voxels, 1);
Vp = zeros(num_voxels, 1);
Ve = zeros(num_voxels, 1);
offset = zeros(num_voxels, 1);

%If we're fitting in concentration space, pre-convert signals to
%concentrations here
if ~fit_to_signal
    dyn_concentrations = ...
        signal_to_concentration(dyn_signals, args.FA, args.TR, T1, args.relax_coeff, args.num_baseline);
end

for i_vox = 1:num_voxels
    %Set T0 for this voxel
    params.T1 = T1(i_vox);
    params.S0 = baseline_S0(i_vox);
    
    %Set initial values for extended Kety model params
    x_init = [
        log10(args.Ktrans_init(i_vox))
        log10(args.Vp_init(i_vox))
        log10(args.Ve_init(i_vox))
        args.offset_init(i_vox)];
    
    %Make sure initial values are valid (Vp=0 and offset=0 both allowed, I'm assuming on non-negative values given)
    zero_vals = ~isfinite(x_init);
    x_init(zero_vals) = args.lower_bounds(zero_vals);
    
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
            fit_model_to_signal(dyn_times, aif, x_init(1), x_init(2), x_init(3), x_init(4), params);
        plot(model_signals, 'g', 'linewidth', 2);
    end
    
    if fit_to_signal
        %objective function, optimising model params in signal space 
        obj_fun = @(x)compute_signal_SSE(x,...
            dyn_signals(i_vox,:), dyn_times, aif,... %fixed variables in the objective function
            params); %auxilliary variables required to compute objective
    else
        %objective function, optimising model params in concentration space 
        obj_fun = @(x)compute_concentration_SSE(x,...
            dyn_concentrations(i_vox,:), dyn_times, aif);%fixed variables in the objective function
    end

    x_fit = fmincon(obj_fun, x_init, [], [], [], [], args.lower_bounds, args.upper_bounds, [],...
        opt_options); 
    
    Ktrans(i_vox) = 10^x_fit(1);
    Vp(i_vox) = 10^x_fit(2);
    Ve(i_vox) = 10^x_fit(3);
    offset(i_vox) = x_fit(4);
    
    if args.debug
        [model_signals] = ...
            fit_model_to_signal(dyn_times, aif, Ktrans(i_vox), Vp(i_vox), Ve(i_vox), offset(i_vox), params);
        plot(model_signals, 'r', 'linewidth', 2);
    end
end

if ~isempty(args.save_path)
    create_folder(fileparts(args.save_path))
    save(args.save_path, 'Ktrans', 'Ve', 'Vp', 'offset');
end

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
function [model_signals] = fit_model_to_signal(dyn_times, aif, Ktrans, Vp, Ve, offset, params)
[model_concentrations] = extended_kety_model(dyn_times, aif, Ktrans, Vp, Ve, offset);
[model_signals] = concentration_to_signal(model_concentrations, ...
    params.FA, params.TR, params.T1, params.S0, params.relax_coeff, params.num_baseline);
% [model_signals_unscaled] = concentration_to_signal(model_concentrations, ...
%     params.FA, params.TR, params.T1, 1, params.relax_coeff);
% 
% model_baseline = mean(model_signals_unscaled(1:params.num_baseline,:));
% model_signals = params.S0*model_signals_unscaled / model_baseline;

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
function [SSE] = compute_signal_SSE(x, dyn_signals, dyn_times, aif, params)

Ktrans = 10^x(1);
Vp = 10^x(2);
Ve = 10^x(3);
offset = x(4);

[model_signals] = fit_model_to_signal(dyn_times, aif, Ktrans, Vp, Ve, offset, params);

SSE = sum((dyn_signals(1:end,:) - model_signals(1:end,:)).^2);%params.num_baseline+

if params.debug > 1
    plot(model_signals, '--');
end

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
function [SSE] = compute_concentration_SSE(x, dyn_concentrations, dyn_times, aif)

Ktrans = 10^x(1);
Vp = 10^x(2);
Ve = 10^x(3);
offset = x(4);

[model_concentrations] = extended_kety_model(dyn_times, aif, Ktrans, Vp, Ve, offset);
SSE = sum((dyn_concentrations - model_concentrations).^2,2);

    





