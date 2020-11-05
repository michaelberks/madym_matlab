function [F_p, v_ecs, k_i, k_ef, f_a, aoffset, voffset] =...
    fit_gadexotate_model(dyn_signals, dyn_times, aif, pif, varargin)
%FIT_GADEXOTATE_MODEL *Insert a one line summary here*
%   [model_signals] = extended_kety_model(dyn_times, aif, Ktrans, Vp, Ve)
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
%     'F_p_init', 1, Initial values for parameters, either scalar or 1 per vox
%     'v_ecs_init', 0.2,...
%     'k_i_init', 0.2,...
%     'k_ef_init', 0.015,...
%     'f_a_init', 0.1,...
%     'aoffset_init', 0.0,...
%     'voffset_init', 0.0,...', 0.0,...
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
    'T1', [],...
    'S0',[],...
    'FA', 20,...
    'TR', 4,...
    'fit_to', 'signal',...
    'input', 'signal',...
    'relax_coeff', 8.0e-3,...
    'F_p_init', 1,...
    'v_ecs_init', 0.2,...
    'k_i_init', 0.2,...
    'k_ef_init', 0.015,...
    'f_a_init', 0.1,...
    'aoffset_init', 0.0,...
    'voffset_init', 0.0,...
    'lower_bounds', [-5 -10 -5 -2 -inf 0 -0.5],...
    'upper_bounds', [1 0 3 -0.3 1 0.5 0.5],...
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
if fit_to_signal
    params.FA = args.FA;
    params.TR = args.TR;
    params.relax_coeff = args.relax_coeff;
    params.num_baseline = args.num_baseline;
    params.debug = args.debug;
end

%Expand up the initial values if only single value shave been given
num_voxels = size(dyn_signals, 2);
if length(args.F_p_init) == 1
    args.F_p_init = ones(num_voxels,1)*args.F_p_init;
end
if length(args.v_ecs_init) == 1
    args.v_ecs_init = ones(num_voxels,1)*args.v_ecs_init;
end
if length(args.k_i_init) == 1
    args.k_i_init = ones(num_voxels,1)*args.k_i_init;
end
if length(args.k_ef_init) == 1
    args.k_ef_init = ones(num_voxels,1)*args.k_ef_init;
end
if length(args.f_a_init) == 1
    args.f_a_init = ones(num_voxels,1)*args.f_a_init;
end
if length(args.aoffset_init) == 1
    args.aoffset_init = ones(num_voxels,1)*args.aoffset_init;
end
if length(args.voffset_init) == 1
    args.voffset_init = ones(num_voxels,1)*args.voffset_init;
end

%Get mean signal at baseline (pre-bolus)
if isempty(args.S0)
    baseline_S0 = mean(dyn_signals(1:args.num_baseline,:));
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
F_p = zeros(num_voxels, 1);
v_ecs = zeros(num_voxels, 1);
k_i = zeros(num_voxels, 1);
k_ef = zeros(num_voxels, 1);
f_a = zeros(num_voxels, 1);
aoffset = zeros(num_voxels, 1);
voffset = zeros(num_voxels, 1);

%If we're fitting in concentration space, pre-convert signals to
%concentrations here
if ~fit_to_signal && strcmpi('input', 'signal')
    dyn_concentrations = signal_to_concentration(...
        dyn_signals, args.FA, args.TR, args.T1, args.relax_coeff, ...
        args.num_baseline, args.S0);
else
    dyn_concentrations = dyn_signals;
end

for i_vox = 1:num_voxels
    if fit_to_signal
        %Set T0 for this voxel
        params.T1 = args.T1(i_vox);
        params.S0 = baseline_S0(i_vox);
    end
    
    %Set initial values for extended Kety model params
    x_init = [
        log10(args.F_p_init(i_vox))
        log10(args.v_ecs_init(i_vox))
        log10(args.k_i_init(i_vox))
        log10(args.k_ef_init(i_vox))
        log10(args.f_a_init(i_vox))
        args.aoffset_init(i_vox)
        args.voffset_init(i_vox)];
    
    %Make sure initial values are valid... check what values are valid?
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
        plot(dyn_signals(:,i_vox));
        axis([0 75 0 400]);
        hold all;
        %[model_signals] = ...
        %    fit_model_to_signal(dyn_times, aif, x_init(1), x_init(2), x_init(3), x_init(4), params);
        plot(model_signals, 'g', 'linewidth', 2);
    end
    
    if fit_to_signal
        %objective function, optimising model params in signal space 
        obj_fun = @(x)compute_signal_SSE(x,...
            dyn_signals(i_vox,:), dyn_times, aif, pif,... %fixed variables in the objective function
            params); %auxilliary variables required to compute objective
    else
        %objective function, optimising model params in concentration space 
        obj_fun = @(x)compute_concentration_SSE(x,...
            dyn_concentrations(i_vox,:), dyn_times, aif, pif);%fixed variables in the objective function
    end

    x_fit = fmincon(obj_fun, x_init, [], [], [], [], args.lower_bounds, args.upper_bounds, [],...
        opt_options); 
    
    F_p(i_vox) = 10^x_fit(1);
    v_ecs(i_vox) = 10^x_fit(2);
    k_i(i_vox) = 10^x_fit(3);
    k_ef(i_vox) = 10^x_fit(4);
    f_a(i_vox) = 10^x_fit(5);
    aoffset(i_vox) = x_fit(6);
    voffset(i_vox) = x_fit(7);
    
    if args.debug
        [model_signals] = ...
            fit_model_to_signal(dyn_times, aif, pif, ...
            F_p(i_vox), v_ecs(i_vox), k_i(i_vox), k_ef(i_vox), f_a(i_vox),...
            aoffset(i_vox), voffset(i_vox), params);
        plot(model_signals, 'r', 'linewidth', 2);
    end
end

if ~isempty(args.save_path)
    create_folder(fileparts(args.save_path))
    save(args.save_path, ...
        'F_p', 'v_ecs', 'k_i', 'k_ef', 'f_a', 'aoffset', 'voffset');
end

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
function [model_signals] = fit_model_to_signal(...
    dyn_times, aif, pif, F_p, v_ecs, k_i, k_ef, f_a, aoffset, voffset, params)
[model_concentrations] = gadoxetate_model(...
    F_p, v_ecs, k_i, k_ef, f_a, aoffset, aif, pif, dyn_times, [], voffset);
    
[model_signals] = concentration_to_signal(model_concentrations, ...
    params.FA, params.TR, params.T1, params.S0, params.relax_coeff, params.num_baseline);
% [model_signals_unscaled] = concentration_to_signal(model_concentrations, ...
%     params.FA, params.TR, params.T1, 1, params.relax_coeff);
% 
% model_baseline = mean(model_signals_unscaled(1:params.num_baseline,:));
% model_signals = params.S0*model_signals_unscaled / model_baseline;

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
function [SSE] = compute_signal_SSE(x, dyn_signals, dyn_times, aif, pif, params)

F_p = 10^x(1);
v_ecs = 10^x(2);
k_i = 10^x(3);
k_ef = 10^x(4);
f_a = 10^x(5);
aoffset = x(6);
voffset = x(7);

[model_signals] = fit_model_to_signal(...
    dyn_times, aif, pif, F_p, v_ecs, k_i, k_ef, f_a, aoffset, voffset, params);

SSE = sum((dyn_signals(1:end,:) - model_signals(1:end,:)).^2);

if params.debug > 1
    plot(model_signals, '--');
end

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
function [SSE] = compute_concentration_SSE(x, dyn_concentrations, dyn_times, aif, pif)

F_p = 10^x(1);
v_ecs = 10^x(2);
k_i = 10^x(3);
k_ef = 10^x(4);
f_a = 10^x(5);
aoffset = x(6);
voffset = x(7);

[model_concentrations] = gadoxetate_model(...
    F_p, v_ecs, k_i, k_ef, f_a, aoffset, aif, pif, dyn_times, [], voffset);
SSE = sum((dyn_concentrations(1:end,:) - model_concentrations(1:end,:)).^2);

    





