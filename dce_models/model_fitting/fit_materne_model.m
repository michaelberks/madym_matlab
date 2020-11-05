function [k1a, k1v, k2, aoffset, voffset, model_signals] = fit_materne_model(...
    dyn_signals, dyn_times, aif, pif, T1, varargin)
%EXTENDED_KETY_MODEL Fit the Materne DCE model to an array of input DCE
%time-series data
%   [model_signals] = fit_materne_model(dyn_times, aif, T1, varargin)
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
%     'k1a_init', 0.2, Initial values for parameters, either scalar or 1 per vox
%     'k1v_init', 0.2,...
%     'k2_init', 0.2,...
%     'aoffset_init', 0.0,...
%     'voffset_init', 0.0,...
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
%
% Outputs:
%
%      k1a, k1v, k2, aoffset, voffset - Model parameters, 1 per voxel of
%      the Materne model
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
    'relax_coeff', 14.0e-3,...
    'k1a_init', 0.2,...
    'k1v_init', 0.2,...
    'k2_init', 0.2,...
    'aoffset_init', 0.0,...
    'voffset_init', 0.0,...
    'lower_bounds', [-10 -10 -10 0 0],...
    'upper_bounds', [1 1 1 0.5 0.5],...
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
if length(args.k1a_init) == 1
    args.k1a_init = ones(num_voxels,1)*args.k1a_init;
end
if length(args.k1v_init) == 1
    args.k1v_init = ones(num_voxels,1)*args.k1v_init;
end
if length(args.k2_init) == 1
    args.k2_init = ones(num_voxels,1)*args.k2_init;
end
if length(args.aoffset_init) == 1
    args.aoffset_init = ones(num_voxels,1)*args.aoffset_init;
end
if length(args.voffset_init) == 1
    args.voffset_init = ones(num_voxels,1)*args.voffset_init;
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
k1a = zeros(num_voxels, 1);
k1v = zeros(num_voxels, 1);
k2 = zeros(num_voxels, 1);
aoffset = zeros(num_voxels, 1);
voffset = zeros(num_voxels, 1);

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
        log10(args.k1a_init(i_vox))
        log10(args.k1v_init(i_vox))
        log10(args.k2_init(i_vox))
        args.aoffset_init(i_vox)
        args.voffset_init(i_vox)];
    
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
            fit_model_to_signal(dyn_times, aif,...
            x_init(1), x_init(2), x_init(3), x_init(4), x_init(5), params);
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
    
    k1a(i_vox) = 10^x_fit(1);
    k1v(i_vox) = 10^x_fit(2);
    k2(i_vox) = 10^x_fit(3);
    aoffset(i_vox) = x_fit(4);
    voffset(i_vox) = x_fit(5);
    
    [model_signals] = ...
        fit_model_to_signal(dyn_times, aif, pif,...
        k1a(i_vox), k1v(i_vox), k2(i_vox), aoffset(i_vox), voffset(i_vox), params);
    
    if args.debug    
        plot(model_signals, 'r', 'linewidth', 2);
    end
end

if ~isempty(args.save_path)
    create_folder(fileparts(args.save_path))
    save(args.save_path, 'k1a', 'k1v', 'k2', 'aoffset', 'voffset');
end

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
function [model_signals] = fit_model_to_signal(...
    dyn_times, aif, pif, k1a, k1v, k2, aoffset, voffset, params)
[model_concentrations] = materne_model(...
    k1a, k1v, k2, aoffset, aif, pif, dyn_times, 0, voffset);
[model_signals] = concentration_to_signal(model_concentrations, ...
    params.FA, params.TR, params.T1, params.S0, params.relax_coeff, params.num_baseline);


%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
function [SSE] = compute_signal_SSE(x, dyn_signals, dyn_times, aif, pif, params)

k1a = 10^x(1);
k1v = 10^x(2);
k2 = 10^x(3);
aoffset = x(4);
voffset = x(5);

[model_signals] = fit_model_to_signal(...
    dyn_times, aif, pif, k1a, k1v, k2, aoffset, voffset, params);

SSE = sum((dyn_signals - model_signals).^2,2);%params.num_baseline+

if params.debug > 1
    plot(model_signals, '--');
end

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
function [SSE] = compute_concentration_SSE(x, dyn_concentrations, dyn_times, aif, pif)

k1a = 10^x(1);
k1v = 10^x(2);
k2 = 10^x(3);
aoffset = x(4);
voffset = x(5);

[model_concentrations] = materne_model(...
    k1a, k1v, k2, aoffset, aif, pif, dyn_times, 0, voffset);
SSE = sum((dyn_concentrations - model_concentrations).^2,2);

    





