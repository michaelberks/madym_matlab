function [aic0, error0, aic1, error1, aic_linear, error_linear] =...
    fit_baseline_model(dyn_input, varargin)
%FIT_BASELINE_MODEL *Insert a one line summary here*
%   [baseline_aic] = fit_baseline_model(dyn_signals)
%
% Inputs:
%      dyn_input - dynamic time-series either signals or concentration. 
%       If 2D input, assumed time-points on 1st dimension (so dyn_input is
%       n_times x n_vox). If >2D, assumed time is on last dimensions - eg
%       a stack of volumes n_y x n_x x n_z x n_times
%
%
% Outputs:
%      baseline_aic - *Insert description of input variable here*
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 08-Aug-2018
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 
args = u_packargs(varargin, 0, ...
    'T1', [],... 
    'S0',[],...
    'FA', 20,...
    'TR', 4,...
    'input_type', 'signal',...
    'relax_coeff', 8.0e-3,...
    'temporal_noise', [],...
    'num_baseline', 8,...
    'plot_rows', 3,...
    'plot_cols', 5,...
    'fig_dir', [],...
    'save_path', [],...
    'debug', false);
clear varargin;

num_dims = ndims(dyn_input);
if num_dims == 2
    reshape_output = false;
else
    %Move input to first dimension, then make 2D
    reshape_output = true;
    orig_sz = size(dyn_input);
    dyn_input = permute(dyn_input, [num_dims 1:num_dims-1]);
    dyn_input = dyn_input(:,:);
end

%If we have raw signals convert to concentrations
if strcmpi(args.input_type, 'signal')
    dyn_concentrations = signal_to_concentration(...
        dyn_input, args.FA, args.TR, args.T1, args.relax_coeff, ...
        args.num_baseline, args.S0);
else
    dyn_concentrations =  dyn_input;
end
[n_t, n_v] = size(dyn_concentrations);

if isempty(args.temporal_noise)
    temporal_noise = ones(n_t,1);
else
    temporal_noise = args.temporal_noise(:);
end

do_1 = nargout > 2;
do_linear = nargout > 4;

%Get overall mean and error
overall_mean = mean(dyn_concentrations);
error0 = sum(((dyn_concentrations - overall_mean).^2) ./ temporal_noise);
k0 = 1;
aic0 = compute_AIC(error0, n_t, k0, true);

if do_1
    %Get pre and post contrast mean concentration
    baseline_conc = dyn_concentrations(1:args.num_baseline,:);
    post_conc = dyn_concentrations(args.num_baseline+1:end,:);

    baseline_mean = mean(baseline_conc);
    post_mean = mean(post_conc);
    
    baseline_error = sum(((baseline_conc - baseline_mean).^2) ...
    ./ temporal_noise(1:args.num_baseline));
    post_error = sum(((post_conc - post_mean).^2) ...
        ./ temporal_noise(args.num_baseline+1:end));
    error1 = baseline_error + post_error;
    
    k1 = 2;
    aic1 = compute_AIC(error1, n_t, k1, true);
end

if do_linear
    %Fit a straight line to mean concentration
    linear_fit = zeros(n_t, n_v);
    x = (1:n_t)';
    for i_v = 1:n_v
        y = dyn_concentrations(:,n_v);
        if any(isnan(y))
            continue;
        end
        p = polyfit(x, y, 1);
        linear_fit(:,n_v) = p(1)*x + p(2);
    end
    error_linear = sum(((dyn_concentrations - linear_fit).^2) ./ temporal_noise);

    k_linear = 3;
    aic_linear = compute_AIC(error_linear, n_t, k_linear, true);
end

if reshape_output
    output_sz = orig_sz(1:num_dims-1);
    
    aic0 = reshape(aic0, output_sz);
    error0 = reshape(error0, output_sz);
    
    if do_1
        aic1 = reshape(aic1, output_sz);
        error1 = reshape(error1, output_sz);
    end
    if do_linear
        aic_linear = reshape(aic_linear, output_sz);
        error_linear = reshape(error_linear, output_sz);
    end
end



