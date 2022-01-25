function [model_params, model_fit, iauc, error_codes, Ct_m, Ct_s] =...
    run_madym_lite(model, input_data, varargin)
%RUN_MADYM_LITE wrapper function to call C++ tool Madym-lite. Fits
%   tracer-kinetic models to DCE time-series, returning the model
%   parameters and modelled concentration time-series
%   [model_params, model_fit, error_codes, Ct_m, Ct_s] = ...
%       run_madym_lite(model, input_data, varargin)
%
% RUN_MADYM_LITE uses the U_PACKARGS interface function
% and so arguments to the function may be passed as name-value pairs
% or as a struct with fields with corresponding names. The field names
% are defined as:
%
% Mandatory Arguments:
%       model (str) - Model to fit, specified by its name in CAPITALS,
%            see notes for options
%
%       input_data 2D array (Nsamples x Ntimes) - Either signals or 
%            concentration, 1 time-series per row
%
% Optional Arguments: This wrapper uses the upackargs interface to allow
%       optional arguments to be entered as name/value pairs. See below for
%       a full description of all options.
%
% Outputs:
%      model_params (2D array, Nsamples x Nparams) - each row contains 
%       the estimated model parameters for the corresponding time-series 
%       in the input_data. The number of columns depends on the model 
%       parameters. The parameters are returned in the same order as the 
%       defined in Madym. See notes for details for each model.
%
%      model_fit (1D array, Nsamples x 1) - sum-of-squared model 
%       residuals for each time-series
%
%      error_codes (2D array, Nsamples x 2) - error codes returned by MaDym
%       for fitting each sample. 0 implies no errors or warnings. For all
%       non-zero values refer to Madym documentation for details
%
%      Ct_m (2D array, Nsamples x Ntimes) - the modelled
%       concentration time-series for each input
%
%      Ct_s (2D array, Nsamples x Ntimes) - the signal-derived
%       concentration time-series for each input
%
% Examples:
%   Fitting to concentration time-series. If using a population AIF, you
%   must supply a vector of dynamic times. A population AIF (PIF) is used
%   if the aif_name (pif_name) option is left empty.
%   [model_params, model_fit, error_codes, Ct_m] = ...
%       run_madym_lite("2CXM", Ct_s, 'dyn_times', t)
%
%   Fitting to concentration time-series using a patient specific AIF. The
%   AIF should be defined in a text file with two columns containing the
%   dynamic times and associated AIF value at each times respectively. Pass
%   the full filepath as input
%   [model_params, model_fit, error_codes, Ct_m] = ...
%       run_madym_lite("2CXM", Ct_s, 'aif_name', 'C:\DCE_data\pt_AIF.txt')
%
%   Fitting to signals - Set input_Ct to false and use options to supply
%   T1 values (and TR, FA, relax_coeff etc) to convert signals to
%   concentration.
%   [model_params, model_fit, error_codes, Ct_m, Ct_s] = ...
%       run_madym_lite("2CXM", dyn_signals, 'dyn_times', t,...
%           'input_Ct', 0, 'T1', T1_vals, 'TR', TR, 'FA', FA)
%
%   Fixing values in a model - eg to fit a TM instead of ETM, set Vp (the
%   3rd parameter in the ETM to 0)
%   [model_params, model_fit, error_codes, Ct_m] = ...
%       run_madym_lite("ETM", Ct_s, 'dyn_times', t,...
%           'fixed_params', 3, 'fixed_values', 0.0)
%
% Notes:
%   Tracer-kinetic models:
% 
%   All models available in the main MaDym and MaDym-Lite C++ tools are
%   available to fit. Currently these are:
% 
%   See the madym_cxx project wiki for more details:
% https://gitlab.com/manchester_qbi/manchester_qbi_public/madym_cxx/-/wikis/dce_models
% 
%   Run: system([local_madym_root 'madym_DCE_lite --help']); to see full set of
%   input options to C++ tool
% 
% See also: RUN_MADYM, TRUN_MADYM_T1
%
% Created: 20-Feb-2019
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 

if ~nargin
    run_test();
    return;
end

% Unpack the arguments:
args = u_packargs(varargin, 1, ... 
    'cmd_exe', [local_madym_root 'madym_DCE_lite'],...
    'output_dir', [], ...Output path, will use temp dir if empty;
  	'output_name', 'madym_analysis.dat', ... Name of output file
    'dyn_times', [], ... Time associated with each dynamic signal (in mins), must be supplied if using population AIF
    'input_Ct', true, ... Flag specifying input dynamic sequence are concentration (not signal) maps
    'output_Ct_sig', nargout > 5,... Flag requesting concentration (derived from signal) are saved to output
    'output_Ct_mod', nargout > 4, ...Flag requesting modelled concentration maps are saved to output
    'no_optimise', NaN, ...Flag to switch off optimising, will just fit initial parameters values for model
    'B1_correction', NaN, ... Apply B1 correction
... The below are all only required if we're converting from signals
    'T1', [], ...Baseline T1 values (in ms)
    'M0', [], ...Baseline M0 values, required if not using ratio method
 	'TR', NaN, ... TR of dynamic series (in ms), must be >0 if converting signals
    'FA', NaN, ... Flip angle of dynamic series (degrees), must be set if converting signals
    'r1_const', NaN, ...Relaxivity constant of concentration in tissue (in ms)
    'M0_ratio', true, ... Flag to use ratio method to scale signal instead of supplying M0
...
    'injection_image', NaN, ... Injection image
    'hct', NaN, ...Haematocrit correction
  	'first_image', NaN, ... First image used to compute model fit
    'last_image', NaN, ... Last image used to compute model fit
    'dose', NaN, ...Concentration dose (mmole/kg)
...
    'aif_name', [], ... Path to precomputed AIF if not using population AIF
  	'pif_name', [], ... Path to precomputed PIF if not deriving from AIF
...
    'IAUC_times', [60 90 120],... "_times (in s) at which to compute IAUC values
    'init_params', [],... Initial values for model parameters to be optimised, either as single vector, or 2D array NSamples x N_params
    'fixed_params', [],...Parameters fixed to their initial values (ie not optimised)
    'fixed_values', [],... _values for fixed parameters (overrides default initial parameter values)"
    'upper_bounds', [],...Upper bounds for each parameter during optimisation
    'lower_bounds', [],...Lower bounds for each parameter during optimisation
    'relative_limit_params', [],...Parameters with relative limits on their optimisation bounds
    'relative_limit_values', [], ..._values for relative bounds, sets lower/upper bound as init param -/+ relative limit
...
    'dyn_noise_values', [],...Varying temporal noise in model fit
  	'max_iter', NaN,... Maximum number of iterations in model fit
    'opt_type', '',... Type of optimisation to run
    'test_enhancement', NaN, ...Set test-for-enhancement flag
    'quiet', NaN,... Suppress output to stdout
    'working_directory', '',...Sets the current working directory for the system call, allows setting relative input paths for data
    'dummy_run', false ...Don't run any thing, just print the cmd we'll run to inspect
    );
clear varargin;

%Get size number of dynamic values - if we've been given a single column
%vector, transpose into a row. Should n_voxels x n_dyns
if size(input_data, 2) == 1
    input_data = input_data';
end
[nSamples, nDyns] = size(input_data);

%If we're converting from signal to concentration, append T1 (and M0 if not
%using ratio method) to pinput data
if ~args.input_Ct
    input_data = [input_data args.T1(:)];
    
    if ~args.M0_ratio
        input_data = [input_data args.M0(:)];
    end
end

%Check for bad samples, these can screw up Madym as the lite version
%of the software doesn't do the full range of error checks Madym proper
%does. So chuck them out now and warn the user
discard_samples = any(isnan(input_data), 2);
if any(discard_samples)
    warning(['Samples with NaN values found,'...
        'these will be set to zero for model-fitting']);
    input_data(discard_samples,:) = 0;
end

discard_samples = any(~isfinite(input_data), 2);
if any(discard_samples)
    warning(['Samples with NaN values found,'...
        'these will be set to zero for model-fitting']);
    input_data(discard_samples,:) = 0;
end

%Get a name for the temporary file we'll write input data to (we'll hold
%off writing anything until we know this isn't a dummy run
input_file = tempname;
    
deleteOutput = 0;
if isempty(args.output_dir)
    args.output_dir = tempdir;
    deleteOutput = 1;
elseif args.output_dir(end) ~= '\' && args.output_dir(end) ~= '/'
    args.output_dir = [args.output_dir '/'];
end
fullOutPath = [args.output_dir model '_' args.output_name];

%Set up initial string using all the values we can directly from the Matlab
%args structure
cmd = sprintf(...
    '%s -m %s --data %s -n %d -o %s -O %s',...
    args.cmd_exe,...
    model,...
    input_file,...
    nDyns,...
    args.output_dir,...
    args.output_name);

%Set the working dir
cmd = add_option('string', cmd, '--cwd', args.working_directory);

%Now set any args that require option inputs
if args.input_Ct
    cmd = sprintf('%s --Ct', cmd);
    
else %Below are only needed if input is signals
    cmd = add_option('float', cmd, '--TR', args.TR);

    cmd = add_option('float', cmd, '--FA', args.FA);

    cmd = add_option('float', cmd, '--r1', args.r1_const);

    cmd = add_option('bool', cmd, '--M0_ratio', args.M0_ratio);
end

cmd = add_option('float', cmd, '-D', args.dose);

cmd = add_option('float', cmd, '-H', args.hct);

cmd = add_option('int', cmd, '-i', args.injection_image);

cmd = add_option('int', cmd, '--first', args.first_image);

cmd = add_option('int', cmd, '--last', args.last_image);

cmd = add_option('bool', cmd, '--Ct_sig', args.output_Ct_sig);

cmd = add_option('bool', cmd, '--Ct_mod', args.output_Ct_mod);

cmd = add_option('bool', cmd, '--no_opt', args.no_optimise);

cmd = add_option('bool', cmd, '--B1_correction', args.B1_correction);

cmd = add_option('int', cmd, '--max_iter', args.max_iter);

cmd = add_option('string', cmd, '--opt_type', args.opt_type);

cmd = add_option('bool', cmd, '--test_enh', args.test_enhancement);

cmd = add_option('bool', cmd, '--quiet', args.quiet);

cmd = add_option('string', cmd, '--aif', args.aif_name);

cmd = add_option('string', cmd, '--pif', args.pif_name);

if ~isempty(args.dyn_times)
    %Get a name for the temporary file we'll write times to (we'll hold
    %off writing anything until we know this isn't a dummy run
    dyn_times_file = tempname;
    cmd = add_option('string', cmd, '-t', dyn_times_file);
end

if ~isempty(args.dyn_noise_values)
    %Get a name for the temporary file we'll write noise to (we'll hold
    %off writing anything until we know this isn't a dummy run
    dynNoise_file = tempname;
    cmd = add_option('string', cmd, '--dyn_noise', dynNoise_file);
end

cmd = add_option('float_list', cmd, '--iauc', args.IAUC_times);

load_params = false;
if ~isempty(args.init_params)
    if nSamples > 1 && size(args.init_params,1) == nSamples
        input_params_file = tempname;
        load_params = true;
        cmd = add_option('string', cmd, '--init_params_file', input_params_file);
    else
        cmd = add_option('float_list', cmd, '--init_params', args.init_params);
    end
end

cmd = add_option('int_list', cmd, '--fixed_params', args.fixed_params);
cmd = add_option('float_list', cmd, '--fixed_values', args.fixed_values);

cmd = add_option('float_list', cmd, '--upper_bounds', args.upper_bounds);
cmd = add_option('float_list', cmd, '--lower_bounds', args.lower_bounds);

cmd = add_option('int_list', cmd, '--relative_limit_params', ...
    args.relative_limit_params);
cmd = add_option('float_list', cmd, '--relative_limit_values', ...
    args.relative_limit_values);

if args.dummy_run
    %Don't actually run anything, just print the command
    fprintf('%s\n', cmd);
    
    model_params = [];
    model_fit = [];
    iauc = [];
    error_codes = [];
    Ct_m = [];
    Ct_s = [];
    return;
    
end

%Write input values to a temporary file
fid = fopen(input_file, 'wt');
for i_row = input_data'
    fprintf(fid, '%6.5f ', i_row);
    fprintf(fid, '\n');
end
fclose(fid);

%Write input params to a temporary file
if load_params
    discard = ~isfinite(args.init_params) | isnan(args.init_params);
    args.init_params(discard) = 0;
    
    fid = fopen(input_params_file, 'wt');
    for i_row = args.init_params'
        fprintf(fid, '%6.5f ', i_row);
        fprintf(fid, '\n');
    end
    fclose(fid);
end

%Write dynamic times to a temporary file if we need to
if ~isempty(args.dyn_times)
    fid = fopen(dyn_times_file, 'wt');
    fprintf(fid, '%6.5f ', args.dyn_times(:));
    fclose(fid);
end

%Write noise vals to a temporary file if we need to
if ~isempty(args.dyn_noise_values)
    fid = fopen(dynNoise_file, 'wt');
    fprintf(fid, '%6.5f ', args.dyn_noise_values(:));
    fclose(fid);
end

%At last.. we can run the command
system(cmd);

%Now load the output from madym lite and extract data to match this
%functions outputs
outputData = load(fullOutPath);

error_codes = outputData(:,1:2);
model_fit = outputData(:,3);

n_iauc = length(args.IAUC_times);% + (args.IAUC_at_peak > 0);
iauc = outputData(:, 3 + (1:n_iauc));

%Workout which columns hold parameter values
n_cols = size(outputData, 2);

param_col1 = n_iauc + 4;
param_col2 = n_cols;
if args.output_Ct_mod
    if args.output_Ct_sig
        %Have another nDyns of dynamic concentration values
        param_col2 = param_col2 - nDyns;
        Ct_s = outputData(:, param_col2+(1:nDyns));
    end
    
    %After the params we have nDyns of model concentration values
    param_col2 = param_col2 - nDyns;
    Ct_m = outputData(:, param_col2+(1:nDyns));
end
model_params = outputData(:, param_col1:param_col2);


%Finally, tidy up any temporary files
delete(input_file);
if load_params
    delete(input_params_file);
end
if ~isempty(args.dyn_times)
    delete(dyn_times_file);
end
if ~isempty(args.dyn_noise_values)
    delete(dynNoise_file);
end
if deleteOutput
    delete(fullOutPath);
end

%--------------------------------------------------------------------------

%%
%Test function to run if no inputs
function run_test()

%Generate an concentration time-series using the ETM
ktrans = 0.25;
ve = 0.2;
vp = 0.1;
tau = 0;
injection_img = 8;
t = linspace(0, 5, 100);
Ca_t = population_aif(t, injection_img);
C_t = extended_tofts_model(ktrans, ve, vp, tau, Ca_t, t);

%Add some noise and rescale so baseline mean is 0
C_tn = C_t + randn(1,100)/100;
C_tn = C_tn - mean(C_tn(1:injection_img),2);

%Use madym lite to fit this data
[model_params_C, model_fit_C, ~,~,CmC_t] = run_madym_lite(...
    'ETM', C_tn, 'dyn_times', t);

%Convert the concentrations to signals with some notional T1 values and
%refit using signals as input
FA = 20;
TR = 3.5;
T1_0 = 1000;
r1_const = 3.4;
S_t0 = 100;
S_tn = concentration_to_signal(...
    C_tn, FA, TR, T1_0, S_t0, r1_const, injection_img);

[model_params_S, model_fit_S, ~,~,CmS_t,Cm_t] = run_madym_lite(...
    'ETM', S_tn, 'dyn_times', t,...
    'input_Ct', 0,...
    'T1', T1_0,...
    'TR', TR,...
    'FA', FA,...
    'r1_const', r1_const,...
    'injection_image', injection_img);

%Convert the modelled concentrations back to signal space
Sm_t = concentration_to_signal(...
    CmS_t, FA, TR, T1_0, S_t0, r1_const, injection_img);

%Display plots of the fit
figure('Name', 'madym_lite test applied');
subplot(2,2,[1 3]);
plot(t, C_tn); hold all;
plot(t, CmC_t);
legend({'C(t)', 'ETM model fit'});
xlabel('Time (mins)');
ylabel('Voxel concentration');
title(sprintf('Input C(t): Model fit SSE = %4.3f', model_fit_C));

subplot(2,2,2);
plot(t, C_tn); hold all;
plot(t, Cm_t, '--');
plot(t, CmS_t);
legend({'C(t)', 'C(t) (output from MaDym)', 'ETM model fit'});
xlabel('Time (mins)');
ylabel('Voxel concentration');
title(sprintf('Input S_t: Model fit SSE = %4.3f', model_fit_S));

subplot(2,2,4);
plot(t, S_tn); hold all;
plot(t, Sm_t);
legend({'S(t)', 'ETM model fit - converted to signal'});
xlabel('Time (mins)');
ylabel('Voxel signal');
title(sprintf('Input S(t): Signal SSE = %4.3f', sum((S_tn-Sm_t).^2)));

fprintf('Parameter estimation (actual, fitted concentration, fitted signal)\n');
fprintf('Ktrans: (%3.2f, %3.2f, %3.2f)\n', ktrans,...
    model_params_C(1), model_params_S(1));
fprintf('Ve: (%3.2f, %3.2f, %3.2f)\n', ve,...
    model_params_C(2), model_params_S(2));
fprintf('Vp: (%3.2f, %3.2f, %3.2f)\n', vp,...
    model_params_C(3), model_params_S(3));
fprintf('Tau: (%3.2f, %3.2f, %3.2f)\n', tau,...
    model_params_C(4), model_params_S(4));




