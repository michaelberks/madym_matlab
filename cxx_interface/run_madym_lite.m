function [model_params, model_fit, iauc, error_codes, model_conc, dyn_conc] =...
    run_madym_lite(model, input_data, varargin)
%RUN_MADYM_LITE wrapper function to call C++ tool Madym-lite. Fits
%   tracer-kinetic models to DCE time-series, returning the model
%   parameters and modelled concentration time-series
%   [model_params, model_fit, error_codes, model_conc, dyn_conc] = ...
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
%      model_conc (2D array, Nsamples x Ntimes) - the modelled
%       concentration time-series for each input
%
%      dyn_conc (2D array, Nsamples x 2) - the signal-derived
%       concentration time-series for each input
%
% Examples:
%   Fitting to concentration time-series. If using a population AIF, you
%   must supply a vector of dynamic times. A population AIF (PIF) is used
%   if the aif_name (pif_name) option is left empty.
%   [model_params, model_fit, error_codes, model_conc] = ...
%       run_madym_lite("2CXM", dyn_conc, 'dyn_times', t)
%
%   Fitting to concentration time-series using a patient specific AIF. The
%   AIF should be defined in a text file with two columns containing the
%   dynamic times and associated AIF value at each times respectively. Pass
%   the full filepath as input
%   [model_params, model_fit, error_codes, model_conc] = ...
%       run_madym_lite("2CXM", dyn_conc, 'aif_name', 'C:\DCE_data\pt_AIF.txt')
%
%   Fitting to signals - Set input_Ct to false and use options to supply
%   T1 values (and TR, FA, relax_coeff etc) to convert signals to
%   concentration.
%   [model_params, model_fit, error_codes, model_conc, dyn_conc] = ...
%       run_madym_lite("2CXM", dyn_signals, 'dyn_times', t,...
%           'input_Ct', 0, 'T1', T1_vals, 'TR', TR, 'FA', FA)
%
%   Fixing values in a model - eg to fit a TM instead of ETM, set Vp (the
%   3rd parameter in the ETM to 0)
%   [model_params, model_fit, error_codes, model_conc] = ...
%       run_madym_lite("ETM", dyn_conc, 'dyn_times', t,...
%           'fixed_params', 3, 'fixed_values', 0.0)
%
% Notes:
%   Tracer-kinetic models:
%
%   All models available in the main MaDym and MaDym-Lite C++ tools are
%   available to fit. Currently these are:
% 
%   "ETM"
%   Extended-Tofts model. Requires single input AIF. 
%   Outputs 4 params (= initial values for optimisation):
%   {Ktrans=0.2, Ve=0.2, Vp=0.2, tau_a=0.0*} *Arterial offset delay
%   To run a standard Tofts model, set Vp = 0 by using options 
%   'fixed_params', [3], 'fixed_values', [0.0],...
% 
%   "DIETM"
%   Extended-Tofts model with dual-input supply. Requires both AIF and PIF.
%   Outputs 6 parameters:
%   {Ktrans=0.2, Ve=0.2, Vp=0.2, fa=0.5, tau_a=0.0, tau_v=0.0*} *Venous delay 
% 
%   "GADOXETATE"
%   Leo's model for gaodxetate contrast in the liver. Requires both AIF and PIF.
%   Outputs 7 parameters:
%   {Fp=0.6, ve=0.2, ki=0.2, kef=0.1, fa=0.5, tau_a=0.025, tau_v=0}
%
%   "MATERNE"
%   Dual-input single compartment model
%   {Fp=0.6, fa=0.5, k2=1.0, tau_a=0.025, tau_v=0}
% 
%   "2CXM"
%   2-compartment exchange model. Single AIF input.
%   Outputs 5 parameters:
%   { Fp=0.6, PS=0.2, v_e=0.2, v_p=0.2, tau_a=0}
% 
%   "DI2CXM"
%   {Fp=0.6, PS=0.2, v_e=0.2, v_p=0.2, fa0.5, tau_a=0, tau_v=0 }
% 
%   "DIIRF"
%   Dual-input, bi-exponential model that fits the functional form of the
%   IRF directly, and can be reduced to any of the (DI)2CXM, GADOXETATE, 
%   MATERNE or TM (but not ETM) models through appropriate fixing of 
%   parameters. See DI-IRF notes on Matlab repository wiki for further
%   explanation, and see functions TWO_CXM_PARAMS_MODEL_TO_PHYS and
%   ACTIVE_PARAMS_MODEL_TO_PHYS for converting the DIRRF outputs into 
%   physiologically meaningful parameters.
%   Outputs 7 parameters: 
%   {Fpos=0.2, Fneg=0.2, Kpos=0.5, Kneg=4.0, fa=0.5, tau_a=0.025, tau_v=0}
%
% See also: RUN_MADYM, TWO_CXM_MODEL, GADOXETATE_MODEL, MATERNE_MODEL,
% EXTENDED_TOFTS_MODEL, TWO_CXM_PARAMS_MODEL_TO_PHYS, ACTIVE_PARAMS_MODEL_TO_PHYS
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
args = u_packargs(varargin, 0, ... 
    'cmd_exe', [local_madym_root 'madym_DCE_lite'],...
    'output_dir', [], ...Output path, will use temp dir if empty;
  	'output_name', 'madym_analysis.dat', ... Name of output file
    'dyn_times', [], ... Time associated with each dynamic signal (in mins), must be supplied if using population AIF
    'input_Ct', 1, ... Flag specifying input dynamic sequence are concentration (not signal) maps
    'output_Ct_sig', nargout > 5,... Flag requesting concentration (derived from signal) are saved to output
    'output_Ct_mod', nargout > 4, ...Flag requesting modelled concentration maps are saved to output
    'no_optimise', 0, ...Flag to switch off optimising, will just fit initial parameters values for model
... The below are all only required if we're converting from signals
    'T1', [], ...Baseline T1 values (in ms)
    'S0', [], ...Baseline S0 values, required if not using ratio method
 	'TR', 0, ... TR of dynamic series (in ms), must be >0 if converting signals
    'FA', NaN, ... Flip angle of dynamic series (degrees), must be set if converting signals
    'r1_const', NaN, ...Relaxivity constant of concentration in tissue (in ms)
    'M0_ratio', true, ... Flag to use ratio method to scale signal instead of supplying S0
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
    'relative_limit_params', [],...Parameters with relative limits on their optimisation bounds
    'relative_limit_values', [], ..._values for relative bounds, sets lower/upper bound as init param -/+ relative limit
...
    'dyn_noise_values', [],...Varying temporal noise in model fit
  	'test_enhancement', false, ...Set test-for-enhancement flag
    'dummy_run', false ...Don't run any thing, just print the cmd we'll run to inspect
    );
clear varargin;

%Get size number of dynamic values - if we've been given a single column
%vector, transpose into a row. Should n_voxels x n_dyns
if size(input_data, 2) == 1
    input_data = input_data';
end
[nSamples, nDyns] = size(input_data);

%If we're converting from signal to concentration, append T1 (and S0 if not
%using ratio method) to pinput data
if ~args.input_Ct
    input_data = [input_data args.T1(:)];
    
    if ~args.M0_ratio
        input_data = [input_data args.S0(:)];
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

%Now set any args that require option inputs
if args.input_Ct
    cmd = sprintf('%s --Ct', cmd);
    
else %Below are only needed if input is signals
    if args.TR
        cmd = sprintf('%s --TR %4.3f', cmd, args.TR);
    end

    if ~isnan(args.FA)
        cmd = sprintf('%s --FA %4.3f', cmd, args.FA);
    end

    if isfinite(args.r1_const)
        cmd = sprintf('%s --r1 %4.3f', cmd, args.r1_const);
    end

    if args.M0_ratio
        cmd = sprintf('%s --M0_ratio', cmd);
    end
end

if isfinite(args.dose)
    cmd = sprintf('%s -D %5.4f', cmd, args.dose);
end

if isfinite(args.hct)
    cmd = sprintf('%s -H %5.4f', cmd, args.hct);
end

if isfinite(args.injection_image)
    cmd = sprintf('%s -i %d', cmd, args.injection_image);
end

if isfinite(args.first_image)
    cmd = sprintf('%s --first %d', cmd, args.first_image);
end

if isfinite(args.last_image)
    cmd = sprintf('%s --last %d', cmd, args.last_image);
end

if args.output_Ct_sig
    cmd = sprintf('%s --Ct_sig', cmd);
end

if args.output_Ct_mod
    cmd = sprintf('%s --Ct_mod', cmd);
end

if args.no_optimise
    cmd = sprintf('%s --no_opt', cmd);
end

if args.test_enhancement
    cmd = sprintf('%s --test_enh', cmd);
end

if ~isempty(args.aif_name)
    cmd = sprintf('%s --aif %s', cmd, args.aif_name);
end

if ~isempty(args.pif_name)
    cmd = sprintf('%s --pif %s', cmd, args.pif_name);
end

if ~isempty(args.dyn_times)
    %Get a name for the temporary file we'll write times to (we'll hold
    %off writing anything until we know this isn't a dummy run
    dyn_times_file = tempname;
    cmd = sprintf('%s -t %s', cmd, dyn_times_file);
end

if ~isempty(args.dyn_noise_values)
    %Get a name for the temporary file we'll write noise to (we'll hold
    %off writing anything until we know this isn't a dummy run
    dynNoise_file = tempname;
    cmd = sprintf('%s --dyn_noise %s', cmd, dynNoise_file);
end

if ~isempty(args.IAUC_times)
    IAUC_str = sprintf('%3.2f', args.IAUC_times(1));
    for i_t = 2:length(args.IAUC_times)
        IAUC_str = sprintf('%s,%3.2f', IAUC_str, args.IAUC_times(i_t));
    end
    cmd = sprintf('%s --iauc %s', cmd, IAUC_str);
end

load_params = false;
if ~isempty(args.init_params)
    if nSamples > 1 && size(args.init_params,1) == nSamples
        input_params_file = tempname;
        load_params = true;
        cmd = sprintf('%s --init_params_file %s', cmd, input_params_file);
    else
        init_str = sprintf('%5.4f', args.init_params(1));
        for i_t = 2:length(args.init_params)
            init_str = sprintf('%s,%5.4f', init_str, args.init_params(i_t));
        end
        cmd = sprintf('%s --init_params %s', cmd, init_str);
    end
end

if ~isempty(args.fixed_params)
    fixed_str = sprintf('%d', args.fixed_params(1));
    for i_t = 2:length(args.fixed_params)
        fixed_str = sprintf('%s,%d', fixed_str, args.fixed_params(i_t));
    end
    cmd = sprintf('%s --fixed_params %s', cmd, fixed_str);
    
    if ~isempty(args.fixed_values)
        fixed_str = sprintf('%5.4f', args.fixed_values(1));
        for i_t = 2:length(args.fixed_values)
            fixed_str = sprintf('%s,%5.4f', fixed_str, args.fixed_values(i_t));
        end
        cmd = sprintf('%s --fixed_values %s', cmd, fixed_str);
    end
end

if ~isempty(args.relative_limit_params)
    relative_str = sprintf('%d', args.relative_limit_params(1));
    for i_t = 2:length(args.relative_limit_params)
        relative_str = sprintf('%s,%d', relative_str, args.relative_limit_params(i_t));
    end
    cmd = sprintf('%s --relative_limit_params %s', cmd, relative_str);
    
    if ~isempty(args.relative_limit_values)
        relative_str = sprintf('%5.4f', args.relative_limit_values(1));
        for i_t = 2:length(args.relative_limit_values)
            relative_str = sprintf('%s,%5.4f', relative_str, args.relative_limit_values(i_t));
        end
        cmd = sprintf('%s --relative_limit_values %s', cmd, relative_str);
    end
end

if args.dummy_run
    %Don't actually run anything, just print the command
    fprintf('%s\n', cmd);
    
    model_params = [];
    model_fit = [];
    iauc = [];
    error_codes = [];
    model_conc = [];
    dyn_conc = [];
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

n_iauc = length(args.IAUC_times);
iauc = outputData(:, 3 + (1:n_iauc));

%Workout which columns hold parameter values
n_cols = size(outputData, 2);

param_col1 = n_iauc + 4;
param_col2 = n_cols;
if args.output_Ct_mod
    if args.output_Ct_sig
        %Have another nDyns of dynamic concentration values
        param_col2 = param_col2 - nDyns;
        dyn_conc = outputData(:, param_col2+(1:nDyns));
    end
    
    %After the params we have nDyns of model concentration values
    param_col2 = param_col2 - nDyns;
    model_conc = outputData(:, param_col2+(1:nDyns));
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




