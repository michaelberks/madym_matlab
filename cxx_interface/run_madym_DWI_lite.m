function [model_params, model_fit, error_codes] =...
    run_madym_DWI_lite(model, signals, B_values, varargin)
%RUN_MADYM_DWI_LITE wrapper function to call C++ tool madym_DWI_lite. Fits
%   diffusion models to DWU signals inputs, returning the fitted parameters
%   and model sum-of-square residuals
%   [model_params, model_fit, error_codes] = ...
%       run_madym_DWI_lite(model, input_data, varargin)
%
% RUN_MADYM_DWI_LITE uses the U_PACKARGS interface function
% and so arguments to the function may be passed as name-value pairs
% or as a struct with fields with corresponding names. The field names
% are defined as:
%
% Mandatory Arguments:
%       model (str) - Model to fit, specified by its name in CAPITALS,
%            see notes for options
%
%       input_data 2D array (Nsamples x NBvalues) - Array of signal data, one
%       sample per row.
%
%       B_values (NBvalues or Nsamples x NBvalues) - Array of B-values.
%       Either a single 1D array with length matching the number of columns
%       in signals, in which case the same B-values are used for all
%       samples. Or an array the same size as signals, in which case
%       separate B-values can be specified for each sample.
%
% Optional Arguments: This wrapper uses the upackargs interface to allow
%       optional arguments to be entered as name/value pairs. See below for
%       a full description of all options.
%
% Outputs:
%      model_params (2D array, Nsamples x Nparams) - each row contains 
%       the estimated model parameters for the corresponding signal-series 
%       in the input_data. The number of columns depends on the model 
%       parameters. The parameters are returned in the same order as the 
%       defined in Madym. See notes for details for each model.
%
%      model_fit (1D array, Nsamples x 1) - sum-of-squared model 
%       residuals for each signal-series
%
%      error_codes (2D array, Nsamples x 2) - error codes returned by MaDym
%       for fitting each sample. 0 implies no errors or warnings. For all
%       non-zero values refer to Madym documentation for details
%
% Examples:
%
% Notes:
%   Diffusion models:
% 
%   All models available in the main madym_DWI C++ tools are
%   available to fit.
%
%   See the madym_cxx project wiki for more details:
% https://gitlab.com/manchester_qbi/manchester_qbi_public/madym_cxx/-/wikis/dwi_models
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
    'cmd_exe', [local_madym_root 'madym_DWI_lite'],...
    'Bvals_thresh', [], ...Thresholds used in IVIM fitting
    'output_dir', [], ...Output path, will use temp dir if empty;
  	'output_name', 'madym_analysis.dat', ... Name of output file
    'quiet', NaN,... Suppress output to stdout
    'dummy_run', false ...Don't run any thing, just print the cmd we'll run to inspect
    );
clear varargin;

%Fit directly supplied FA and signal data using calculate_T1_lite
[nSamples, nSignals] = size(signals);

%Do error checking on required inputs 
if numel(B_values) == nSignals
    B_values = repmat(B_values(:)', nSamples, 1);
elseif ~all(size(B_values)==[nSamples, nSignals])
    error('Size of B values array does not match size of signals array');
end

%Set up temporary files for B-values and signals (we'll hold
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

cmd = sprintf(...
    '%s --DWI_model %s --data %s --n_DWI %d -o %s -O %s',...
    args.cmd_exe,...
    model, ...
    input_file,...
    nSignals,...
    args.output_dir,...
    args.output_name);

cmd = add_option('float_list', cmd, '--Bvals_thresh', args.Bvals_thresh);

%Check for bad samples, these can screw up Madym as the lite version
%of the software doesn't do the full range of error checks Madym proper
%4.3foes. So chuck them out now and warn the user
discard_samples = ...
    any( isnan(B_values) |...
    ~isfinite(B_values), 2) |...
    any(isnan(signals) |...
    ~isfinite(signals), 2);

if any(discard_samples)
    warning(['Samples with NaN values found,'...
        'these will be set to zero for model-fitting']);
    B_values(discard_samples,:) = 0;
    signals(discard_samples,:) = 0;
end    

if args.dummy_run
    %Don't actually run anything, just print the command
    fprintf('%s\n', cmd);
    
    model_params = [];
    model_fit = [];
    error_codes = [];
    return;
    
end

%Write input values to a temporary file
fid = fopen(input_file, 'wt');
for i_row = 1:nSamples
    fprintf(fid, '%6.5f ', B_values(i_row,:));
    fprintf(fid, '%6.5f ', signals(i_row,:));
    fprintf(fid, '\n');
end
fclose(fid);

%At last.. we can run the command
system(cmd, '-echo');

%Now load the output from madym and extract data to match this
%functions outputs
outputData = load(fullOutPath);
model_params = outputData(:,1:end-2);
model_fit = outputData(:,end-1);
error_codes = outputData(:,end);

%Tidy up temporary files
delete(input_file);   
if deleteOutput
    delete(fullOutPath);
end


%--------------------------------------------------------------------------

%%
%Test function to run if no inputs
function run_test()

%Run tests for ADC and IVIM models
sigma = 1;
B_vals = [0.0, 20.0, 40.0, 60.0, 80.0, 100.0, 300.0, 500.0, 800.0];

%Generate ADC test data with Rician noise added
S0 = 100;
ADC = 0.8e-3;
signals = ADC_model(B_vals, S0, ADC);
signals_n = add_rician_noise(signals, sigma);

%Use madym lite to fit this data
[model_params, model_fit] = run_madym_DWI_lite(...
    'ADC', signals_n, B_vals);
S0_f = model_params(1);
ADC_f = model_params(2);
signals_f = ADC_model(B_vals, S0_f, ADC_f );

%Display plots of the fit
figure('Name', 'madym_DWI_lite test applied');
subplot(1,2,1);

plot(B_vals, signals, 'k--', 'linewidth', 2); hold all;
plot(B_vals, signals_n, 'rx', 'markersize', 10);
plot(B_vals, signals_f, 'b-', 'linewidth', 2);
plot(B_vals, signals_f, 'go', 'markersize', 10);
legend({'Noise free signals' 
    'Signals + Rician noise' 
    'Fitted ADC model signals'
    ''});
xlabel('B-values (msecs)');
ylabel('Signal');
title({'ADC: Parameter estimates (actual,fit)';...
    sprintf('S0: (%4.3f, %4.3f), ADC: (%4.3f, %4.3f)',...
    S0, S0_f, 1e3*ADC, 1e3*ADC_f)});

%Generate IVIM test data with Rician noise added
S0 = 100;
D = 0.8e-3;
f = 0.2;
D_star = 15e-3;
Bvals_thresh = [40.0,60.0,100.0,150.0];

signals = IVIM_model(B_vals, S0, D, f, D_star);
signals_n = add_rician_noise(signals, sigma);

%Use madym lite to fit this data
[model_params, model_fit] = run_madym_DWI_lite(...
    'IVIM', signals_n, B_vals,...
    'Bvals_thresh', Bvals_thresh);
S0_f = model_params(1);
D_f = model_params(2);
f_f = model_params(3);
D_star_f = model_params(4);
signals_f = IVIM_model(B_vals, S0_f, D_f, f_f, D_star_f);

%Display plots of the fit
subplot(1,2,2);
plot(B_vals, signals, 'k--', 'linewidth', 2); hold all;
plot(B_vals, signals_n, 'rx', 'markersize', 10);
plot(B_vals, signals_f, 'b-', 'linewidth', 2);
plot(B_vals, signals_f, 'go', 'markersize', 10);
legend({'Noise free signals' 
    'Signals + Rician noise' 
    'Fitted IVIM model signals'
    ''});
xlabel('B-values (msecs)');
ylabel('Signal');
title({'IVIM: Parameter estimates (actual,fit)';...
    sprintf('S0: (%4.3f, %4.3f), D: (%4.3f, %4.3f), f: (%4.3f, %4.3f), D*: (%4.3f, %4.3f)',...
    S0, S0_f, 1e3*D, 1e3*D_f, f, f_f, 1e3*D_star, 1e3*D_star_f)});




