function [status, result] = run_madym_DCE(varargin)
%RUN_MADYM_DCE wrapper function to call C++ tool Madym. Fits
%   tracer-kinetic models to DCE time-series stored in Analyze format images,
%   saving the model parameters and modelled concentration time-series also
%   in Analyze format. The status and result of the system call to Madym is
%   returned.
%   [status, result] = ...
%       run_madym_DCE(varargin)
%
% RUN_MADYM_DCE uses the U_PACKARGS interface function
% and so arguments to the function may be passed as name-value pairs
% or as a struct with fields with corresponding names. See below for
%       a full description of all options.
%
% Outputs:
%      [status,result] returned by the system call to the Madym executable.
%      These may be operating system dependent, however status=0 should
%      mean an error free operation. If status is non-zero an error has
%      probably occurred, the result of which may be set in result. In any
%      event, it is best to check the output_dir, and any program logs that
%      have been saved there.
%
% Examples:
%   Fitting to concentration time-series. If using a population AIF, you
%   must supply a vector of dynamic times. A population AIF (PIF) is used
%   if the aifName (pifName) option is left empty.
%   [status,result] = ...
%       run_madym_DCE('model', "2CXM", 'output_dir', 'C:\QBI\mdm_analysis\')
%
% Notes:
%   Tracer-kinetic models:
%
%   All models available in the main MaDym and MaDym-Lite C++ tools are
%   available to fit.
% 
%   See the madym_cxx project wiki for more details:
%   https://gitlab.com/manchester_qbi/manchester_qbi_public/madym_cxx/-/wikis/dce_models
% 
%   Run: system([local_madym_root 'madym_DCE --help']); to see full set of
%   input options to C++ tool
% 
%
% See also: RUN_MADYM_LITE, RUN_MADYM_T1
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
    'cmd_exe', [local_madym_root 'madym_DCE'],...
    'config', '',... Path to a config file to set default options
    'model', '',... Model to fit, specified by its name in CAPITALS, see notes for options
    'output_dir', '',...Folder in which output maps will be saved 
    'no_optimise', false, ...Flag to switch off optimising, will just fit initial parameters values for model
    'T1_vols', [], ..._file names of signal volumes to from which baseline T1 is mapped
    'T1_method', '',...Method for mapping baseline T1 (eg VFA)
    'dynamic_basename', [], ...Template name for dynamic sequences eg. dynamic/dyn_
    'dyn_dir', '',...Dir to dynamic series (can be included in dynamic basename)
    'sequence_format', '',...Format for converting dynamic series index to string, eg %01u
    'sequence_start', NaN,...Start index of dynamic series
    'sequence_step', NaN,...Step size between indexes in dynamic series
    'n_dyns', 0,...Number of dynamic sequence maps to load. If <=0, loads all maps in dynamic dir matching -dyn pattern
    'input_Ct', false,...Flag specifying input dynamic sequence are concentration (not signal) maps
    'output_Ct_sig', false,...Flag requesting concentration (derived from signal) are saved to output
    'output_Ct_mod', false,...Flag requesting modelled concentration maps are saved to output
    'T1_name', '',...Path to T1 map
    'M0_name', '',...Path to M0 map
    'B1_name', '',...Path to B1 correction map
    'B1_correction', false, ... Apply B1 correction
    'B1_scaling', NaN, ... Scaling factor to use with B1 map
    'r1_const', NaN,...Relaxivity constant of concentration in tissue (in ms)
    'TR', NaN, ... TR of dynamic sequence, only required if T1 method is IR
    'M0_ratio', false,...Flag to use ratio method to scale signal instead of supplying M0
    'dose', NaN,...Concentration dose (mmole/kg)
    'injection_image', NaN,...Injection image
    'hct', NaN,...Haematocrit correction
    'T1_noise', NaN,...PD noise threshold
    'first_image', NaN,...First image used to compute model fit
    'last_image', NaN,...Last image used to compute model fit
    'roi_name', '',...Path to ROI map
    'aif_name', '',...Path to precomputed AIF if not using population AIF
    'aif_map', '',...Path to to mask from which AIF will be computed on the fly
    'pif_name', '',...Path to precomputed PIF if not deriving from AIF
    'IAUC_times', [],..._times (in s) at which to compute IAUC values
    'IAUC_at_peak', false,...Flag requesting IAUC computed at peak signal
    'param_names', [],...Names of model parameters to be optimised, used to name the output parameter maps
    'init_params', [],...Initial values for model parameters to be optimised, either as single vector, or 2D array NSamples x N_params
    'fixed_params', [],...Parameters fixed to their initial values (ie not optimised)
    'fixed_values', [],..._values for fixed parameters (overrides default initial parameter values)
    'upper_bounds', [],...Upper bounds for each parameter during optimisation
    'lower_bounds', [],...Lower bounds for each parameter during optimisation
    'relative_limit_params', [],...Parameters with relative limits on their optimisation bounds
    'relative_limit_values', [],..._values for relative bounds, sets lower/upper bound as init param -/+ relative limit
    'residuals', '',... Path to existing residuals map
    'init_maps_dir', '',...Path to directory containing maps of parameters to initialise fit (overrides init_params)
    'init_map_params', [],...Parameters initialised from maps (if empty and init_maps_dir set, all params from maps)
    'dyn_noise', false,...Set to use varying temporal noise in model fit
    'test_enhancement', false,...Set test-for-enhancement flag
    'max_iter', NaN,... Maximum number of iterations in model fit
    'img_fmt_r', '',...Set image read format
    'img_fmt_w', '',...Set image write format
    'overwrite', false,...Set overwrite existing analysis in output dir
    'no_audit', false,... Turn off audit log
    'no_log', false,... Turn off propgram log
    'quiet', false,... Suppress output to stdout
    'program_log_name', '',...Program log file name
    'audit_name', '',...Audit file name
    'audit_dir', '',...Folder in which audit logs are saved
    'config_out', '',...Name of output config file
    'error_name', '',...Error codes image file name
    'working_directory', '',...Sets the current working directory for the system call, allows setting relative input paths for data
    'dummy_run', false);...Don't run any thing, just print the cmd we'll run to inspect
clear varargin;

%Set up base command
cmd = sprintf('%s ', args.cmd_exe);
    
%Check if a config file exists
if ~isempty(args.config)
    cmd = sprintf('%s --config %s', cmd, args.config);   
end

%Set the working dir
if ~isempty(args.working_directory)
    cmd = sprintf('%s --cwd %s', cmd, args.working_directory);
end

%Set TK model
if ~isempty(args.model)
    cmd = sprintf('%s -m %s', cmd, args.model);
end

%Set output directory
if ~isempty(args.output_dir)
    cmd = sprintf('%s -o %s', cmd, args.output_dir);
end

%Set the dynamic names
if ~isempty(args.dynamic_basename)
    cmd = sprintf('%s -d %s', cmd, args.dynamic_basename);
end

if ~isempty(args.dyn_dir)
    cmd = sprintf('%s --dyn_dir %s', cmd, args.dyn_dir);
end

%Set the format
if ~isempty(args.sequence_format)
    cmd = sprintf('%s --sequence_format %s', cmd, args.sequence_format);
end

if isfinite(args.sequence_start)
    cmd = sprintf('%s --sequence_start %d', cmd, args.sequence_start);
end

if isfinite(args.sequence_step)
    cmd = sprintf('%s --sequence_step %d', cmd, args.sequence_step);
end

%Set the number of dynamics
if args.n_dyns
    cmd = sprintf('%s --n_dyns %d', cmd, args.n_dyns);
end

%Set image formats
if ~isempty(args.img_fmt_r)
    cmd = sprintf('%s --img_fmt_r %s', cmd, args.img_fmt_r);
end
if ~isempty(args.img_fmt_w)
    cmd = sprintf('%s --img_fmt_w %s', cmd, args.img_fmt_w);
end

%Now set any args that require option inputs
if args.input_Ct
    cmd = sprintf('%s --Ct', cmd);
end

%Set T1 options
if ~isempty(args.T1_name)
    cmd = sprintf('%s --T1 %s', cmd, args.T1_name);
end

if ~isempty(args.M0_name)
    cmd = sprintf('%s --M0 %s', cmd, args.M0_name);
end

if ~isempty(args.T1_method)
    cmd = sprintf('%s --T1_method %s', cmd, args.T1_method);
end 

if ~isempty(args.T1_vols)
    %Set VFA files in the options string
    t1_str = sprintf('%s', args.T1_vols{1});
    for i_t = 2:length(args.T1_vols)
        t1_str = sprintf('%s,%s', t1_str, args.T1_vols{i_t});
    end
    cmd = sprintf('%s --T1_vols %s', cmd, t1_str);
end

if isfinite(args.T1_noise)
    cmd = sprintf('%s --T1_noise %5.4f', cmd, args.T1_noise);
end


%Set any other options required to convert signal to concentration
if isfinite(args.r1_const)
    cmd = sprintf('%s --r1 %4.3f', cmd, args.r1_const);
end

if isfinite(args.TR)
    cmd = sprintf('%s --TR %4.3f', cmd, args.TR);
end

if isfinite(args.dose)
    cmd = sprintf('%s -D %4.3f', cmd, args.dose);
end

if args.M0_ratio
    cmd = sprintf('%s --M0_ratio', cmd);
end

%B1 correction options
if ~isempty(args.B1_name)
    cmd = sprintf('%s --B1 %s', cmd, args.B1_name);
end    
if args.B1_correction
    cmd = sprintf('%s --B1_correction', cmd);
end
if isfinite(args.B1_scaling)
    cmd = sprintf('%s --B1_scaling %d', cmd, args.B1_scaling);
end    

%Now go through all the other optional parameters, and if they've been set,
%set the necessary option flag in the cmd string
if args.no_optimise
    cmd = sprintf('%s --no_opt', cmd);
end

if args.output_Ct_sig
    cmd = sprintf('%s --Ct_sig', cmd);
end

if args.output_Ct_mod
    cmd = sprintf('%s --Ct_mod', cmd);
end

if args.test_enhancement
    cmd = sprintf('%s --test_enh', cmd);
end

if isfinite(args.max_iter)
    cmd = sprintf('%s --max_iter %d', cmd, args.max_iter);
end

if args.dyn_noise
    cmd = sprintf('%s --dyn_noise', cmd);
end

if args.overwrite
    cmd = sprintf('%s --overwrite', cmd);
end

if args.no_audit
    cmd = sprintf('%s --no_audit', cmd);
end

if args.no_log
    cmd = sprintf('%s --no_log', cmd);
end

if args.quiet
    cmd = sprintf('%s --quiet', cmd);
end

if isfinite(args.hct)
    cmd = sprintf('%s -H %d', cmd, args.hct);
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

if ~isempty(args.roi_name)
    cmd = sprintf('%s --roi %s', cmd, args.roi_name);
end

if ~isempty(args.residuals)
    cmd = sprintf('%s --residuals %s', cmd, args.residuals);
end

if ~isempty(args.aif_name)
    cmd = sprintf('%s --aif %s', cmd, args.aif_name);
end

if ~isempty(args.aif_map)
    cmd = sprintf('%s --aif_map %s', cmd, args.aif_map);
end

if ~isempty(args.pif_name)
    cmd = sprintf('%s --pif %s', cmd, args.pif_name);
end

if ~isempty(args.program_log_name)
    cmd = sprintf('%s --program_log %s', cmd, args.program_log_name);
end

if ~isempty(args.audit_name)
    cmd = sprintf('%s --audit %s', cmd, args.audit_name);
end

if ~isempty(args.audit_dir)
    cmd = sprintf('%s --audit_dir %s', cmd, args.audit_dir);
end

if ~isempty(args.config_out)
    cmd = sprintf('%s --config_out %s', cmd, args.config_out);
end

if ~isempty(args.error_name)
    cmd = sprintf('%s -E %s', cmd, args.error_name);
end

if ~isempty(args.IAUC_times)
    IAUC_str = sprintf('%3.2f', args.IAUC_times(1));
    for i_t = 2:length(args.IAUC_times)
        IAUC_str = sprintf('%s,%3.2f', IAUC_str, args.IAUC_times(i_t));
    end
    cmd = sprintf('%s --iauc %s', cmd, IAUC_str);
end

if args.IAUC_at_peak
    cmd = sprintf('%s --iauc_peak %s', cmd);
end

if ~isempty(args.init_maps_dir)
    cmd = sprintf('%s --init_maps %s', cmd, init_maps_dir);
end
if ~isempty(args.init_params)
    init_str = sprintf('%d', args.init_params(1));
    for i_t = 2:length(args.init_params)
        init_str = sprintf('%s,%d', init_str, args.init_params(i_t));
    end
    cmd = sprintf('%s --init_params %s', cmd, init_str);
end

if ~isempty(args.param_names)
    param_str = sprintf('%s', args.param_names{1});
    for i_t = 2:length(args.param_names)
        param_str = sprintf('%s,%d', param_str, args.param_names{i_t});
    end
    cmd = sprintf('%s --param_names %s', cmd, param_str);
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

if ~isempty(args.upper_bounds)
    upper_str = sprintf('%d', args.upper_bounds(1));
    for i_t = 2:length(args.upper_bounds)
        upper_str = sprintf('%s,%d', upper_str, args.upper_bounds(i_t));
    end
    cmd = sprintf('%s --upper_bounds %s', cmd, upper_str);
end

if ~isempty(args.lower_bounds)
    lower_str = sprintf('%d', args.lower_bounds(1));
    for i_t = 2:length(args.lower_bounds)
        lower_str = sprintf('%s,%d', lower_str, args.lower_bounds(i_t));
    end
    cmd = sprintf('%s --lower_bounds %s', cmd, lower_str);
end

if ~isempty(args.relative_limit_params)
    relative_str = sprintf('%d', args.relative_limit_params(1));
    for i_t = 2:length(args.relative_limit_params)
        relative_str = sprintf('%s,%d', ...
            relative_str, args.relative_limit_params(i_t));
    end
    cmd = sprintf('%s --relative_limit_params %s', cmd, relative_str);
    
    if ~isempty(args.relative_limit_values)
        relative_str = sprintf('%5.4f', args.relative_limit_values(1));
        for i_t = 2:length(args.relative_limit_values)
            relative_str = sprintf('%s,%5.4f',...
                relative_str, args.relative_limit_values(i_t));
        end
        cmd = sprintf('%s --relative_limit_values %s', cmd, relative_str);
    end
end

if args.dummy_run
    %Don't actually run anything, just print the command
    status = []; result = [];
    fprintf('%s\n', cmd);  
    return;  
end

%Otherwise we can run the command:
[status, result] = system(cmd, '-echo');

%--------------------------------------------------------------------------
function run_test()

%Generate an concentration time-series using the ETM
nDyns = 100;
ktrans = 0.25;
ve = 0.2;
vp = 0.1;
tau = 0;
injection_img = 8;
t = linspace(0, 5, nDyns);
Ca_t = population_aif(t, injection_img);
C_t = extended_tofts_model(ktrans, ve, vp, tau, Ca_t, t);

%Add some noise and rescale so baseline mean is 0
C_tn = C_t + randn(1,100)/100;
C_tn = C_tn - mean(C_tn(1:injection_img),2);

%Convert the concentrations to signals with some notional T1 values and
%refit using signals as input
FA = 20;
TR = 3.5;
T1_0 = 1000;
r1_const = 3.4;
S_t0 = 100;
S_tn = concentration_to_signal(...
    C_tn, FA, TR, T1_0, S_t0, r1_const, injection_img);

dyn_dir = [tempdir 'dynamics/'];
create_folder(dyn_dir);

St_names = cell(nDyns,1);
Ct_names = cell(nDyns,1);
for i_dyn = 1:nDyns
    
    %Write out 1x1 concentration maps and xtr files
    Ct_name = [dyn_dir 'Ct_' num2str(i_dyn)];
    Ct_names{i_dyn} =  [Ct_name '.hdr'];
    Ct_xtr_name = [Ct_name '.xtr'];
    
    %Convert dynamic time (in minutes) to xtr timestamp
    t_in_secs = 60*t(i_dyn);
    hh = floor(t_in_secs / (3600));
    mm = floor((t_in_secs - 3600*hh) / 60);
    ss = t_in_secs - 3600*hh - 60*mm;
    timestamp = 10000*hh + 100*mm + ss;

    save_img_volume(C_tn(i_dyn), Ct_names{i_dyn}, [1 1 1], [], [], 0);
    write_xtr_file(Ct_xtr_name, 0, ...
        'FlipAngle', FA,...
        'TR', TR,...
        'TimeStamp', timestamp);
    
    %Write out 1x1 signal maps and xtr files
    St_name = [dyn_dir 'St_' num2str(i_dyn)];
    St_names{i_dyn} =  [St_name '.hdr'];
    St_xtr_name = [St_name '.xtr'];

    save_img_volume(S_tn(i_dyn), St_names{i_dyn}, [1 1 1], [], [], 0);
    write_xtr_file(St_xtr_name, 0, ...
        'FlipAngle', FA,...
        'TR', TR,...
        'TimeStamp', timestamp);  
    
end
T1_name = [tempdir 'T1.hdr'];
save_img_volume(T1_0, T1_name, [1 1 1], [], [], 0);
%%
%Apply Madym to concentration maps
Ct_output_dir = [tempdir 'mdm_analysis_Ct/'];
[status, result] = run_madym_DCE(...
    'model', 'ETM', ...
    'output_dir', Ct_output_dir,...
    'dynamic_basename', [dyn_dir 'Ct_'],...
    'input_Ct', 1,...
    'output_Ct_sig',1,...
    'output_Ct_mod',1,...
    'injection_image', injection_img,...
    'img_fmt_r', 'ANALYZE',...
    'img_fmt_w', 'ANALYZE',...
    'overwrite', 1,...
    'no_audit', 1);

fprintf('status: %d\nresult: %s\n', status, result);

ktrans_fit = load_img_volume([Ct_output_dir 'Ktrans.hdr']);
ve_fit = load_img_volume([Ct_output_dir 'v_e.hdr']);
vp_fit = load_img_volume([Ct_output_dir 'v_p.hdr']);
tau_fit = load_img_volume([Ct_output_dir 'tau_a.hdr']);

fprintf('Parameter estimation (actual, fitted concentration, fitted signal)\n');
fprintf('Ktrans: (%3.2f, %3.2f)\n', ktrans, ktrans_fit);
fprintf('Ve: (%3.2f, %3.2f)\n', ve, ve_fit);
fprintf('Vp: (%3.2f, %3.2f)\n', vp, vp_fit);
fprintf('Tau: (%3.2f, %3.2f)\n', tau, tau_fit);

Cs_t = zeros(1,100);
Cm_t = zeros(1,100);
for i_dyn = 1:nDyns
    Cs_t(i_dyn) = ...
        load_img_volume([Ct_output_dir 'Ct_sig/Ct_sig' num2str(i_dyn) '.hdr']);
    Cm_t(i_dyn) = ...
        load_img_volume([Ct_output_dir 'Ct_mod/Ct_mod' num2str(i_dyn) '.hdr']);
end

figure('Name', 'madym test applied');
subplot(1,2,1);
plot(t, C_tn); hold all;
plot(t, Cs_t, '--');
plot(t, Cm_t);
legend({'C(t)', 'C(t) (output from MaDym)', 'ETM model fit'});
xlabel('Time (mins)');
ylabel('Voxel concentration');
title(sprintf('Input C(t): Model fit SSE = %4.3f', sum((C_tn-Cm_t).^2)));
%%
%Apply Madym to signal maps
St_output_dir = [tempdir 'mdm_analysis_St/'];
[status, result] = run_madym_DCE(...
    'model', 'ETM', ...
    'output_dir', St_output_dir,...
    'dynamic_basename', [dyn_dir 'St_'],...
    'input_Ct', 0,...
    'output_Ct_sig',1,...
    'output_Ct_mod',1,...
    'T1_name', T1_name,...
    'r1_const', r1_const,...
    'injection_image', injection_img,...
    'img_fmt_r', 'ANALYZE',...
    'img_fmt_w', 'ANALYZE',...
    'overwrite', 1,...
    'no_audit', 1);

fprintf('status: %d\nresult: %s\n', status, result);

ktrans_fit = load_img_volume([St_output_dir 'Ktrans.hdr']);
ve_fit = load_img_volume([St_output_dir 'v_e.hdr']);
vp_fit = load_img_volume([St_output_dir 'v_p.hdr']);
tau_fit = load_img_volume([St_output_dir 'tau_a.hdr']);

fprintf('Parameter estimation (actual, fitted concentration, fitted signal)\n');
fprintf('Ktrans: (%3.2f, %3.2f)\n', ktrans, ktrans_fit);
fprintf('Ve: (%3.2f, %3.2f)\n', ve, ve_fit);
fprintf('Vp: (%3.2f, %3.2f)\n', vp, vp_fit);
fprintf('Tau: (%3.2f, %3.2f)\n', tau, tau_fit);

Cs_t = zeros(1,100);
Cm_t = zeros(1,100);
for i_dyn = 1:nDyns
    Cs_t(i_dyn) = ...
        load_img_volume([St_output_dir 'Ct_sig/Ct_sig' num2str(i_dyn) '.hdr']);
    Cm_t(i_dyn) = ...
        load_img_volume([St_output_dir 'Ct_mod/Ct_mod' num2str(i_dyn) '.hdr']);
end

subplot(1,2,2);
plot(t, C_tn); hold all;
plot(t, Cs_t, '--');
plot(t, Cm_t);
legend({'C(t)', 'C(t) (output from MaDym)', 'ETM model fit'});
xlabel('Time (mins)');
ylabel('Voxel concentration');
title(sprintf('Input S(t): Model fit SSE = %4.3f', sum((C_tn-Cm_t).^2)));
%%
%Delete input and output files
for i_dyn = 1:nDyns
    delete(Ct_names{i_dyn});
    delete([Ct_names{i_dyn}(1:end-4) '.img']);
    delete([Ct_names{i_dyn}(1:end-4) '.xtr']);
    
    delete(St_names{i_dyn});
    delete([St_names{i_dyn}(1:end-4) '.img']);
    delete([St_names{i_dyn}(1:end-4) '.xtr']);
end  
delete(T1_name);

delete_dirs = {...
    [Ct_output_dir 'Ct_sig'],...
    [Ct_output_dir 'Ct_mod'],...
    Ct_output_dir,...
    [St_output_dir 'Ct_sig'],...
    [St_output_dir 'Ct_mod'],...
    St_output_dir};

for i_d = 1:length(delete_dirs)
    del_dir = delete_dirs{i_d};

    output_files = dir(del_dir);
    for i_f = 3:length(output_files)
        delete(fullfile(del_dir, output_files(i_f).name));
    end
    rmdir(del_dir);
end

%% ------------------------------------------------------------------------