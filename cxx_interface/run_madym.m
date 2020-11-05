function [status, result] = run_madym(model, output_dir, varargin)
%RUN_MADYM wrapper function to call C++ tool Madym. Fits
%   tracer-kinetic models to DCE time-series stored in Analyze format images,
%   saving the model parameters and modelled concentration time-series also
%   in Analyze format. The status and result of the system call to Madym is
%   returned.
%   [model_params, model_fit, error_codes, model_conc, dyn_conc] = ...
%       run_madym(model, output_dir, varargin)
%
% RUN_MADYM uses the U_PACKARGS interface function
% and so arguments to the function may be passed as name-value pairs
% or as a struct with fields with corresponding names. The field names
% are defined as:
%
% Mandatory Arguments:
%       model (str) - Model to fit, specified by its name in CAPITALS,
%            see notes for options
%
%       output_dir (str) - Folder in which output maps will be saved 
%       If output_dir doesn't exist, it will be created. If it exists and
%       is not empty, will generate error unless 'overwrite' flag set
%
% Optional Arguments: This wrapper uses the upackargs interface to allow
%       optional arguments to be entered as name/value pairs. See below for
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
%       run_madym("2CXM", 'C:\QBI\mdm_analysis\')
%
%   Fitting to concentration time-series using a patient specific AIF. The
%   AIF should be defined in a text file with two columns containing the
%   dynamic times and associated AIF value at each times respectively. Pass
%   the full filepath as input
%   [model_params, model_fit, error_codes, model_conc] = ...
%       run_madym_lite("2CXM", dyn_conc, 'aifName', 'C:\DCE_data\pt_AIF.txt')
%
%   Fitting to signals - Set inputCt to false and use options to supply
%   T1 values (and TR, FA, relax_coeff etc) to convert signals to
%   concentration.
%   [model_params, model_fit, error_codes, model_conc, dyn_conc] = ...
%       run_madym_lite("2CXM", dyn_signals, 'dyn_times', t,...
%           'inputCt', 0, 'T1', T1_vals, 'TR', TR, 'FA', FA)
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
%   To run a standard Tofts model, use 
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
% See also: RUN_MADYM_LITE, TWO_CXM_MODEL, GADOXETATE_MODEL, MATERNE_MODEL,
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
    'cmd_exe', [local_madym_root 'madym_DCE'],...
    'no_optimise', false, ...Flag to switch off optimising, will just fit initial parameters values for model
    'T1_vols', [], ..._file names of signal volumes to from which baseline T1 is mapped
    'dynamic_basename', [], ...Template name for dynamic sequences eg. dynamic/dyn_
    'dynamic_name_format', '',...Format for converting dynamic series index to string, eg %01u
    'n_dyns', 0,...Number of dynamic sequence maps to load. If <=0, loads all maps in dynamic dir matching -dyn pattern
    'input_Ct', true,...Flag specifying input dynamic sequence are concentration (not signal) maps
    'output_Ct_sig', true,...Flag requesting concentration (derived from signal) are saved to output
    'output_Ct_mod', true,...Flag requesting modelled concentration maps are saved to output
    'T1_name', '',...Path to T1 map
    'M0_name', '',...Path to M0 map
    'r1_const', NaN,...Relaxivity constant of concentration in tissue (in ms)
    'M0_ratio', true,...Flag to use ratio method to scale signal instead of supplying M0
    'dose', NaN,...Concentration dose (mmole/kg)
    'injection_image', NaN,...Injection image
    'hct', NaN,...Haematocrit correction
    'T1_noise', NaN,...PD noise threshold
    'first_image', NaN,...First image used to compute model fit
    'last_image', NaN,...Last image used to compute model fit
    'roi_name', '',...Path to ROI map
    'aif_name', '',...Path to precomputed AIF if not using population AIF
    'pif_name', '',...Path to precomputed PIF if not deriving from AIF
    'IAUC_times', [],..._times (in s) at which to compute IAUC values
    'param_names', [],...Names of model parameters to be optimised, used to name the output parameter maps
    'init_params', [],...Initial values for model parameters to be optimised, either as single vector, or 2D array NSamples x N_params
    'fixed_params', [],...Parameters fixed to their initial values (ie not optimised)
    'fixed_values', [],..._values for fixed parameters (overrides default initial parameter values)
    'relative_limit_params', [],...Parameters with relative limits on their optimisation bounds
    'relative_limit_values', [],..._values for relative bounds, sets lower/upper bound as init param -/+ relative limit
    'init_maps_dir', '',...Path to directory containing maps of parameters to initialise fit (overrides init_params)
    'init_map_params', [],...Parameters initialised from maps (if empty and init_maps_dir set, all params from maps)
    'dyn_noise', false,...Set to use varying temporal noise in model fit
    'test_enhancement', false,...Set test-for-enhancement flag
    'sparse_write', false,...Set sparseWrite to save output map sin sparse mode
    'overwrite', false,...Set overwrite existing analysis in output dir
    'program_log_name', '',...Program log file name
    'audit_name', '',...Audit file name
    'error_name', '',...Error codes image file name
    'working_directory', '',...Sets the current working directory for the system call, allows setting relative input paths for data
    'dummy_run', false);...Don't run any thing, just print the cmd we'll run to inspect
clear varargin;

cmd = sprintf('%s -m %s -o %s ',...
    args.cmd_exe, model, output_dir);  

%Set the working dir
if ~isempty(args.working_directory)
    cmd = sprintf('%s --cwd %s', cmd, args.working_directory);
end

%Set the dynamic names
if ~isempty(args.dynamic_basename)
    cmd = sprintf('%s -d %s', cmd, args.dynamic_basename);
end

%Now set any args that require option inputs
if args.input_Ct
    cmd = sprintf('%s --Ct', cmd);
    
else %Below are only needed if input is signals
    
    %Check if we have VFA files
    if isempty(args.T1_vols)
        
        %If not we need a T1 map
        if isempty(args.T1_name)
            error(['You must supply either T1 input files ' ... 
                'to compute baseline T1 or a T1 map']);
        else
            cmd = sprintf('%s --T1 %s', cmd, args.T1_name);
        end
        
        %And if we're not using the ratio method, an M0 map
        if ~args.M0_ratio
            if isempty(args.M0_name)
                error('If M0_name is off, you must supply an M0 map');
            else
                cmd = sprintf('%s --M0 %s', cmd, args.M0_name);
            end
            
        end

    else
        %Set VFA files in the options string
        t1_str = sprintf('%s', args.T1_vols{1});
        for i_t = 2:length(args.T1_vols)
            t1_str = sprintf('%s,%s', t1_str, args.T1_vols{i_t});
        end
        cmd = sprintf('%s --T1_vols %s', cmd, t1_str);
        
        if args.T1_noise
            cmd = sprintf('%s --T1_noise %5.4f', cmd, args.T1_noise);
        end

    end 
    
    %Set any other options required to convert signal to concentration
    if isfinite(args.r1_const)
        cmd = sprintf('%s --r1 %4.3f', cmd, args.r1_const);
    end

    if isfinite(args.dose)
        cmd = sprintf('%s -D %4.3f', cmd, args.dose);
    end

    if args.M0_ratio
        cmd = sprintf('%s --M0_ratio', cmd);
    end
    
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

if args.dyn_noise
    cmd = sprintf('%s --dyn_noise', cmd);
end

if args.overwrite
    cmd = sprintf('%s --overwrite', cmd);
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

if ~isempty(args.aif_name)
    cmd = sprintf('%s --aif %s', cmd, args.aif_name);
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

if ~isempty(args.init_maps_dir)
    cmd = sprintf('%s --init_maps %s', cmd, init_maps_dir);
elseif ~isempty(args.init_params)
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
[status, result] = system(cmd);

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
[status, result] = run_madym(...
    'ETM', ...
    Ct_output_dir,...
    'dynamic_basename', [dyn_dir 'Ct_'],...
    'input_Ct', 1,...
    'output_Ct_sig',1,...
    'output_Ct_mod',1,...
    'injection_image', injection_img,...
    'overwrite', 1);

fprintf('status: %d\nresult: %s\n', status, result);

ktrans_fit = load_img_volume([Ct_output_dir 'ktrans.hdr']);
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
        load_img_volume([Ct_output_dir 'Ct_sig' num2str(i_dyn) '.hdr']);
    Cm_t(i_dyn) = ...
        load_img_volume([Ct_output_dir 'Ct_mod' num2str(i_dyn) '.hdr']);
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
[status, result] = run_madym(...
    'ETM', ...
    St_output_dir,...
    'dynamic_basename', [dyn_dir 'St_'],...
    'input_Ct', 0,...
    'output_Ct_sig',1,...
    'output_Ct_mod',1,...
    'T1_name', T1_name,...
    'r1_const', r1_const,...
    'injection_image', injection_img,...
    'overwrite', 1);

fprintf('status: %d\nresult: %s\n', status, result);

ktrans_fit = load_img_volume([St_output_dir 'ktrans.hdr']);
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
        load_img_volume([St_output_dir 'Ct_sig' num2str(i_dyn) '.hdr']);
    Cm_t(i_dyn) = ...
        load_img_volume([St_output_dir 'Ct_mod' num2str(i_dyn) '.hdr']);
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
output_files = dir(Ct_output_dir);
for i_f = 3:length(output_files)
    delete([Ct_output_dir output_files(i_f).name]);
end
rmdir(Ct_output_dir);
output_files = dir(St_output_dir);
for i_f = 3:length(output_files)
    delete([St_output_dir output_files(i_f).name]);
end
rmdir(St_output_dir);

%% ------------------------------------------------------------------------
% Below are the full C++ options for madym lite
%      madym options_:
%   -c [ --config ] arg (="")             Read input parameters from a 
%                                         configuration file
%   --cwd arg (="")                       Set the working directory
% 
% madym config options_:
%   --Ct arg (=0)                         Flag specifying input dynamic sequence 
%                                         are concentration (not signal) maps
%   -d [ --dyn ] arg (=dyn_)              Root name for dynamic volumes
%   --dyn_dir arg (="")                   Folder containing dynamic volumes, can 
%                                         be left empty if already included in 
%                                         option --dyn
%   --dyn_name_format arg (=%01u)        Number format for suffix specifying 
%                                         temporal index of dynamic volumes
%   -n [ --n_dyns ] arg (=0)              Number of DCE volumes, if 0 uses all 
%                                         images matching file pattern
%   -i [ --inj ] arg (=8)                 Injection image
%   --roi arg (="")                       Path to ROI map
%   -T [ --T1_method ] arg (=VFA)         Method used for baseline T1 mapping
%   --T1_vols arg (=[])                   _filepaths to input signal volumes (eg 
%                                         from variable flip angles)
%   --T1_noise arg (=0)                   Noise threshold for fitting baseline T1
%   --n_T1 arg (=0)                       Number of input signals for baseline T1
%                                         mapping
%   --M0_ratio arg (=1)                   Flag to use ratio method to scale 
%                                         signal instead of precomputed M0
%   --T1 arg (="")                        Path to precomputed T1 map
%   --M0 arg (="")                        Path to precomputed M0 map
%   --r1 arg (=3.3999999999999999)        Relaxivity constant of concentration in
%                                         tissue
%   --aif arg (="")                       Path to precomputed AIF, if not set 
%                                         uses population AIF
%   --pif arg (="")                       Path to precomputed AIF, if not set 
%                                         derives from AIF
%   --aif_slice arg (=0)                  Slice used to automatically measure AIF
%   -D [ --dose ] arg (=0.10000000000000001)
%                                         Contrast-agent dose
%   -H [ --hct ] arg (=0.41999999999999998)
%                                         Haematocrit level
%   -m [ --model ] arg (="")              Tracer-kinetic model
%   --init_params arg (=[ ])              Initial values for model parameters to 
%                                         be optimised
%   --init_maps arg (="")                 Path to folder containing parameter 
%                                         maps for per-voxel initialisation
%   --init_map_params arg (=[ ])          Index of parameters sampled from maps
%   --param_names arg (=[])               Names of model parameters, used to 
%                                         override default output map names
%   --fixed_params arg (=[ ])             Index of parameters fixed to their 
%                                         initial values
%   --fixed_values arg (=[ ])             _values for fixed parameters
%   --relative_limit_params arg (=[ ])    Index of parameters to which relative 
%                                         limits are applied
%   --relative_limit_values arg (=[ ])    _values for relative limits - optimiser 
%                                         bounded to range initParam +/- rel_limit
%   --first arg (=0)                      First image used in model fit cost 
%                                         function
%   --last arg (=0)                       Last image used in model fit cost 
%                                         function
%   --no_opt arg (=0)                     Flag to turn-off optimisation, default 
%                                         false
%   --dyn_noise arg (=0)                  Flag to use varying temporal noise in 
%                                         model fit, default false
%   --test_enh arg (=1)                   Flag to test for enhancement before 
%                                         fitting model, default true
%   --max_iter arg (=0)                   Max iterations per voxel in 
%                                         optimisation - 0 for no limit
%   --Ct_sig arg (=0)                     Flag to save signal-derived dynamic 
%                                         concentration maps
%   --Ct_mod arg (=0)                     Flag to save modelled dynamic 
%                                         concentration maps
%   -I [ --iauc ] arg (=[ 60 90 120 ])    _times (in s, post-bolus injection) at 
%                                         which to compute IAUC
%   -o [ --output ] arg (="")             Output folder
%   --overwrite arg (=0)                  Flag to overwrite existing analysis in 
%                                         output dir, default false
%   --sparse arg (=0)                     Flag to write output in sparse Analyze 
%                                         format
%   -E [ --err ] arg (=error_codes)       _filename of error codes map
%   --program_log arg (=ProgramLog.txt)   _filename of program log, will be 
%                                         appended with datetime
%   --config_out arg (=config.txt)        _filename of output config file, will be
%                                         appended with datetime
%   --audit arg (=AuditLog.txt)           _filename of audit log, will be appended
%                                         with datetime
%   --audit_dir arg (=C:\isbe\code\obj_msvc2015\manchester_qbi_public\bin\Release\audit_logs/)
%                                         Folder in which audit output is saved
% 
%   -h [ --help ]                         Print options and quit
%   -v [ --version ]                      Print version and quit
        