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
    'n_dyns', NaN,...Number of dynamic sequence maps to load. If <=0, loads all maps in dynamic dir matching -dyn pattern
    'input_Ct', NaN,...Flag specifying input dynamic sequence are concentration (not signal) maps
    'output_Ct_sig', NaN,...Flag requesting concentration (derived from signal) are saved to output
    'output_Ct_mod', NaN,...Flag requesting modelled concentration maps are saved to output
    'T1_name', '',...Path to T1 map
    'M0_name', '',...Path to M0 map
    'B1_name', '',...Path to B1 correction map
    'B1_correction', NaN, ... Apply B1 correction
    'B1_scaling', NaN, ... Scaling factor to use with B1 map
    'r1_const', NaN,...Relaxivity constant of concentration in tissue (in ms)
    'TR', NaN, ... TR of dynamic sequence, only required if T1 method is IR
    'M0_ratio', NaN,...Flag to use ratio method to scale signal instead of supplying M0
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
    'nifti_scaling', NaN, ... If set, applies intensity scaling and offset when reading/writing NIFTI images
    'nifti_4D', NaN, ... If set, reads NIFTI 4D images for T1 mapping and dynamic inputs
    'IAUC_times', [],..._times (in s) at which to compute IAUC values
    'IAUC_at_peak', NaN,...Flag requesting IAUC computed at peak signal
    'param_names', [],...Names of model parameters to be optimised, used to name the output parameter maps
    'init_params', [],...Initial values for model parameters to be optimised, either as single vector, or 2D array NSamples x N_params
    'fixed_params', [],...Parameters fixed to their initial values (ie not optimised)
    'fixed_values', [],... Values for fixed parameters (overrides default initial parameter values)
    'repeat_param', NaN,...Index of parameter at which repeat fits will be made
    'repeat_values', [],... Values for repeat parameter
    'upper_bounds', [],...Upper bounds for each parameter during optimisation
    'lower_bounds', [],...Lower bounds for each parameter during optimisation
    'relative_limit_params', [],...Parameters with relative limits on their optimisation bounds
    'relative_limit_values', [],..._values for relative bounds, sets lower/upper bound as init param -/+ relative limit
    'residuals', '',... Path to existing residuals map
    'init_maps_dir', '',...Path to directory containing maps of parameters to initialise fit (overrides init_params)
    'init_map_params', [],...Parameters initialised from maps (if empty and init_maps_dir set, all params from maps)
    'dyn_noise', NaN,...Set to use varying temporal noise in model fit
    'test_enhancement', NaN,...Set test-for-enhancement flag
    'voxel_size_warn_only', NaN,... If set, only logs warning when voxel sizes don't match
    'max_iter', NaN,... Maximum number of iterations in model fit
    'opt_type', '',... Type of optimisation to run
    'img_fmt_r', '',...Set image read format
    'img_fmt_w', '',...Set image write format
    'overwrite', NaN,...Set overwrite existing analysis in output dir
    'no_audit', NaN,... Turn off audit log
    'no_log', NaN,... Turn off propgram log
    'quiet', NaN,... Suppress output to stdout
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
cmd = add_option('string', cmd, '--config', args.config);

%Set the working dir
cmd = add_option('string', cmd, '--cwd', args.working_directory);

%Set TK model
cmd = add_option('string', cmd, '-m', args.model);

%Set output directory
cmd = add_option('string', cmd, '-o', args.output_dir);

%Set the dynamic names
cmd = add_option('string', cmd, '-d', args.dynamic_basename);

cmd = add_option('string', cmd, '--dyn_dir', args.dyn_dir);

%Set the format
cmd = add_option('string', cmd, '--sequence_format', args.sequence_format);

cmd = add_option('int', cmd, '--sequence_start', args.sequence_start);

cmd = add_option('int', cmd, '--sequence_step', args.sequence_step);

%Set the number of dynamics
cmd = add_option('int', cmd, '--n_dyns', args.n_dyns);

%Set image formats
cmd = add_option('string', cmd, '--img_fmt_r', args.img_fmt_r);

cmd = add_option('string', cmd, '--img_fmt_w', args.img_fmt_w);

%Now set any args that require option inputs
cmd = add_option('bool', cmd, '--Ct', args.input_Ct);

%Set T1 options
cmd = add_option('string', cmd, '--T1', args.T1_name);

cmd = add_option('string', cmd, '--M0', args.M0_name);

cmd = add_option('string', cmd, '--T1_method', args.T1_method); 

cmd = add_option('string_list', cmd, '--T1_vols', args.T1_vols);

cmd = add_option('float', cmd, '--T1_noise', args.T1_noise);


%Set any other options required to convert signal to concentration
cmd = add_option('float', cmd, '--r1', args.r1_const);

cmd = add_option('float', cmd, '--TR', args.TR);

cmd = add_option('float', cmd, '-D', args.dose);

cmd = add_option('bool', cmd, '--M0_ratio', args.M0_ratio);

%B1 correction options
cmd = add_option('bool', cmd, '--B1', args.B1_name);   

cmd = add_option('bool', cmd, '--B1_correction', args.B1_correction);

cmd = add_option('float', cmd, '--B1_scaling', args.B1_scaling);    

%Now go through all the other optional parameters, and if they've been set,
%set the necessary option flag in the cmd string
cmd = add_option('bool', cmd, '--no_opt', args.no_optimise);

cmd = add_option('bool', cmd, '--Ct_sig', args.output_Ct_sig);

cmd = add_option('bool', cmd, '--Ct_mod', args.output_Ct_mod);

cmd = add_option('bool', cmd, '--test_enh', args.test_enhancement);

cmd = add_option('int', cmd, '--max_iter', args.max_iter);

cmd = add_option('string', cmd, '--opt_type', args.opt_type);

cmd = add_option('bool', cmd, '--dyn_noise', args.dyn_noise);

cmd = add_option('bool', cmd, '--overwrite', args.overwrite);

cmd = add_option('bool', cmd, '--no_audit', args.no_audit);

cmd = add_option('bool', cmd, '--no_log', args.no_log);

cmd = add_option('bool', cmd, '--quiet', args.quiet);

cmd = add_option('float', cmd, '-H', args.hct);

cmd = add_option('int', cmd, '-i', args.injection_image);

cmd = add_option('int', cmd, '--first', args.first_image);

cmd = add_option('int', cmd, '--last', args.last_image);

cmd = add_option('string', cmd, '--roi', args.roi_name);

cmd = add_option('string', cmd, '--residuals', args.residuals);

cmd = add_option('string', cmd, '--aif', args.aif_name);

cmd = add_option('string', cmd, '--aif_map', args.aif_map);

cmd = add_option('string', cmd, '--pif', args.pif_name);

cmd = add_option('string', cmd, '--program_log', args.program_log_name);

cmd = add_option('string', cmd, '--audit', args.audit_name);

cmd = add_option('string', cmd, '--audit_dir', args.audit_dir);

cmd = add_option('string', cmd, '--config_out', args.config_out);

cmd = add_option('string', cmd, '-E', args.error_name);

cmd = add_option('float_list', cmd, '--iauc', args.IAUC_times);

cmd = add_option('bool', cmd, '--iauc_peak', args.IAUC_at_peak);

cmd = add_option('bool', cmd, '--nifti_scaling', args.nifti_scaling);

cmd = add_option('bool', cmd, '--nifti_4D', args.nifti_4D);

cmd = add_option('string', cmd, '--init_maps', args.init_maps_dir);

cmd = add_option('int_list', cmd, '--init_map_params', args.init_map_params);

cmd = add_option('float_list', cmd, '--init_params', args.init_params);

cmd = add_option('string_list', cmd, '--param_names', args.param_names);

cmd = add_option('int_list', cmd, '--fixed_params', args.fixed_params);
cmd = add_option('float_list', cmd, '--fixed_values', args.fixed_values);

cmd = add_option('int', cmd, '--repeat_param', args.repeat_param);
cmd = add_option('float_list', cmd, '--repeat_values', args.repeat_values);

cmd = add_option('float_list', cmd, '--upper_bounds', args.upper_bounds);
cmd = add_option('float_list', cmd, '--lower_bounds', args.lower_bounds);

cmd = add_option('int_list', cmd, '--relative_limit_params', ...
    args.relative_limit_params);
cmd = add_option('float_list', cmd, '--relative_limit_values', ...
    args.relative_limit_values);

cmd = add_option('bool', cmd, '--voxel_size_warn_only',...
    args.voxel_size_warn_only);

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