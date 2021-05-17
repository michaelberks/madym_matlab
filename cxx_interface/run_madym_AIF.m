function [status, result] = run_madym_AIF(varargin)
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
% or as a struct with fields with corresponding names. See below for
% a full description of all options.
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
%
% Notes:
% 
%
% See also: RUN_MADYM, RUN_MADYM_LITE, RUN_MADYM_T1
%
% Created: 20-Feb-2019
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 

% Unpack the arguments:
args = u_packargs(varargin, 1, ...
    'cmd_exe', [local_madym_root 'madym_AIF'],...
    'config', '',... Path to a config file to set default options
    'output_dir', [], ...Output path, will use temp dir if empty;
  	'T1_vols', [], ..._file names of signal volumes to from which baseline T1 is mapped
    'T1_method', '',...Method for mapping baseline T1 (eg VFA)
    'dynamic_basename', [], ...Template name for dynamic sequences eg. dynamic/dyn_
    'dyn_dir', '',...Dir to dynamic series (can be included in dynamic basename)
    'sequence_format', '',...Format for converting dynamic series index to string, eg %01u
    'sequence_start', NaN,...Start index of dynamic series
    'sequence_step', NaN,...Step size between indexes in dynamic series
    'n_dyns', 0,...Number of dynamic sequence maps to load. If <=0, loads all maps in dynamic dir matching -dyn pattern
    'input_Ct', true,...Flag specifying input dynamic sequence are concentration (not signal) maps
    'T1_name', '',...Path to T1 map
    'M0_name', '',...Path to M0 map
    'B1_name', '',...Path to B1 correction map
    'B1_correction', false, ... Apply B1 correction
    'B1_scaling', NaN, ... Scaling factor to use with B1 map
    'r1_const', NaN,...Relaxivity constant of concentration in tissue (in ms)
    'M0_ratio', true,...Flag to use ratio method to scale signal instead of supplying M0
    'injection_image', NaN,...Injection image
    'T1_noise', NaN,...PD noise threshold
    'roi_name', '',...Path to ROI map
    'aif_map', '',...Path to to mask from which AIF will be computed on the fly
    'TR', NaN,... TR of dynamic series
    'aif_slices', '', ... Slices used to automatically measure AIF
    'aif_x_range', '', ... Range of voxels to consider as AIF candidates in x-axis
  	'aif_y_range', '', ... Range of voxels to consider as AIF candidates in y-axis
  	'min_T1_blood', NaN, ... Minimum T1 to be considered as potential blood voxel
  	'peak_time', NaN,... Time window post bolus for peak to arrive
  	'prebolus_noise', NaN,... Estimate of noise on the prebolus signal used when unable to compute
  	'prebolus_min_images', NaN,... Minimum number of images required to estimate prebolus noise
  	'select_pct', NaN, ... Percentage of candidates to select
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
        
        if ~isempty(args.T1_method)
            cmd = sprintf('%s --T1_method %s', cmd, args.T1_method);
        end

    end 
    
    %Set any other options required to convert signal to concentration
    if isfinite(args.r1_const)
        cmd = sprintf('%s --r1 %4.3f', cmd, args.r1_const);
    end

    if args.M0_ratio
        cmd = sprintf('%s --M0_ratio', cmd);
    end
    
end

%B1 correction options
if ~isempty(args.B1_name)
    cmd = sprintf('%s --B1 %s', cmd, args.B1_name);
end    
if args.B1_correction
    cmd = sprintf('%s --B1_correction', cmd);
end
if isfinite(args.B1_scaling)
    cmd = sprintf('%s --B1_scaling %4.3f', cmd, args.B1_scaling);
end    

%Now go through all the other optional parameters, and if they've been set,
%set the necessary option flag in the cmd string
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

if isfinite(args.injection_image)
    cmd = sprintf('%s -i %d', cmd, args.injection_image);
end

if ~isempty(args.roi_name)
    cmd = sprintf('%s --roi %s', cmd, args.roi_name);
end

if ~isempty(args.aif_map)
    cmd = sprintf('%s --aif_map %s', cmd, args.aif_map);
end

if isfinite(args.TR)
    cmd = sprintf('%s --TR %4.3f', cmd, args.TR);
end

if ~isempty(args.aif_slices)
    cmd = sprintf('%s --aif_slices %s', cmd, args.aif_slices);
end

if ~isempty(args.aif_x_range)
    cmd = sprintf('%s --aif_x_range %s', cmd, args.aif_x_range);
end

if ~isempty(args.aif_y_range)
    cmd = sprintf('%s --aif_y_range %s', cmd, args.aif_y_range);
end

if isfinite(args.min_T1_blood)
    cmd = sprintf('%s --min_T1_blood %4.3f', cmd, args.min_T1_blood);
end

if isfinite(args.peak_time)
    cmd = sprintf('%s --peak_time %4.3f', cmd, args.peak_time);
end

if isfinite(args.prebolus_noise)
    cmd = sprintf('%s --prebolus_noise %4.3f', cmd, args.prebolus_noise);
end

if isfinite(args.prebolus_min_images)
    cmd = sprintf('%s --prebolus_min_images %d', cmd,...
        args.prebolus_min_images);
end

if isfinite(args.select_pct)
    cmd = sprintf('%s --select_pct %4.3f', cmd, args.select_pct);
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



if args.dummy_run
    %Don't actually run anything, just print the command
    status = []; result = [];
    fprintf('%s\n', cmd);  
    return;  
end

%Otherwise we can run the command:
[status, result] = system(cmd, '-echo');

%--------------------------------------------------------------------------


%% ------------------------------------------------------------------------