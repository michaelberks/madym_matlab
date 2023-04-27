function [status, result] = run_madym_AIF(varargin)
%RUN_MADYM_AIF wrapper function to call C++ tool madym_AIF. Auto-detects an
%arterial input function from time-series of DCE-MRI image volumes
%   [status, result] = ...
%       run_madym_AIF(varargin)
%
% RUN_MADYM_AIF uses the U_PACKARGS interface function
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
    'nifti_scaling', NaN, ... If set, applies intensity scaling and offset when reading/writing NIFTI images
    'nifti_4D', NaN, ... If set, reads NIFTI 4D images for T1 mapping and dynamic inputs
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

cmd = add_option('bool', cmd, '--nifti_scaling', args.nifti_scaling);

cmd = add_option('bool', cmd, '--nifti_4D', args.nifti_4D);

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
        cmd = add_option('string_list', cmd, '--T1_vols', args.T1_vols);
        
        cmd = add_option('float', cmd, '--T1_noise', args.T1_noise);
        
        cmd = add_option('string', cmd, '--T1_method', args.T1_method);

    end 
    
    %Set any other options required to convert signal to concentration
    cmd = add_option('float', cmd, '--r1', args.r1_const);

    cmd = add_option('bool', cmd, '--M0_ratio', args.M0_ratio);
    
end

%B1 correction options
cmd = add_option('string', cmd, '--B1', args.B1_name);

cmd = add_option('bool', cmd, '--B1_correction', args.B1_correction);

cmd = add_option('float', cmd, '--B1_scaling', args.B1_scaling);

%Now go through all the other optional parameters, and if they've been set,
%set the necessary option flag in the cmd string
cmd = add_option('bool', cmd, '--overwrite', args.overwrite);

cmd = add_option('bool', cmd, '--no_audit', args.no_audit);

cmd = add_option('bool', cmd, '--no_log', args.no_log);

cmd = add_option('bool', cmd, '--quiet', args.quiet);

cmd = add_option('int', cmd, '-i', args.injection_image);

cmd = add_option('string', cmd, '--roi', args.roi_name);

cmd = add_option('string', cmd, '--aif_map', args.aif_map);

cmd = add_option('float', cmd, '--TR', args.TR);

cmd = add_option('string', cmd, '--aif_slices', args.aif_slices);

cmd = add_option('string', cmd, '--aif_x_range', args.aif_x_range);

cmd = add_option('string', cmd, '--aif_y_range', args.aif_y_range);

cmd = add_option('float', cmd, '--min_T1_blood', args.min_T1_blood);

cmd = add_option('float', cmd, '--peak_time', args.peak_time);

cmd = add_option('float', cmd, '--prebolus_noise', args.prebolus_noise);

cmd = add_option('int', cmd, '--prebolus_min_images', args.prebolus_min_images);

cmd = add_option('float', cmd, '--select_pct', args.select_pct);

cmd = add_option('string', cmd, '--program_log', args.program_log_name);

cmd = add_option('string', cmd, '--audit', args.audit_name);

cmd = add_option('string', cmd, '--audit_dir', args.audit_dir);

cmd = add_option('string', cmd, '--config_out', args.config_out);

cmd = add_option('string', cmd, '-E', args.error_name);


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