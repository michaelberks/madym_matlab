function [status, result] = run_madym_DicomConvert(varargin)
%RUN_MADYM_DICOMCONVERT MADYM wrapper function to call C++ tool madym_DicomConvert. Converts original dicom
%image slices into 3D Nifti/Analyze format images for input to the main Madym
%analysis tools, additionally create the XTR files of meta-information 
%required by Madym.
%   [status, result] = ...
%       run_madym_DicomConvert(varargin)
%
% RUN_MADYM_DICOMCONVERT uses the U_PACKARGS interface function
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
    'cmd_exe', [local_madym_root 'madym_DicomConvert'],...
    'config', '',... Path to a config file to set default options
    'output_dir', [], ...Output path, will use temp dir if empty;
    'T1_vols', [], ...,File names of generated signal volumes to from which baseline T1 is mapped
    'T1_dir', '', ...Folder to which T1 input volumes saved
    'dynamic_basename', '', ...Template name for dynamic sequences eg. dynamic/dyn_
    'dyn_dir', '',...Dir to dynamic series (can be included in dynamic basename)
    'sequence_format', '', ...Format for converting dynamic series index to string, eg %01u
    'sequence_start', NaN, ...Start index for dynamic series file names
    'sequence_step', NaN, ...Step between indexes of filenames in dynamic series
    'n_dyns', NaN, ...Number of dynamic sequence maps to load. If <=0, loads all maps in dynamic dir matching -dyn pattern
    'img_fmt_w', '', ...Image format for writing output
    'dicom_dir', '', ...Folder containing DICOM data
    'dicom_series_file', '', ...Filename to/from which dicom series data is written/read
    'T1_input_series', [], ...Indices of the dicom series for each T1 input
    'dyn_series', NaN, ...,Index of the dicom series for the dynamic DCE time-series
    'single_series', NaN, ...,Index of the dicom series for converting a generic single volume
    'dicom_filter', '', ...File filter for dicom sort (eg IM_)
    'vol_name', '', ...Output filename for converting a single dicom volume
    'sort', NaN, ...,Sort the files in Dicom dir into separate series, writing out the series information
    'make_t1', NaN, ...,Make T1 input images from dicom series
    'make_single', NaN, ...,Make single 3D image from dicom series
    'make_dyn', NaN, ...,Make dynamic images from dynamic series
    'make_t1_means', NaN, ...,Make mean of each set of T1 input repeats
    'make_dyn_mean', NaN, ...,Make temporal mean of dynamic images
    'flip_x', NaN, ...,Flip dicom slices horizontally before copying into 3D image volume
    'flip_y', NaN, ...,Flip dicom slices vertically before copying into 3D image volume
    'scale_tag', '', ...Dicom tag key (group,element) for rescale slope, in hexideciaml form - for Philips this is (0x2005, 0x100e)
    'offset_tag', '', ...Dicom tag key (group,element) for rescale intercept, in hexideciaml form - for Philips this is (0x2005, 0x100d)
    'dicom_scale', NaN, ...,Additional scaling factor applied to the dicom data
    'dicom_offset', NaN, ...,Additional offset factor applied to the dicom data
    'acquisition_time_tag', '', ...Dicom tag key (group,element) for acquisition time, if empty uses DCM_AcquisitionTime
    'temp_res', NaN, ...,Time in seconds between volumes in the DCE sequence, used to fill acquisition time not set in dynTimeTag
    'repeat_prefix', '', ...Prefix of image name for repeats in DICOM series, appended with sequence_format index and stored in series name folder
    'mean_suffix', '', ...Suffix of image name for mean of repeats in DICOM series, appended to series name
    'overwrite', NaN, ...,Set overwrite existing analysis in output dir
    'program_log_name', '', ...Program log file name
    'audit_dir', '', ...Folder in which audit output is saved
    'audit_name', '', ...Audit file name
    'config_out', '', ...Error codes image file name
    'error_name', '', ...Filename of output config file, will be appended with datetime
    'no_log', NaN, ...,Switch off program logging
    'no_audit', NaN, ...,Switch off audit logging
    'quiet', NaN, ...,Do not display logging messages in cout
    'working_directory', '', ...Sets the current working directory for the system call, allows setting relative input paths for data
    'dummy_run', false); %Don't run any thing, just print the cmd we'll run to inspect

%Set up base command
cmd = sprintf('%s ', args.cmd_exe);
    
%Set all the other commands
cmd = add_option('string', cmd, '--config', args.config);

cmd = add_option('string', cmd, '--cwd', args.working_directory);

cmd = add_option('string', cmd, '-o', args.output_dir);

cmd = add_option('string', cmd, '-d', args.dynamic_basename);

cmd = add_option('string', cmd, '--dyn_dir', args.dyn_dir);

cmd = add_option('string', cmd, '--sequence_format', args.sequence_format);

cmd = add_option('int', cmd, '--sequence_start', args.sequence_start);

cmd = add_option('int', cmd, '--sequence_step', args.sequence_step);

cmd = add_option('int', cmd, '--n_dyns', args.n_dyns);

cmd = add_option('string', cmd, '--img_fmt_w', args.img_fmt_w);

cmd = add_option('string_list', cmd, '--T1_vols', args.T1_vols);

cmd = add_option('string', cmd, '--T1_dir', args.T1_dir);

cmd = add_option('string', cmd, '--dicom_dir', args.dicom_dir);

cmd = add_option('string', cmd, '--dicom_series_file', args.dicom_series_file);

cmd = add_option('int_list', cmd, '--T1_series', args.T1_input_series);

cmd = add_option('int', cmd, '--dyn_series', args.dyn_series);

cmd = add_option('int', cmd, '--single_series', args.single_series);

cmd = add_option('string', cmd, '--dicom_filter', args.dicom_filter);

cmd = add_option('string', cmd, '--vol_name', args.vol_name);

cmd = add_option('bool', cmd, '--sort', args.sort);

cmd = add_option('bool', cmd, '--make_t1', args.make_t1);

cmd = add_option('bool', cmd, '--make_single', args.make_single);

cmd = add_option('bool', cmd, '--make_dyn', args.make_dyn);

cmd = add_option('bool', cmd, '--make_t1_means', args.make_t1_means);

cmd = add_option('bool', cmd, '--make_dyn_mean', args.make_dyn_mean);

cmd = add_option('bool', cmd, '--flip_x', args.flip_x);

cmd = add_option('bool', cmd, '--flip_y', args.flip_y);

cmd = add_option('string', cmd, '--scale_tag', args.scale_tag);

cmd = add_option('string', cmd, '--offset_tag', args.offset_tag);

cmd = add_option('float', cmd, '--dicom_scale', args.dicom_scale);

cmd = add_option('float', cmd, '--dicom_offset', args.dicom_offset);

cmd = add_option('string', cmd, '--acquisition_time_tag', args.acquisition_time_tag);

cmd = add_option('float', cmd, '--temp_res', args.temp_res);

cmd = add_option('string', cmd, '--repeat_prefix', args.repeat_prefix);

cmd = add_option('string', cmd, '--mean_suffix', args.mean_suffix);

cmd = add_option('bool', cmd, '--overwrite', args.overwrite);

cmd = add_option('bool', cmd, '--no_audit', args.no_audit);

cmd = add_option('bool', cmd, '--no_log', args.no_log);

cmd = add_option('bool', cmd, '--quiet', args.quiet);

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
    

    

