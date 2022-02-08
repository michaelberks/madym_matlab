function [status, result] = run_madym_MakeXtr(varargin)
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
    'cmd_exe', [local_madym_root 'madym_MakeXtr'],...
    'config', '',... Path to a config file to set default options
    'T1_method', '',...T1 method to use to fit, see notes for options
    'T1_vols', [], ...,File names of generated signal volumes to from which baseline T1 is mapped
    'T1_dir', '', ...Folder to which T1 input volumes saved
    'DWI_vols', [], ...,File names of generated signal volumes for DWI models
    'DWI_dir', '', ...Folder to which DWI input volumes saved
    'dynamic_basename', '', ...Template name for dynamic sequences eg. dynamic/dyn_
    'dyn_dir', '',...Dir to dynamic series (can be included in dynamic basename)
    'sequence_format', '', ...Format for converting dynamic series index to string, eg %01u
    'sequence_start', NaN, ...Start index for dynamic series file names
    'sequence_step', NaN, ...Step between indexes of filenames in dynamic series
    'n_dyns', NaN, ...Number of dynamic sequence maps to load. If <=0, loads all maps in dynamic dir matching -dyn pattern 
    'make_t1', NaN, ...,Make T1 input images from dicom series    
    'make_dyn', NaN, ...,Make dynamic images from dynamic series    
    'temp_res', NaN, ...,Time in seconds between volumes in the DCE sequence, used to fill acquisition time not set in dynTimeTag
    'TR', NaN, ... TR of dynamic series (in ms)
    'FA', NaN, ... Flip angle of dynamic series (degrees)
    'VFAs', [], ... Flip angle of T1 inputs (degrees)
    'TIs', [], ... Inversion times of T1 inputs (msecs)
    'Bvalues', [], ... B-values of DWI model inputs
    'dyn_times', [], ... Time associated with each dynamic signal (in mins), must be supplied if using population AIF
    'working_directory', '', ...Sets the current working directory for the system call, allows setting relative input paths for data
    'dummy_run', false); %Don't run any thing, just print the cmd we'll run to inspect

%Set up base command
cmd = sprintf('%s ', args.cmd_exe);
    
%Set all the other commands
cmd = add_option('string', cmd, '--config', args.config);

cmd = add_option('string', cmd, '--cwd', args.working_directory);

cmd = add_option('string', cmd, '--T1_method', args.T1_method);

cmd = add_option('string_list', cmd, '--T1_vols', args.T1_vols);

cmd = add_option('string', cmd, '--T1_dir', args.T1_dir);

cmd = add_option('string_list', cmd, '--DWI_vols', args.DWI_vols);

cmd = add_option('string', cmd, '--DWI_dir', args.DWI_dir);

cmd = add_option('string', cmd, '-d', args.dynamic_basename);

cmd = add_option('string', cmd, '--dyn_dir', args.dyn_dir);

cmd = add_option('string', cmd, '--sequence_format', args.sequence_format);

cmd = add_option('int', cmd, '--sequence_start', args.sequence_start);

cmd = add_option('int', cmd, '--sequence_step', args.sequence_step);

cmd = add_option('int', cmd, '--n_dyns', args.n_dyns);

cmd = add_option('bool', cmd, '--make_t1', args.make_t1);

cmd = add_option('bool', cmd, '--make_dyn', args.make_dyn);

cmd = add_option('float', cmd, '--temp_res', args.temp_res);

cmd = add_option('float', cmd, '--TR', args.TR);

cmd = add_option('float', cmd, '--FA', args.FA);

cmd = add_option('float_list', cmd, '--VFAs', args.VFAs);

cmd = add_option('float_list', cmd, '--TIs', args.TIs);

cmd = add_option('float_list', cmd, '--Bvalues', args.Bvalues);

if ~isempty(args.dyn_times)
    %Get a name for the temporary file we'll write times to (we'll hold
    %off writing anything until we know this isn't a dummy run
    dyn_times_file = tempname;
    fid = fopen(dyn_times_file, 'wt');
    fprintf(fid, '%6.5f ', args.dyn_times(:));
    fclose(fid);
    cmd = add_option('string', cmd, '-t', dyn_times_file);
end

if args.dummy_run
    %Don't actually run anything, just print the command
    status = []; result = [];
    fprintf('%s\n', cmd);  
    return;  
end

%Otherwise we can run the command:
[status, result] = system(cmd, '-echo');

%Tidy up dyn times file
if exist('dyn_times_file', 'var')
    delete(dyn_times_file);
end

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
    

