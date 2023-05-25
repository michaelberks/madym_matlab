function [status, result] = run_madym_DWI(varargin)
%RUN_MADYM_DWI wrapper function to call C++ tool madym_DWI. Fits diffusion
%models (ADC or IVIM) to DWI 3D volumes
%   [status, result] = ...
%       run_madym_DWI(varargin)
%
% RUN_MADYM_DWI uses the U_PACKARGS interface function
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
if ~nargin
    run_test();
    return;
end

% Unpack the arguments:
args = u_packargs(varargin, 1, ...
    'cmd_exe', [local_madym_root 'madym_DWI'],...
    'config', '',... Path to a config file to set default options
    'output_dir', [], ...Output path, will use temp dir if empty;
  	'roi_name', '',...Path to ROI map
    'DWI_vols', [], ..._file names of signal volumes to from which baseline T1 is mapped
    'model', '',...Diffusion model to fit
    'Bvals_thresh', [], ...Thresholds used in IVIM fitting
    'img_fmt_r', '',...Set image read format
    'img_fmt_w', '',...Set image write format
    'nifti_scaling', NaN, ... If set, applies intensity scaling and offset when reading/writing NIFTI images
    'nifti_4D', NaN, ... If set, reads NIFTI 4D images for T1 mapping and dynamic inputs
    'use_BIDS', NaN, ... If set, writes images using BIDS json meta info
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

cmd = add_option('string', cmd, '--roi', args.roi_name);

%Set image formats
cmd = add_option('string', cmd, '--img_fmt_r', args.img_fmt_r);

cmd = add_option('string', cmd, '--img_fmt_w', args.img_fmt_w);

cmd = add_option('bool', cmd, '--nifti_scaling', args.nifti_scaling);

cmd = add_option('bool', cmd, '--nifti_4D', args.nifti_4D);

cmd = add_option('bool', cmd, '--use_BIDS', args.use_BIDS);

cmd = add_option('string_list', cmd, '--DWI_vols', args.DWI_vols);
        
cmd = add_option('string', cmd, '--DWI_model', args.model);

cmd = add_option('float_list', cmd, '--Bvals_thresh', args.Bvals_thresh);

%Now go through all the other optional parameters, and if they've been set,
%set the necessary option flag in the cmd string
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

%--------------------------------------------------------------------------

%% ------------------------------------------------------------------------
%Test function to run if no inputs
function run_test()

%Run tests for ADC and IVIM models
sigma = 1;
B_vals = [0.0, 20.0, 40.0, 60.0, 80.0, 100.0, 300.0, 500.0, 800.0];
nBvals = length(B_vals);

%Generate ADC test data with Rician noise added
S0 = 100;
ADC = 0.8e-3;
signals = ADC_model(B_vals, S0, ADC);
signals_n = add_rician_noise(signals, sigma);

%Write these signals to Analyze images
Bvals_names = cell(nBvals,1);
DWI_dir = [tempdir 'DWI_test/'];
create_folder(DWI_dir);
for i_b = 1:nBvals
    Bval_name = [DWI_dir 'Bval_' num2str(i_b)];
    Bvals_names{i_b} =  [Bval_name '.hdr'];
    xtr_name = [Bval_name '.xtr'];

    save_img_volume(signals(:,i_b), Bvals_names{i_b}, [1 1 1], [], [], 0);
    write_xtr_file(xtr_name, 0, ...
        'B', B_vals(i_b),...
        'TimeStamp', 120000.0);
end

%Use madym_DWI to fit this data
run_madym_DWI(...
    'DWI_vols', Bvals_names,... 
    'model', 'ADC',...DWI model to use to fit, see notes for options
    'output_dir', DWI_dir, ...
    'overwrite', 1,...Set overwrite existing analysis in output dir ON
    'img_fmt_r', 'ANALYZE',...
    'img_fmt_w', 'ANALYZE',...
    'dummy_run', 0,...
    'no_audit', 1);

S0_f = load_img_volume([DWI_dir '/S0.hdr']);
ADC_f = load_img_volume([DWI_dir '/ADC.hdr']);
signals_f = ADC_model(B_vals, S0_f, ADC_f );

%Delete the test data
delete_list = dir([DWI_dir '/*']);
for i_d = 3:length(delete_list)
    delete([DWI_dir '/' delete_list(i_d).name]);
end
rmdir(DWI_dir);

%Display plots of the fit
figure('Name', 'madym_DWI test applied');
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

%Write these signals to Analyze images
Bvals_names = cell(nBvals,1);
DWI_dir = [tempdir 'DWI_test/'];
create_folder(DWI_dir);
for i_b = 1:nBvals
    Bval_name = [DWI_dir 'Bval_' num2str(i_b)];
    Bvals_names{i_b} =  [Bval_name '.hdr'];
    xtr_name = [Bval_name '.xtr'];

    save_img_volume(signals(:,i_b), Bvals_names{i_b}, [1 1 1], [], [], 0);
    write_xtr_file(xtr_name, 0, ...
        'B', B_vals(i_b),...
        'TimeStamp', 120000.0);
end

%Use madym_DWI to fit this data
run_madym_DWI(...
    'DWI_vols', Bvals_names,... 
    'model', 'IVIM',...DWI model to use to fit, see notes for options
    'Bvals_thresh', Bvals_thresh,...
    'output_dir', DWI_dir, ...
    'overwrite', 1,...Set overwrite existing analysis in output dir ON
    'img_fmt_r', 'ANALYZE',...
    'img_fmt_w', 'ANALYZE',...
    'dummy_run', 0,...
    'no_audit', 1);

S0_f = load_img_volume([DWI_dir '/S0.hdr']);
D_f = load_img_volume([DWI_dir '/d.hdr']);
f_f = load_img_volume([DWI_dir '/f.hdr']);
D_star_f = load_img_volume([DWI_dir '/dstar.hdr']);
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