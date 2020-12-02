function [T1, M0, errorCodes, status, result] =...
    run_madym_T1(varargin)
%RUN_MADYM_T1 wrapper function to call C++ T1 calculator, applying the
%variable flip angle method. Inputs can be paths to analyze-format images,
%or numeric arrays
%   [model_params, model_fit, error_codes, model_conc, dyn_conc] = ...
%       run_madym_T1(model, inputData, varargin)
%
% RUN_MADYM_T1 uses the U_PACKARGS interface function
% and so arguments to the function may be passed as name-value pairs
% or as a struct with fields with corresponding names. The field names
% are defined as:
%
% Mandatory Arguments:
%
% Optional Arguments: This wrapper uses the upackargs interface to allow
%       optional arguments to be entered as name/value pairs. See below for
%       a full description of all options.
%
% Outputs:
%      
%      T1 (1D array, Nsamples x 1 or []) - if fitting to numeric data,
%      vector of T1 values computed for each input sample.
%
%      M0 (1D array, Nsamples x 1 or []) - if fitting to numeric data,
%      vector of M0 values computed for each input sample.
%
%      [status,result] - returned by the system call to the Madym executable.
%      These may be operating system dependent, however status=0 should
%      mean an error free operation. If status is non-zero an error has
%      probably occurred, the result of which may be set in result. In any
%      event, it is best to check the output_dir, and any program logs that
%      have been saved there.
%
% Examples:
%   Fitting to full volumes:
%
%   Fitting to numeric data:
%
%
% Notes:
%
%   All T1 methods implemented in the main MaDym and MaDym-Lite C++ tools are
%   available to fit. Currently only the variable flip angle method is available:
% 
%   "VFA"
%
% See also: RUN_MADYM, RUN_MADYM_LITE
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
    'cmd_exe', [local_madym_root 'madym_T1'],...
    'cmd_exe_lite', [local_madym_root 'madym_T1_lite'],...
    'T1_vols', [],... Cell array of variable flip angle file paths
    'FAs', ... FAs, either single vector used for all samples, or 2D array, 1 row per sample
    'signals', ...Signals associated with each FA, 1 row per sample
  	'TR', [],... TR in msecs, required if directly fitting (otherwise will be taken from FA map headers);
    'method', 'VFA',...T1 method to use to fit, see notes for options
    'output_dir', [], ...Output path, will use temp dir if empty;
  	'output_name', 'madym_analysis.dat', ... Name of output file
    'noise_thresh', 100, ... PD noise threshold
    'roi_name', [],...Path to ROI map
    'error_name', [],... Name of error codes image
    'overwrite', 0,...Set overwrite existing analysis in output dir ON
    'dummy_run', 0 ...Don't run any thing, just print the cmd we'll run to inspect
    );
clear varargin;

%%
if isempty(args.T1_vols) && (isempty(args.FAs) || isempty(args.signals))
    error('Must supply either map names, or both FA and signal data');
end

%Set up output directory
deleteOutput = 0;
if isempty(args.output_dir)
    args.output_dir = tempdir;
    deleteOutput = 1;
elseif args.output_dir(end) ~= '\' && args.output_dir(end) ~= '/'
    args.output_dir = [args.output_dir '/'];
end

%Check if fitting to full volumes saved on disk, or directly supplied data
if ~isempty(args.T1_vols)
    %Use calculate_T1 to fit full volumes
    use_lite = false;
    
    %Set up FA map names
    nFAs = length(args.T1_vols);
    if nFAs < 3
        error('Only %d FA maps supplied, require at least 3 for T1 fitting',...
            nFAs);
    end  
    fa_str = sprintf('%s', args.T1_vols{1});
    for i_t = 2:nFAs
        fa_str = sprintf('%s,%s', fa_str, args.T1_vols{i_t});
    end
   
    %Initialise command argument
    cmd = sprintf(...
        '%s -T %s --T1_vols %s -o %s',...
        args.cmd_exe,...
        args.method,...
        fa_str,...
        args.output_dir);
    
    if ~isempty(args.error_name)
        cmd = sprintf('%s -E %s', cmd, args.error_name);
    end
    if ~isempty(args.roi_name)
        cmd = sprintf('%s --roi %s', cmd, args.roi_name);
    end
    if args.overwrite
        cmd = sprintf('%s --overwrite', cmd);
    end
    cmd = sprintf('%s --T1_noise %5.4f', cmd, args.noise_thresh);
    
else
    %Fit directly supplied FA and signal data using calculate_T1_lite
    use_lite = true;
    [nSamples, nFAs] = size(args.signals);
    
    %Do error checking on required inputs
    if nFAs < 3
        error('Only %d FAs supplied, require at least 3 for T1 fitting',...
            nFAs);
    end
    
    if numel(args.FAs) == nFAs
        args.FAs = repmat(args.FAs(:)', nSamples, 1);
    elseif ~all(size(args.FAs)==[nSamples, nFAs])
        error('Size of FAs array does not match size of signals array');
    end
    
    %TR must be supplied
    if isempty(args.TR) || ~isnumeric(args.TR)
        error('You must supply a numeric TR value (in msecs) to fit directly to data');
    end
    
    %Set up temporary files for FAs and signals (we'll hold
    %off writing anything until we know this isn't a dummy run
    input_file = tempname;
    
    cmd = sprintf(...
        '%s -T %s --data %s --n_T1 %d --TR %4.3f -o %s -O %s',...
        args.cmd_exe_lite,...
        args.method,...
        input_file,...
        nFAs,...
        args.TR,...
        args.output_dir,...
        args.output_name);
    
    %Check for bad samples, these can screw up Madym as the lite version
    %of the software doesn't do the full range of error checks Madym proper
    %does. So chuck them out now and warn the user
    discard_samples = any(...
        isnan(args.FAs) |...
        ~isfinite(args.FAs) |...
        isnan(args.signals) |...
        ~isfinite(args.signals), 2);
    
    if any(discard_samples)
        warning(['Samples with NaN values found,'...
            'these will be set to zero for model-fitting']);
        args.FAs(discard_samples,:) = 0;
        args.signals(discard_samples,:) = 0;
    end    
    
end

if args.dummy_run
    %Don't actually run anything, just print the command
    fprintf('%s\n', cmd);
    
    T1 = [];
    M0 = [];
    errorCodes = [];
    status = [];
    result = [];
    return;
    
end

%For the lite method, no we can write the temporary files
if use_lite
    %Write input values to a temporary file
    fid = fopen(input_file, 'wt');
    for i_row = 1:nSamples
        fprintf(fid, '%6.5f ', args.FAs(i_row,:));
        fprintf(fid, '%6.5f ', args.signals(i_row,:));
        fprintf(fid, '\n');
    end
    fclose(fid);
end

%At last.. we can run the command
[status, result] = system(cmd);

if use_lite
    %Now load the output from calculate T1 lite and extract data to match this
    %functions outputs
    
    if nargout
        fullOutPath = [args.output_dir args.method '_' args.output_name];
        outputData = load(fullOutPath);
        T1 = outputData(:,1);
        M0 = outputData(:,2);
        errorCodes = outputData(:,3);
    end
    
    %Tidy up temporary files
    delete(input_file);   
    if deleteOutput
        delete(fullOutPath);
    end
else
    if nargout
        T1Path = [args.output_dir 'T1.hdr'];
        M0Path = [args.output_dir 'M0.hdr'];
        errorPath = [args.output_dir 'error_tracker.hdr'];
        T1 = load_img_volume(T1Path);
        M0 = load_img_volume(M0Path);
        try
            errorCodes = load_img_volume(errorPath);
        catch
            errorCodes = [];
        end
    end
    if deleteOutput
        delete([args.output_dir 'T1.hdr']);
        delete([args.output_dir 'T1.img']);
        delete([args.output_dir 'T1.xtr']);
        delete([args.output_dir 'M0.hdr']);
        delete([args.output_dir 'M0.img']);
        delete([args.output_dir 'M0.xtr']);
        delete([args.output_dir args.error_name, '.hdr']);
        delete([args.output_dir args.error_name, '.img']);
    end
    
end   

%%
%Test function to run if no inputs
function run_test()

    %Generate some signals from sample FA, TR, T1 and M0 values
    T1 = [500 1000 1500 500 1000 1500]';
    M0 = [1000 1000 1000 2000 2000 2000]';
    TR = 3.5;
    FAs = [2 10 18];
    
    signals = signal_from_T1(T1, M0, FAs, TR);
    
    %First run this in data mode using calculate_T1_lite:    
    [T1_fit, M0_fit] = run_madym_T1(...
        'FAs', FAs,... 
        'signals', signals,...
        'TR', TR,... 
        'method', 'VFA',...
        'dummy_run', 0);
    signals_fit = signal_from_T1(T1_fit, M0_fit, FAs, TR);
    
    figure('Name', 'calculate_T1_lite applied to numeric data');
    for i_sample = 1:6
        subplot(2,3,i_sample); hold all;
        plot(FAs, signals(i_sample,:), 'r*');
        plot(FAs, signals_fit(i_sample,:), 'go');
        plot(FAs, signals_fit(i_sample,:), 'b-');
        title({'Parameter estimates (actual,fit)';...
            sprintf('T1: (%d, %4.1f), M0: (%d, %4.1f)',...
            T1(i_sample), T1_fit(i_sample), M0(i_sample), M0_fit(i_sample))});
        
        if i_sample == 1
            legend({'Signals', 'Fit to signals',' '});
        end
        xlabel('Flip angle (degrees)');
        ylabel('Signal intensity');        
    end 
    
    %Now save the flip-angle data at Analyze images and apply the full
    %volume method
    FA_names = cell(3,1);
    T1_dir = [tempdir 'T1_test/'];
    create_folder(T1_dir);
    for i_fa = 1:3
        FA_name = [T1_dir 'FA_' num2str(i_fa)];
        FA_names{i_fa} =  [FA_name '.hdr'];
        xtr_name = [FA_name '.xtr'];

        save_img_volume(signals(:,i_fa), FA_names{i_fa}, [1 1 1], [], [], 0);
        write_xtr_file(xtr_name, 0, ...
            'FlipAngle', FAs(i_fa),...
            'TR', TR,...
            'TimeStamp', 120000.0);
    end
    roi_name = [tempdir 'ROI.hdr'];
    save_img_volume(ones(size(T1)), roi_name, [1 1 1]);
    
    [T1_fit, M0_fit] = run_madym_T1(...
        'T1_vols', FA_names,... 
        'method', 'VFA',...T1 method to use to fit, see notes for options
        'output_dir', T1_dir, ...
        'noise_thresh', 0, ... PD noise threshold
        'roi_name', [],...Path to ROI map
        'overwrite', 1,...Set overwrite existing analysis in output dir ON
        'dummy_run', 0);
    signals_fit = signal_from_T1(T1_fit, M0_fit, FAs, TR);
    
    for i_fa = 1:3
        delete(FA_names{i_fa});
        delete([FA_names{i_fa}(1:end-4) '.img']);
        delete([FA_names{i_fa}(1:end-4) '.xtr']);
    end  
    delete(roi_name);
    delete([roi_name(1:end-3) 'img']);
    
    figure('Name', 'calculate_T1 applied to Analyze volumes');
    for i_sample = 1:6
        subplot(2,3,i_sample); hold all;
        plot(FAs, signals(i_sample,:), 'r*');
        plot(FAs, signals_fit(i_sample,:), 'go');
        plot(FAs, signals_fit(i_sample,:), 'b-');
        title({'Parameter estimates (actual,fit)';...
            sprintf('T1: (%d, %4.1f), M0: (%d, %4.1f)',...
            T1(i_sample), T1_fit(i_sample), M0(i_sample), M0_fit(i_sample))});
        
        if i_sample == 1
            legend({'Signals', 'Fit to signals',' '});
        end
        xlabel('Flip angle (degrees)');
        ylabel('Signal intensity');        
    end 
    
    return;
    
%% -----------------------------------------------------------------------
% Below are the full C++ options for madym_T1 and madym_T1_lite
% 
% ********************* madym_T1 **************************************
%       madym_T1 options_:
%   -c [ --config ] arg (="")             Read input parameters from a 
%                                         configuration file
%   --cwd arg (="")                       Set the working directory
% 
% madym_T1 config options_:
%   --roi arg (="")                       Path to ROI map
%   -T [ --T1_method ] arg (=VFA)         Method used for baseline T1 mapping
%   --T1_vols arg (=[])                   _filepaths to input signal volumes (eg 
%                                         from variable flip angles)
%   --T1_noise arg (=0)                   Noise threshold for fitting baseline T1
%   --n_T1 arg (=0)                       Number of input signals for baseline T1
%                                         mapping
%   -o [ --output ] arg (="")             Output folder
%   --sparse arg (=0)                     Flag to write output in sparse Analyze 
%                                         format
%   --overwrite arg (=0)                  Flag to overwrite existing analysis in 
%                                         output dir, default false
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
%   -v [ --version ]
% 
% **************** madym_T1_lite **************************************
%     madym_T1_lite config options_:
%   --cwd arg (="")                       Set the working directory
%   --data arg (="")                      Input data filename, see notes for 
%                                         options
%   -T [ --T1_method ] arg (=VFA)         Method used for baseline T1 mapping
%   --FA arg (=0)                         FA of dynamic series
%   --TR arg (=0)                         TR of dynamic series
%   --T1_noise arg (=0)                   Noise threshold for fitting baseline T1
%   --n_T1 arg (=0)                       Number of input signals for baseline T1
%                                         mapping
%   -o [ --output ] arg (="")             Output folder
%   -O [ --output_name ] arg (=madym_analysis.dat)
%                                         Name of output data file
% 
%   -h [ --help ]                         Print options and quit
%   -v [ --version ]                      Print version and quit '''



