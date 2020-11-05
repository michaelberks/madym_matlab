function [] = install_madym(qbi_share_path, operating_system, madym_root, version, run_tests)
%INSTALL_MADYM install latest version of Madym C++ tools to your local
%machine. REQUIRES that you have set your local madym root first (see
%below)
%   [] = install_madym(qbi_share_path)
%
% Inputs:
%      qbi_share_path (str) - path to the root of the shared QBI drive (eg Q:\)
%
%      operating_system (str) - must match {'Windows', 'Ubuntu', MacOs'}
%
%      version (str - []) - version to install, if empty will install
%      latest version. If set, must match a version on the shared drive
%
%      version (int - 1) - flag to run tests on successful install. Set 0
%      not to run tests. _values >1 run advanced tests (currently not
%      implemented)
%
% Outputs: None
%
% Example: install_madym('Q:', 'Windows');
%          install_madym('Volumes/qbi', 'MacOs');
%
% Notes:
%
% See also: RUN_MADYM_TESTS, LOCAL_MADYM_ROOT, LATEST_MADYM_VERSION,
% LOCAL_MADYM_VERSION
%
% Created: 01-May-2019
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 
if ~exist('qbi_share_path', 'var') || isempty(qbi_share_path)
    h = helpdlg('Use finder to select root of share drive',...
        'Step 1: locate share drive');
    uiwait(h);
    
    qbi_share_path = uigetdir();
    if ~qbi_share_path
        return;
    end
end

if ~exist('operating_system', 'var') || isempty(operating_system)
    operating_system = questdlg('What is your operating system?', ...
                         'Choose OS', ...
                         'MacOS', 'Windows', 'Ubuntu', 'MacOS');
    if isempty(operating_system)
        return;
    end
end
if ~ismember(operating_system, {'Windows', 'Ubuntu', 'MacOS'})
    error('Operating systm %s not recognised. Must match Windows, Ubuntu or MacOS',...
        operating_system);
end

if ~exist('madym_root', 'var') || isempty(madym_root)
    madym_root = local_madym_root(false);    
end

%Set the madym root - if it's empty above, the set function will ask the
%user to choose it.
set_madym_root(madym_root);

%Finally, call local_madym_root again, as we dont know if in the set
%function above, the user may have chosen not to overwrite an existing
%madym root
madym_root = local_madym_root(false);

if ~exist('version', 'var') || isempty(version)
    version = latest_madym_version(qbi_share_path);
end
if isempty(version)
    error(['Could not retrieve madym version from the share drive.',...
        ' Check the shared paths: %s.'], qbi_share_path);
end

if ~exist('run_tests', 'var') || isempty(run_tests)
    run_tests = 1;
end

%Make path to the latest madym version and get file list
madym_path = ...
    [qbi_share_path '/software/madym_versions/' version '/' operating_system '/'];
madym_files = dir(madym_path);

%Copy files
create_folder(madym_root);
for i_f = 3:length(madym_files)
    copyfile([madym_path madym_files(i_f).name],...
        [madym_root madym_files(i_f).name]);
    fprintf('Copied %s from %s to %s\n',...
        madym_files(i_f).name, madym_path, madym_root);
end

%Check version now matches
local_version = local_madym_version();

if strcmpi(local_version, version)
    fprintf('*****************************************************\n');
    fprintf('Successfully installed version %s on this machine\n',...
        version);
    fprintf('*****************************************************\n');
else
    fprintf('*****************************************************\n');
    fprintf(['Problem installing: local version returned %s, '...
        'this suggests version %s from the share drive did not install\n'],...
        local_version, version);
    fprintf('*****************************************************\n');
end

%Run tests if set
if run_tests
    run_madym_tests(run_tests);
end
    

