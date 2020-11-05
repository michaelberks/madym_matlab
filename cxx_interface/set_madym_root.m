function [] = set_madym_root(madym_root)
%SET_MADYM_ROOT sets customised path to Madym executables on users local machine
%   [] = local_madym_root(madym_root)
%
% Inputs:
%      madym_root (str-[]) - Path to the madym binaries on this machine. If
%      empty, will generate a UI folder select window to choose.
%
% Outputs:
%
%
% Example: set_madym_root('~code/cxx_bin/madym');
%
% Notes:
%
% See also: STARTUP, LOCAL_MADYM_ROOT
%
% Created: 27-Feb-2019
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 
if ~exist('madym_root', 'var') || isempty(madym_root)
    h = helpdlg('Use finder to select folder for your madym tools',...
        'Locate madym root');
    uiwait(h);
    madym_root = uigetdir();
    if ~madym_root
        return;
    end
    madym_root = [madym_root '/'];
else
    %Check root is filesep terminated
    if madym_root(end) ~= '\' && madym_root(end) ~= '/'
        madym_root = [madym_root '/'];
    end
end

%First check if madym root already exists
old_madym_root = local_madym_root(false);

if ~isempty(old_madym_root)
    if strcmpi(madym_root, old_madym_root)
        fprintf('Madym root already set to %s, nothing to do\n',...
            madym_root);
        return;
    else
        choice = questdlg(...
            sprintf('Replace %s as madym root?', old_madym_root), ...
            'Madym root already set', ...
            'Yes', 'No', 'Yes');
        if strcmpi(choice, 'No')
            fprintf('Leaving madym root already set as %s\n',...
                old_madym_root);
            return;
        end
    end
end

%Append setenv to startup (create startup if it doesn't exist) so it
%persists in future sesssions
if ~exist('startup.m', 'file')
    fid = fopen('startup.m', 'wt');
else
    fid = fopen(which('startup'), 'at');
end
fprintf(fid, '\nsetenv(''MADYM_ROOT'', ''%s'');\n', madym_root);
fclose(fid);

%Set the variable for this session
setenv('MADYM_ROOT', madym_root);
fprintf('New madym root set to %s\n', madym_root);