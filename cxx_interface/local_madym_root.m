function [madym_root] = local_madym_root(empty_check)
%LOCAL_MADYM_ROOT customises a path based on users local machine
%   [madym_root] = local_madym_root
%
% Inputs:
%       empty_check (bool - true) - if true warns if madym root not set
%
% Outputs:
%      madym_root - Path to the madym binaries on this machine
%
%
% Example: cmd_exe = [local_madym_root 'madym_lite'];
%
% Notes: Please set your local root in an environment variable MADYM_ROOT.
% You can do this either by setting a system (or user) wide environment 
% variable, or (the easier and preferred way), adding a line
% SETENV('MADYM_ROOT', '<replace_with_your_local_path>'); into your Matlab
% startup script. You don't have a startup script? Set one up now!! ;o)
%
% See also: STARTUP
%
% Created: 27-Feb-2019
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 
if ~exist('empty_check', 'var') || isempty(empty_check)
    empty_check = true;
end

madym_root = getenv('MADYM_ROOT');
if empty_check && isempty(madym_root)
    warning(['MADYM_ROOT environment variable not set. '...
        'Add "setenv(''MADYM_ROOT'' your_local_root)" in your startup to set']);
end

%Deal with spaces by quote enclosing
if any(madym_root == ' ')
    %Get rid of existing quotes then add new
    madym_root(madym_root=='"' | madym_root == '''') = '';
    madym_root = ['"' madym_root '"'];
end

if ~isempty(madym_root) && madym_root(end) ~= '\' && madym_root(end) ~= '/'
    madym_root = [madym_root '/'];
end

            
