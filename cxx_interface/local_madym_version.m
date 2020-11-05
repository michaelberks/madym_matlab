function [version] = local_madym_version(cmd_exe)
%LOCAL_MADYM_VERSION return version of madym installed on this machine
%   [version] = local_madym_version()
%
% Inputs:
%       cmd_exe (str - []) - can be set if non-default installation of
%       madym, otherwise defaults to madym_lite at local_madym_root
%
% Outputs:
%       version (str) - Madym version in (Major).(Minor).(Patch)
%
% Example:
%
% Notes:
%
% See also: LATEST_MADYM_VERSION, LOCAL_MADYM_ROOT
%
% Created: 01-May-2019
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 
if ~exist('cmd_exe', 'var') || isempty(cmd_exe)
    cmd_exe = [local_madym_root(false) 'madym_lite'];
end

v_file = tempname;
status = system([cmd_exe ' -version > ' v_file]);

if status
    version = [];
    return;
end


fid = fopen(v_file);
t = textscan(fid,'%s');
fclose(fid);
delete(v_file);
version = t{1}{1};
fprintf('Current Madym version on this computer is %s\n', version);
