function [up_to_date] = check_madym_updates(qbi_share_path)
%CHECK_MADYM_UPDATES checks whether the local madym version matches the 
% latest version available from the share drive.
%   [] = check_madym_updates(qbi_share_path)
%
% Inputs:
%      qbi_share_path - path to the root of the shared QBI drive (eg Q:\)
%
% Outputs:
%      up_to_date (bool) - true if local version matched QBI share latest
%      version, false otherwise
%
% Example: [version] = latest_madym_version('Q:');
%
% Notes:
%
% See also: INSTALL_MADYM, LOCAL_MADYM_VERSION, LATEST_MADYM_VERSION
%
% Created: 01-May-2019
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester
if ~exist('qbi_share_path', 'var')
    qbi_share_path = [];
end
latest_version = latest_madym_version(qbi_share_path);
local_version = local_madym_version();

up_to_date = strcmpi(latest_version, local_version);

if up_to_date
    fprintf('Your madym version %s is up-to-date\n', local_version);
else
    fprintf(...
        'Your madym version %s does not match the latest available version %s\n',...
        local_version, latest_version);
    fprintf('Run install_madym to update\n');
end
    



