function [version] = latest_madym_version(qbi_share_path)
%LATEST_MADYM_VERSION return the most current madym version available a string. MB has
%responsibility for maintaining. Note this may NOT be the version you have
%installed locally. Call local_madym_version to check that.
%   [version] = latest_madym_version(qbi_share_path)
%
% Inputs:
%      qbi_share_path - path to the root of the shared QBI drive (eg Q:\).
%      If empty, you will be prompted to select via uigetdir
%
% Outputs:
%      version (str) - version in v(Major).(Minor).(Patch) form of latest
%      madym version available on the QBI shared drive
%
% Example: [version] = latest_madym_version('Q:');
%
% Notes:
%
% See also: LOCAL_MADYM_VERSION
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

if ispc
    madym_version_hdr = ...
        [qbi_share_path '/software/madym_versions/mdm_version_win.h'];
else
    madym_version_hdr = ...
        [qbi_share_path '/software/madym_versions/mdm_version.h'];
end

%Open version header and scan tags to get version
fid = fopen(madym_version_hdr);
t = textscan(fid, '%s');

ver_idx = find(strcmpi(t{1}, 'GIT_TAG')) + 1;
version = t{1}{ver_idx};

%Remove annoying "s
version(version=='"') = [];

fprintf('Latest Madym version on QBI share is %s\n', version);



