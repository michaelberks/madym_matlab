function [secs] = madym_timestamp_to_secs(timestamp)
%MADYM_TIMESTAMP_TO_SECS Helper function to convert time in secs into a 
% timestamp format used by madym
%   [secs] = madym_timestamp_to_secs(timestamp)
%
% Inputs:
%      timestamp - a string - though stored as a double - hhmmss.msecs
%
%
% Outputs:
%      secs - time in seconds from midnight
%
%
% Example:
%
% Notes:
%
% See also: MADYM_TIMESTAMP_FROM_SECS
%
% Created: 16-Mar-2021
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 
hh = floor(timestamp / 10000);
mm = floor( (timestamp - 10000 * hh) / 100);
ss = (timestamp ...
    - 10000 * hh ...
    - 100 * mm);
secs = hh * 60 * 60 ...
    + double(mm) * 60 ...
    + ss;