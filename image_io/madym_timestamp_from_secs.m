function [timestamp] = madym_timestamp_from_secs(secs)
%MADYM_TIMESTAMP_TO_SECS Helper function to convert time in secs into a 
% timestamp format used by madym
%   [secs] = madym_timestamp_to_secs(timestamp)
%
% Inputs:
%   secs - time in seconds from midnight
%
% Outputs: 
%   timestamp - a string - though stored as a double - hhmmss.msecs
%
%      
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
hh = floor(secs / (3600));
mm = floor((secs - 3600 * hh) / 60);
ss = secs - 3600 * hh - 60 * mm;
timestamp = 10000 * hh + 100 * mm + ss;