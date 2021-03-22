function [dyn_times, dyn_FA, dyn_TR] = get_dyn_times(root_path, num_vols, index_fmt)
%GET_DYN_TIMES get meta information (scan time, FA, TR) from xtr files for 
% folder of dynamic volumes
%   [times] = get_dyn_times(root_path, index_fmt, num_vols)
%
% Inputs:
%      root_path - folder + filename root where volumes are
%
%      num_vols - number of volumes to load
%
%
%      index_fmt ('%01u') - format that converts indexes into volume suffix
%
%
% Outputs:
%      dyn_times - *Insert description of input variable here*
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 29-Mar-2017
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 

if ~exist('index_fmt', 'var') || isempty(index_fmt)
    index_fmt = '%01u';
end

dyn_times = zeros(num_vols,1);
dyn_TR = zeros(num_vols,1);
dyn_FA = zeros(num_vols,1);

for i_vol = 1:num_vols
    
    vol_path = [root_path sprintf(index_fmt, i_vol) '.xtr'];
    [xtr_fields, xtr_values] = read_xtr_file(vol_path);
    
    %Would be more robust to check fields present than assume they are
    dyn_times(i_vol) =...
        madym_timestamp_to_secs(...
        xtr_values(strcmpi(xtr_fields, 'timestamp')));
    
    dyn_FA(i_vol) = xtr_values(strcmpi(xtr_fields, 'FlipAngle'));
    dyn_TR(i_vol) = xtr_values(strcmpi(xtr_fields, 'TR'));
    
end
dyn_times = (dyn_times - dyn_times(1))/60;