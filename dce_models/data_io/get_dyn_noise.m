function [dyn_noise] = get_dyn_noise(root_path, num_vols, index_fmt)
%GET_DYN_TIMES get volume specific noise parameter from xtr files in folder
%of dynamic volumes
%   [times] = get_dyn_times(root_path, index_fmt, num_vols)
%
% Inputs:
%      root_path - folder + filename root where volumes are
%
%      num_vols - number of volumes to load
%
%      index_fmt ('%01u') - format that converts indexes into volume suffix
%
%
% Outputs:
%      dyn_noise - N_times x 1 array of dynamic noise values
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

dyn_noise = nan(num_vols,1);

for i_vol = 1:num_vols
    xtr_path = [root_path sprintf(index_fmt, i_vol) '.xtr'];
    [fields, values] = read_xtr_file(xtr_path);
    noise_idx = strcmpi(fields, 'NoiseSigma');
    
    if ~any(noise_idx)
        warning(['No NoiseSigma field in ' xtr_path...
            ', returned values will contain NaNs']);
    else
        noise_vals = values(noise_idx);
        if length(noise_vals) > 1
            warning(['Multiple NoiseSigma field in ' xtr_path...
            ', using final value']);
        end
        dyn_noise(i_vol) = noise_vals(end);
    end
end