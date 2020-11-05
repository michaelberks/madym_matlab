function [dyn_signals, signal_baseline] = get_dyn_vals(root_path, num_vols, roi, index_fmt, num_baseline_pts)
%GET_DYN_VALS given directory of volumes and ROI mask, get array of
%time-series for voxels in the ROI
%   [times] = get_dyn_vals(root_path, num_vols, roi, index_fmt)
%
% Inputs:
%      root_path - folder + filename root where volumes are
%
%      num_vols - number of volumes to load
%
%      roi - mask volume, must have same dims as dynamic volumes
%
%      index_fmt ('%01u') - format that converts indexes into volume suffix
%
%      num_baseline_pts (5) - number of volumes use to compute baseline mean
%      for each voxel
%
%
% Outputs:
%      dyn_signals - N_vox x N_times array of time-series for each voxel
%
%      signal_baseline - N_vox x 1 mean baseline signal for each voxel
%
%
% Example: [dyn_signals, signal_baseline] = ...
%               get_dyn_vals('dynamics/dyn_', 50, roi_mask)
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
if ~exist('num_baseline_pts', 'var') || isempty(num_baseline_pts)
    num_baseline_pts = 5;
end

%If ROI is a path, load it from disk
if ischar(roi)
    roi = load_img_volume(roi) > 0;
end

num_pts = nnz(roi);
dyn_signals = zeros(num_pts, num_vols);

load_vols = ischar(root_path);
for i_vol = 1:num_vols
    if load_vols
        vol_path = [root_path sprintf(index_fmt, i_vol) '.hdr'];
        vol = load_img_volume(vol_path);
    else
        vol = root_path(:,:,:,i_vol);
    end
    dyn_signals(:,i_vol) = vol(roi);
end

if nargout > 1
    signal_baseline =  mean(dyn_signals(:,1:num_baseline_pts),2);
end
