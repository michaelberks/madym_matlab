function [iauc] = compute_IAUC(dyn_concentrations, dyn_times, aif_injection_image, num_secs, convert_to_secs, enhancing_mask)
%COMPUTE_IAUC integrate dynamic series concentrations to compute IAUC
%   [iauc] = compute_IAUC(dyn_concentrations, dyn_times, aif_injection_image, end_time)
%
% Inputs:
%      dyn_concentrations - Array of concentrations in dynamic series,
%      assume time sequence is along final dimension as input can either
%      be n_voxels x n_times or n_y x n_x n_z x_ n_times
%
%      dyn_times - vector of times for each volume in dynamic series,
%      length much match size of final dimension in dyn_concentrations
%
%      aif_injection_image - index of image at point bolus was injected
%
%      num_secs - length of time from injection image over which to
%      integrate concentrations
%
%      convert_to_secs - scalar to convert dyn_times to seconds (as by
%      default these are stored in minutes, so convert_to_secs = 60). If
%      instead dyn_times is already is seconds, set convert_to_secs to 1
%      (or 1000 if msecs etc)
%
%      enhancing_mask - mask of voxels that enhanced, IAUC is set to zero
%      otherwise. Default all voxels enhancing
%
%
% Outputs:
%      iauc - Integrated area under concentration curve, will be an array
%      matching the size of the first ndims-1 dimension of
%      dyn_concentrations
%
%
% Example: [iauc90] = compute_IAUC(C_t, dyn_times, 8, 90);
%
% Notes: Uses discrete trapezeoid rule to integrate
%
% See also:
%
% Created: 31-Oct-2017
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester
if ~exist('convert_to_secs', 'var') || isempty(convert_to_secs)
    convert_to_secs = 60;
end
n_dims = ndims(dyn_concentrations);

if length(dyn_times) ~= size(dyn_concentrations, n_dims)
    error('Length of dyn_times vector does not much number of volumes in dyn_concentrations array');
end

dyn_times = dyn_times * convert_to_secs;
iauc_idx = find((dyn_times - dyn_times(aif_injection_image)) > num_secs ,1);

s_idx = aif_injection_image;
switch n_dims
    
    case 2
        time_intervals = reshape(...
            dyn_times(s_idx+1:iauc_idx) - dyn_times(s_idx:iauc_idx-1),...
            1, []);
         iauc = sum( bsxfun(@times,...
            0.5*dyn_concentrations(:,s_idx:iauc_idx-1) + ...
            0.5*dyn_concentrations(:,s_idx+1:iauc_idx) , ...
            time_intervals), 2);
    case 3
        time_intervals = reshape(...
            dyn_times(s_idx+1:iauc_idx) - dyn_times(s_idx:iauc_idx-1),...
            1, 1, []);
        iauc = sum( bsxfun(@times,...
            0.5*dyn_concentrations(:,:,s_idx:iauc_idx-1) + ...
            0.5*dyn_concentrations(:,:,s_idx+1:iauc_idx) , ...
            time_intervals), 3);
    case 4
        time_intervals = reshape(...
            dyn_times(s_idx+1:iauc_idx) - dyn_times(s_idx:iauc_idx-1),...
            1, 1, 1, []);
        iauc = sum( bsxfun(@times,...
            0.5*dyn_concentrations(:,:,:,s_idx:iauc_idx-1) + ...
            0.5*dyn_concentrations(:,:,:,s_idx+1:iauc_idx) , ...
            time_intervals), 4);
end

%Mask out any non-enhancing voxels
if exist('enhancing_mask', 'var')
    iauc(~enhancing_mask) = 0;
end
    
    



