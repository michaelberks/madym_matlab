function [accumulation_map] = make_DIBEM_accumlation_map(...
    C_t, dyn_t, injection_img, period_1_times, period_2_times)
%MAKE_DIBEM_ACCUMLATION_MAP make map of voxels where a DIBEM model
%accumlates CA
%   [accumulation_map] = make_DIBEM_accumlation_map(C_m, dyn_t, injection_img)
%
% Inputs:
%      C_t - 4-D array of contrast time-series, with time on 4th dimension 
%      
%      dyn_t - time of each point in time series
%
%      injection_img - index of image at which CA was injected
%
%      period_1_times - 2 element vector, start and end times (in minutes)
%      of first period. Default [1 2]
%
%      period_2_times - 2 element vector, start and end times (in minutes)
%      of second period. Default [4 5]
%
%
% Outputs:
%      accumulation_map - 3D logical array, matching size of C_m, 1
%      wherever modelled contrast in period 2 is greater than period 1
%
%
% Example:
%
% Notes: accumulation is defined
%
% See also:
%
% Created: 30-Jun-2020
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 
if ~exist('period_1_times', 'var') || isempty(period_1_times)
    period_1_times = [1 2];
end
if ~exist('period_2_times', 'var') || isempty(period_2_times)
    period_2_times = [4 5];
end

%Set dynamic times so the injection image is at time zero
dyn_t = dyn_t - dyn_t(injection_img);
t_period1 = (dyn_t > period_1_times(1)) & (dyn_t < period_1_times(2));
t_period2 = (dyn_t > period_2_times(1)) & (dyn_t < period_2_times(2));
C_period1 = mean(C_t(:,:,:,t_period1),4);
C_period2 = mean(C_t(:,:,:,t_period2),4);
accumulation_map = C_period2 > C_period1;