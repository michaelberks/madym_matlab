function [mse, mae, time_series_error] = compute_time_series_error(time_series, start_t, end_t)
%COMPUTE_TIME_SERIES_ERROR return mean-squared and mean absolute error for
%time-series data, defined for each time-point as the difference between it
%and the mid-point of its two temporal neighbours
%   [mean_error, time_series_error] = compute_time_series_error(time_series, start_t, end_t)
%
% Inputs:
%      time_series - N-d array with time on last dimension
%
%      start_t - firts time-point to include in total error
%
%      end_t - last timepoint to inlcude
%
%
% Outputs:
%      mse - mean-squared error
%
%      mae - mean absolut error
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 04-May-2018
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 
if ~exist('start_t', 'var') || isempty(start_t)
    start_t = 2;
elseif start_t < 2
    start_t = 2;
    warning('start t must be > 1. setting start_t = 2');
end

n_dims = ndims(time_series);
t_dim = n_dims;

n_t = size(time_series, t_dim);
if ~exist('end_t', 'var') || isempty(end_t)
    end_t = n_t-1;
elseif end_t >= n_t
    end_t = n_t-1;
    warning('end_t must be < n_t (%d), setting end_t = %d', n_t, end_t);
    
end

n_times = end_t - start_t +1;
switch n_dims
    
    case 2
        t0 = time_series(:,start_t:end_t);
        t_m1 = time_series(:,(start_t-1):(end_t-1));
        t_p1 = time_series(:,(start_t+1):(end_t+1));
    case 3
        t0 = time_series(:,:,start_t:end_t);
        t_m1 = time_series(:,:,(start_t-1):(end_t-1));
        t_p1 = time_series(:,:,(start_t+1):(end_t+1));
    case 4
        t0 = time_series(:,:,:,start_t:end_t);
        t_m1 = time_series(:,:,:,(start_t-1):(end_t-1));
        t_p1 = time_series(:,:,:,(start_t+1):(end_t+1));
end
        

t_avg = (t_m1 + t_p1) / 2;

time_series_error = t_avg - t0;

mse = sum(time_series_error.^2, t_dim) / n_times;
mae = sum(abs(time_series_error), t_dim) / n_times;





