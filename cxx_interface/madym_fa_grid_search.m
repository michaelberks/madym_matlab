function [local_fit, grid_search, grid_fit] = madym_fa_grid_search(C_t, ...
    madym_args)
%MADYM_FA_GRID_SEARCH *Insert a one line summary here*
%   [local_fit, grid_fit] = madym_fa_grid_search(C_t, varargin)
%
% Inputs:
%      C_t - [n_pts, n_times] Array of concentration time-series
%
%      madym_args - input options to configure madym see RUN_MADYM_LITE
%
%
% Outputs:
%      local_fit - structure containing madym output for locally
%      constrained fa fit. Contains fields
%           model_params [n_pts, n_params]
%           model_fit [n_pts, 1]
%           C_t_m [n_pts, n_times]
%           
%
%      grid_search - structure containing madym output for each fixed fa
%      value in the grid search
%
%      grid_fit - structure containing madym output for the best fitting fa
%      for each point
%
%
% Example:
%
% Notes:
%
% See also: RUN_MADYM_LITE
%
% Created: 18-Jun-2020
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 

%TODO - get model size from selected options
n_params = 7;

%Set default fixed fa (TODO - put in options)
grid_fa_vals = linspace(0,1,11);
n_fas = length(grid_fa_vals);

%Set up containers for model parameters and fit from grid search of fixed
%fa, and for the best fit from the grid
[n_pts, n_t] = size(C_t);
grid_search.model_params = zeros(n_pts, n_params, n_fas);
grid_search.model_fit = zeros(n_pts, n_fas);
grid_search.C_t_m = zeros(n_pts, n_t, n_fas);
                
grid_fit.model_params = zeros(n_pts, n_params);
grid_fit.model_fit = zeros(n_pts, 1);
grid_fit.C_t_m = zeros(n_pts, n_t);

madym_args.fixed_params = [5 7];

%Loop over fixed fa values                    
for i_fa = 1:n_fas
    fprintf('Fitting with fixed fa = %2.1f\n', ...
        grid_fa_vals(i_fa));

    madym_args.fixed_values = [grid_fa_vals(i_fa) 0];

    %Call run madym lite for fixed fa
    [...
        grid_search.model_params(:,:,i_fa), ...
        grid_search.model_fit(:,i_fa) ,~,~,...
        grid_search.C_t_m(:,:,i_fa)] =...
        run_madym_lite('DIBEM', C_t, madym_args);

    % Keep track of best fit          
    if i_fa == 1
       grid_fit.model_fit = grid_search.model_fit(:,i_fa);
       grid_fit.model_params = grid_search.model_params(:,:,i_fa);
       grid_fit.C_t_m = grid_search.C_t_m(:,:,i_fa);
    else
       swap_idx = grid_search.model_fit(:,i_fa) < grid_fit.model_fit;
       grid_fit.model_fit(swap_idx) = ...
           grid_search.model_fit(swap_idx,i_fa);
       grid_fit.model_params(swap_idx,:) = ...
           grid_search.model_params(swap_idx,:,i_fa);
       grid_fit.C_t_m(swap_idx,:) = ...
           grid_search.C_t_m(swap_idx,:,i_fa);

    end
end

%Now perform locally constrained fit - TODO where to source these initial
%params from?
init_params = repmat([0.0,0.6,0.0,0.5,1.0,0.025,0.0], n_pts, 1);
init_params(:,5) = grid_fit.model_params(:,5);

madym_args.relative_limit_values = 0.1;
madym_args.relative_limit_params = 5;
madym_args.init_params = init_params;
madym_args.fixed_values = [];
madym_args.fixed_params = 7;
[...
    local_fit.model_params,...
    local_fit.model_fit,~,~,...
    local_fit.C_t_m] =...
    run_madym_lite('DIBEM', C_t, madym_args);