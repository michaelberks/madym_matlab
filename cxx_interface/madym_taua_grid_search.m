function [local_fit, grid_search, grid_fit] = madym_taua_grid_search(C_t, ...
    madym_args)
%MADYM_TAUA_GRID_SEARCH *Insert a one line summary here*
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
grid_tau_vals = linspace(0,0.5,11);
n_taus = length(grid_tau_vals);

%Set up containers for model parameters and fit from grid search of fixed
%fa, and for the best fit from the grid
[n_pts, n_t] = size(C_t);
grid_search.model_params = zeros(n_pts, n_params, n_taus);
grid_search.model_fit = zeros(n_pts, n_taus);
grid_search.C_t_m = zeros(n_pts, n_t, n_taus);
                
grid_fit.model_params = zeros(n_pts, n_params);
grid_fit.model_fit = zeros(n_pts, 1);
grid_fit.C_t_m = zeros(n_pts, n_t);

madym_args.fixed_params = [3 6 7];

%Loop over fixed fa values                    
for i_tau = 1:n_taus
    fprintf('Fitting with fixed tau_a = %3.2f\n', ...
        grid_tau_vals(i_tau));

    madym_args.fixed_values = [0 grid_tau_vals(i_tau) 0];

    %Call run madym lite for fixed fa
    [...
        grid_search.model_params(:,:,i_tau), ...
        grid_search.model_fit(:,i_tau) ,~,~,...
        grid_search.C_t_m(:,:,i_tau)] =...
        run_madym_lite('DIBEM', C_t, madym_args);

    % Keep track of best fit          
    if i_tau == 1
       grid_fit.model_fit = grid_search.model_fit(:,i_tau);
       grid_fit.model_params = grid_search.model_params(:,:,i_tau);
       grid_fit.C_t_m = grid_search.C_t_m(:,:,i_tau);
    else
       swap_idx = grid_search.model_fit(:,i_tau) < grid_fit.model_fit;
       grid_fit.model_fit(swap_idx) = ...
           grid_search.model_fit(swap_idx,i_tau);
       grid_fit.model_params(swap_idx,:) = ...
           grid_search.model_params(swap_idx,:,i_tau);
       grid_fit.C_t_m(swap_idx,:) = ...
           grid_search.C_t_m(swap_idx,:,i_tau);

    end
end

%Now perform locally constrained fit - TODO where to source these initial
%params from?
init_params = repmat(madym_args.init_params, n_pts, 1);
init_params(:,6) = grid_fit.model_params(:,6);

madym_args.relative_limit_values = 0.05;
madym_args.relative_limit_params = 6;
madym_args.init_params = init_params;
madym_args.fixed_values = [];
madym_args.fixed_params = [3 7];
[...
    local_fit.model_params,...
    local_fit.model_fit,~,~,...
    local_fit.C_t_m] =...
    run_madym_lite('DIBEM', C_t, madym_args);