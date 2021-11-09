function [Cp_t] = separable_compartment_model(v_p, T_p, F_t, T_t, tau_a, Ca_t, dyn_times)
%SEPARABLE_COMPARTMENT_MODEL Compute concentration time-series from model
%parameters for the separable compartment model
%   [C_t] = separable_compartment_model(v_p, T_p, F_t, T_t, tau, Ca_t, dyn_times)
%
% Inputs:
%
%      Ktrans, Ve, VP, offset - PK parameters for the tofts model, can either be
%      scalar or arrays of inputs. If more than one non-scalar, they must
%      have the same number of elements
%
%      Ca_t - Arterial input function
%
%      dyn_times - 1d-array of times associated with each time-point in
%      dynamic series
%
%
% Outputs:
%      C_t - Model concentration time series
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

n_t = size(Ca_t, 1); %n
num_voxels = max([...
    numel(v_p)...
    numel(T_p)...
    numel(F_t)...
    numel(T_t)...
    numel(tau_a)]); %m

%Make sure model params are all m x 1
v_p = double(v_p(:));
T_p = double(T_p(:));
F_t = double(F_t(:));
T_t = double(T_t(:));
tau_a = double(tau_a(:));

%Precompute inverse of T_t and T_p
T_pi = 1 ./ T_p;
T_ti = 1 ./ T_t;

%%
Cp_t = zeros(num_voxels,n_t);
C_t = zeros(num_voxels,n_t);

Ca_t0 = Ca_t(1);
Cp_t0 = Cp_t(1);

%Set the initial integral sums to zero, these will be aggregated as we move
%through the time series, adding new terms at each time using the trapezium
%rule
ep_i = 0;
et_i = 0;
for i_t = 2:n_t
    
    %Get current time, and time change
    t1 = dyn_times(i_t);
    delta_t = t1 - dyn_times(i_t-1);
    
    %Compute (offset) arterial input for this time
    Ca_ti = interp1(dyn_times, Ca_t, t1 - tau_a, 'linear', 'extrap');
    
    %Update the integral for the Cp(t) convolution
    ep_delta = exp(-delta_t * T_pi);
    Ap = delta_t * 0.5 * (Ca_ti + Ca_t0.*ep_delta);
    ep_i = ep_i .* ep_delta + Ap;
    Cp_t(:,i_t) = T_pi .* ep_i;
    
    %Update the integral for the C(t) convolution
    et_delta = exp(-delta_t * T_ti);
    At = delta_t * 0.5 * (Cp_ti + Cp_t0.*et_delta);
    et_i = et_i .* et_delta + At;
    
    C_t(:,i_t) = v_p * Ca_ti + F_t .* et_i;
    
    Ca_t0 = Ca_ti;
    Cp_t0 = Cp_ti;
    
end

