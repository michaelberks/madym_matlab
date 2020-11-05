function [C_t] = extended_tofts_model(Ktrans, v_e, v_p, offset, Ca_t, dyn_times)
%EXTENDED_TOFTS_MODEL Compute concentration time-series from model
%parameters for the extended-Tofts model
%   [C_t] = extended_tofts_model(Ktrans, v_e, v_p, offset, Ca_t, dyn_times)
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
    numel(Ktrans)...
    numel(v_e)...
    numel(v_p)...
    numel(offset)]); %m

%Make sure model params are all m x 1
Ktrans = double(Ktrans(:));
v_p = double(v_p(:));
v_e = double(v_e(:));
offset = double(offset(:));

%Precompute exponent, k_ep
k_ep = Ktrans ./ v_e;

%%
C_t = zeros(num_voxels,n_t);
Ca_t0 = interp1(dyn_times, Ca_t, dyn_times(1) - offset, 'linear', 'extrap');
e_i = 0;
for i_t = 2:n_t
    
    %Get current time, and time change
    t1 = dyn_times(i_t);
    delta_t = t1 - dyn_times(i_t-1);
    
    %Compute (offset) combined arterial and vascular input for this time
    Ca_ti = interp1(dyn_times, Ca_t, t1 - offset, 'linear', 'extrap');
    
    %Update the exponentials for the transfer terms in the two compartments
    e_delta = exp(-delta_t * k_ep);
    
    %Combine the two compartments with the rate constant to get the final
    %concentration at this time point
    A = delta_t * 0.5 * (Ca_ti + Ca_t0.*e_delta);

    e_i = e_i .* e_delta + A;
    C_t(:,i_t) = v_p * Ca_ti + Ktrans .* e_i;
    
    Ca_t0 = Ca_ti;
    
end

