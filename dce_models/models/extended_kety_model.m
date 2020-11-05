function [C_t] = extended_kety_model(dyn_times, Ca_t, Ktrans, Vp, Ve, offset)
%EXTENDED_KETY_MODEL Obsolete - use extended_tofts_model instead
%   [model_signals] = extended_kety_model(dyn_times, aif, Ktrans, Vp, Ve)
%
% Inputs:
%      dyn_times - *Insert description of input variable here*
%
%      Ca_t - Arterial input function
%
%      Ktrans - *Insert description of input variable here*
%
%      Vp - *Insert description of input variable here*
%
%      Ve - *Insert description of input variable here*
%
%
% Outputs:
%      model_signals - *Insert description of input variable here*
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

n_t = size(Ca_t, 2); %n
num_voxels = max([...
    length(Ktrans)...
    length(Ve)...
    length(Vp)...
    length(offset)]); %m

%Make sure model params are all m x 1 
Ktrans = double(Ktrans(:));
Vp = double(Vp(:));
Ve = double(Ve(:));
offset = double(offset(:));

%Precompute exponent, k_ep
k_ep = Ktrans ./ Ve;

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
    C_t(:,i_t) = Vp * Ca_ti + Ktrans .* e_i;
    
    Ca_t0 = Ca_ti;
    
end

