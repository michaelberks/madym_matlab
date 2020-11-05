function [C_t, Cp_t] = ...
    materne_model(k1a, k1v, k2, aoffset, Ca_t, Cv_t, dyn_times, Hct, voffset)
%MATERNE_MODEL compute predicted tissue concentration time-series given
%a set of input PK parameters assuming the Materne DCE model
%   [C_t, Cp_t] = gadoxetate_model(F, ve, vp)
%
% Inputs:
%      k1a - arterial flow-rate constant
%
%      k1v - (modified) hepatic-portal vein flow-rate constant
%
%      k2 - outflow rate constant
%
% 
%      Ca_t, Cv_t - arterial and venous input functions
%
%      dyn_times - time of each time point
%       
%      Hct - haematocrit constant (assumed 0.42)
%
%      aoffset, voffset - offset times of arrival for conccentraion for
%      Ca_t and Cv_t
%
% Outputs:
%      C_t - concentration at time t
%
%
% Example:
%
% Notes:
%
%              |
%     k1a*Ca_t |
%              v
% k1v*Cv_t |--------|  k2   
%  ------> |  Vecs  |------>
%          |________|       

%Fig1. Dual input single compartment uptake and efflux model
%
%
% Concentration model equation
%   C_t = exp(-k2*t) * Cp_t
%
% Where
%   Cp_t = (k1a.Ca_t + k1v.Cv_t) / (1 - Hct)
%
% 
% See paper: Materne R, Van Beers BE, Smith AM, et al. Non-invasive quantification
% of liver perfusion with dynamic computed tomography and
% a dual-input one-compartmental model. Clin Sci (Lond) 2000;99:517–525
%
%
% See also:
%
% Created: 22-May-2018
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester

if ~exist('Hct', 'var') || isempty(Hct)
    Hct = 0.42;
end
if ~exist('aoffset', 'var') || isempty(aoffset)
    aoffset = zeros(size(k1a));
end
if ~exist('voffset', 'var') || isempty(voffset)
    voffset = zeros(size(k1a));
end

%Make sure model params are all m x 1 
k1a = double(k1a(:));
k1v = double(k1v(:));
k2 = double(k2(:));
aoffset = double(aoffset(:));
voffset = double(voffset(:)) + aoffset;

n_t = length(dyn_times);
n_vox = length(k1a);

%Finally convolve the exchange time series with the combined input function
C_t = zeros(n_vox, n_t);
Cp_t = zeros(n_vox, n_t);
for i_t = 2:n_t
    
    %Get current time, and time change
    t1 = dyn_times(i_t);
    delta_t = t1 - dyn_times(i_t-1);
    
    %Compute (offset) combined arterial and vascular input for this time
    Ca_ti = interp1(dyn_times, Ca_t, t1 - aoffset, 'linear', 'extrap');
    Cv_ti = interp1(dyn_times, Cv_t, t1 - voffset,...
        'linear', 'extrap');   
    Cp_t(:,i_t) = (k1a.*Ca_ti + k1v.*Cv_ti) / (1 - Hct);
    
    %Update the exponentials for the transfer terms in the two compartments
    e_delta = exp(-delta_t * k2);
    
    %Combine the two compartments with the rate constant to get the final
    %concentration at this time point
    A = delta_t * 0.5 * (Cp_t(:,i_t) + Cp_t(:,i_t-1).*e_delta);
    C_t(:,i_t) = e_delta .* C_t(:,i_t-1) + A;
    
end







