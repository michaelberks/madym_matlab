function [C_t] = ...
    two_cx_model(F_p, PS, v_e, v_p,  aoffset, Ca_t, dyn_times)
%GADOXETATE_MODEL compute predicted tissue concentration time-series given
%a set of input PK parameters assuming a Gadoxetate specific DCE model
%   [C_t, Cp_t] = gadoxetate_model(F, ve, vp)
%
% Inputs:
%      F_p - flow plasma rate
%
%      PS - extraction flow
%
%      v_e - extravascular, extracellular volume
%
%      v_p - plasma volume 
%
%      Ca_t - arterial input functions
%
%      dyn_times - time of each time point
%
%      aoffset - offset times of arrival for conccentraion for Ca_t
%
% Outputs:
%      C_t - concentration at time t
%
%
% Example:
%
% Notes:
%
%
% Concentration model equation
%    C(t) = F_p*[(1 - E_pos).exp(-t*T_neg) + E_pos.exp(-t*E_pos)] ** Ca(t)
%
% Where
%   T_pos = K_sum + K_root;
%   T_neg = K_sum - K_root;
% 
%   E_pos = (T_neg - Kb) ./ (T_neg + T_pos);
%
% Derived from
%
%   Kp = (F_p + PS) ./ v_p;
%   Ke = PS ./ v_ecs;
%   Kb = F_p ./ v_p;
%
% Where
%
%   F_p - flow plasma rate
%   PS = extraction flow
%   v_e - extra cellular extra vascular volume
%   v_p - plasma vlume
% 
% See paper: Phys Med Bio. 2010;55:6431-6643
%   "Error estimation for perfusion parameters obtained using the 
%   two-compartment exchange model in dynamic contrast-enhanced MRI: a simulation study"
%   R Luypaert, S Sourbron, S Makkat and J de Mey.
%
%
% See also:
%
% Created: 01-Feb-2019
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester
if ~exist('aoffset', 'var') || isempty(aoffset)
    aoffset = zeros(size(F_p));
end

K_max = 1e9;

%Get number of times and number of voxels, noting we allow model params to
%be scalar (so need to check lengths of all inputs) - this makes it easier
%to fix all params but one, and test a range of inputs on the other to
%observe how that parameter alters the model
n_t = length(dyn_times);

F_p = double(F_p(:));
PS = double(PS(:));
v_e = double(v_e(:));
v_p = double(v_p(:));
aoffset = double(aoffset(:));

n_vox = max([...
    length(F_p)
    length(PS)
    length(v_e)
    length(v_p)
    length(aoffset)]);

if length(F_p) ~= n_vox
    if length(F_p) == 1
        F_p = ones(n_vox,1)*F_p;
    else
        error('Incorrect size for F_p, should be 1 x %d', n_vox);
    end
end
if length(PS) ~= n_vox
    if length(PS) == 1
        PS = ones(n_vox,1)*PS;
    else
        error('Incorrect size for PS, should be 1 x %d', n_vox);
    end
end
if length(v_e) ~= n_vox
    if length(v_e) == 1
        v_e = ones(n_vox,1)*v_e;
    else
        error('Incorrect size for v_e, should be 1 x %d', n_vox);
    end
end
if length(v_p) ~= n_vox
    if length(v_p) == 1
        v_p = ones(n_vox,1)*v_p;
    else
        error('Incorrect size for v_p, should be 1 x %d', n_vox);
    end
end
if length(aoffset) ~= n_vox
    if length(aoffset) == 1
        aoffset = ones(n_vox,1)*aoffset;
    else
        error('Incorrect size for aoffset, should be 1 x %d', n_vox);
    end
end

%We derive the params in a standalone function now, this takes care of
%checks on FP, PS to choose the best form of derived parameters
[F_pos, F_neg, K_pos, K_neg] = ...
    two_cx_params_phys_to_model(F_p, PS, v_e, v_p);

%Irf is of form: I(t) = F_pos.exp(-tK_pos) + F_neg.exp(-tK_neg)
%C(t) = I(t) ** Ca(t)
C_t = zeros(n_vox,n_t); 
Ca_t0 = interp1(dyn_times, Ca_t, dyn_times(1) - aoffset, 'linear', 'extrap');
Ft_pos = 0;
Ft_neg = 0;
for i_t = 2:n_t
    
    %Get current time, and time change
    t1 = dyn_times(i_t);
    delta_t = t1 - dyn_times(i_t-1);
    
    %Compute (offset) combined arterial and vascular input for this time
    Ca_ti = interp1(dyn_times, Ca_t, t1 - aoffset, 'linear', 'extrap');
    
    %Update the exponentials for the transfer terms in the two compartments        
    et_pos = exp(-delta_t .* K_pos);
    et_neg = exp(-delta_t .* K_neg);
        
    %Use iterative trick to update the convolutions of transfers with the
    %input function. This only works when the exponent is finite, otherwise
    %the exponential is zero, and the iterative solution is not valid. For
    %these voxels, set A_pos/neg to zero
    A_pos = delta_t * 0.5 * (Ca_ti + Ca_t0.*et_pos);
    A_pos(K_pos > K_max) = 0;    
    
    A_neg = delta_t * 0.5 * (Ca_ti + Ca_t0.*et_neg);
    A_neg(K_neg > K_max) = 0;
    
    Ft_pos = Ft_pos.*et_pos + A_pos;
    Ft_neg = Ft_neg.*et_neg + A_neg;
    
    %Combine the two compartments with the rate constant to get the final
    %concentration at this time point
    C_t(:,i_t) = F_pos.*Ft_pos + F_neg.*Ft_neg;
    Ca_t0 = Ca_ti;
end






