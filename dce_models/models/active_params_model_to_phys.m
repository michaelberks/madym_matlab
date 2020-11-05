function [F_p, v_ecs, k_i, k_ef] = ...
    active_params_model_to_phys(F_pos, F_neg, K_pos, K_neg, using_Fp, warn_mode)
%ACTIVE_PARAMS_MODEL_TO_PHYS starting with the derived parameters fitted in
% the IRF-3 model, convert to the physiological parameters F_p, v_ecs, k_i
% and k_ef
%model given input physiological parameters
%   [F_p, v_ecs, k_i, k_ef] = active_params_model_to_phys(K_pos, K_neg, F_pos, F_neg)
%
% Inputs:
%
%      K_pos, K_neg - exponents in 2CXM model IRF
%
%      F_pos, F_neg - scalars in 2CXM model IRF
% Outputs:
%      F_p - flow plasma rate
%
%      v_ecs - extra-cellular space (v_i = 1 - v_ecs)
%
%      k_i - active-uptake rate
%
%      k_ef - efflux rate
%
% Example:
%
% Concentration model equation
%   Cl_t = F_p.(E_i.exp(-t/Ti) / (1 - T_e/T_i) + (1 - E_i/(1 - T_e / T_i)).exp(-t/Te)) * Cp_t
%
% Where
%   Cp_t = (f_a.Ca_t + f_v.Cv_t) / (1 - Hct)
%
%   F_p - flow plasma rate
%   T_e = v_ecs / (F_p + k_i) - extracellular mean transit time
%   T_i = vi / kef - intracellular mean transit time
%   E_i = ki / (Fp + ki) - the hepatic uptake fraction
%   f_a - the arterial fraction
%   f_v = 1 - fa - estimate of hepatic portal venous fraction
%   v_i = 1 - v_ecs - estimate of intracellular volume
% 
% See paper: Invest Radiol. 2017 Feb;52(2):111-119. doi: 10.1097/RLI.0000000000000316.
%   "Quantitative Assessment of Liver Function Using Gadoxetate-Enhanced Magnetic Resonance Imaging: 
%   Monitoring Transporter-Mediated Processes in Healthy Volunteers"
%   Georgiou L1, Penny J, Nicholls G, Woodhouse N, Blé FX, Hubbard Cristinacce PL, Naish JH.
%
%
% See also: GADOXETATE_MODEL, ACTIVE_PARAMS_PHYS_TO_MODEL
%
% Created: 01-Feb-2019
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester
if ~exist('using_Fp', 'var') || isempty('using_Fp')
    using_Fp = false;
end

if ~exist('warn_mode', 'var') || isempty('warn_mode')
    warn_mode = 'warn';
end

%First get F_p from F_pos and F_neg
if ~using_Fp
    F_p = F_pos + F_neg;
    E_pos = F_pos ./ F_p;
else
    F_p = F_pos;
    E_pos = F_neg;
end

%Derivation is only valid for K_pos < K_neg. If not, the swapping
%F_pos, K_pos for F_neg, K_neg will generate valid active parameters (and
%an indentical concentration time series due to symmetry of the
%bi-exponential). User defines whether swap with warning, quietly or force
%an error if invalid voxels found
swap_idx = K_pos > K_neg;
if any(swap_idx(:))
    switch warn_mode
        case 'warn'
            warning('K_pos > K_neg for %d of %d voxels. Switching these voxels.\n',...
                sum(swap_idx(:)), numel(swap_idx));
        case 'error'
            error(['K_pos > K_neg for %d of %d voxels. ' ...
                'Run with warn_mode = ''quiet'' or ''warn to switch these voxels.\n'],...
                sum(swap_idx(:)), numel(swap_idx));
        case 'quiet'
            %do nothing
        otherwise
            error('Warn mode %s not recognised. Must be ''warn'', ''quiet'' or ''error''',...
                warn_mode);
    end
       
    if ~using_Fp
        %F_p doesn't change it is the sum of F_pos and F_neg
        %E_pos needs to remade from F_neg for the swapped indices
        E_pos(swap_idx) = F_neg(swap_idx) ./ F_p(swap_idx);
    else
        %F_p doesn't change, E_pos needs negating
        E_pos(swap_idx) = 1 - E_pos(swap_idx);       
    end
    
    %K_pos and K_neg are just a straight swap
    K_pos_swap = K_pos(swap_idx);
    K_pos(swap_idx) = K_neg(swap_idx);   
    K_neg(swap_idx) = K_pos_swap;
    
    
end

%Now derive Te, Ti and Ei
Te = 1 ./ K_neg;
Ti = 1 ./ K_pos;
Ei = E_pos .* (1 - Te ./ Ti);

%Can solve for k_i in terms of F_p and Ei
k_i = Ei .* F_p ./ (1 - Ei);

%Solve for v_ecs in terms of Te, F_p and K-i
v_ecs = Te .* (F_p + k_i);

%Finally solve for k_ef in terms of v_ecs and Ti
k_ef = (1 - v_ecs) ./ Ti;











