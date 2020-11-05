function [F_pos, F_neg, K_pos, K_neg] = ...
    active_params_phys_to_model(F_p, v_ecs, k_i, k_ef)
%ACTIVE_PARAMS_PHYS_TO_MODEL compute the derived parameters for the active
%-uptake model given input physiological parameters
%   [K_pos, K_neg, F_pos, F_neg] = active_params_phys_to_model(F_p, v_e, k_i, k_ef)
%
% Inputs:
%      F_p - flow plasma rate
%
%      v_ecs - extra-cellular space (v_i = 1 - v_ecs)
%
%      k_i - active-uptake rate
%
%      k_ef - efflux rate
%
% Outputs:
%
%      K_pos, K_neg - exponents in model IRF
%
%      F_pos, F_neg - scalars in model IRF
%
%
% Example:
%
% Notes:
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
% See also: GADOXETATE_MODEL, ACTIVE_PARAMS_MODEL_TO_PHYS
%
% Created: 01-Feb-2019
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester

%Compute derived parameters from input parameters
T_e = v_ecs ./ (F_p + k_i); % extracellular mean transit time
v_i = 1 - v_ecs; % - etsimate of intracellular volume
T_i = v_i ./ k_ef; % intracellular mean transit time
E_i = k_i ./ (F_p + k_i); % the hepatic uptake fraction

%This can also be precomputed
E_pos = E_i ./ (1 - T_e./T_i);

K_neg = 1 ./ T_e;
F_neg = F_p.*(1 - E_pos);

K_pos = 1 ./ T_i;
F_pos = F_p.*E_pos;







