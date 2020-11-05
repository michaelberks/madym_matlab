function [F_pos, F_neg, K_pos, K_neg] = ...
    two_cx_params_phys_to_model(F_p, PS, v_e, v_p, using_Fp)
%TWO_CX_PARAMS_PHYS_TO_MODEL compute the derived parameters for the 2CXM
%model given input physiological parameters
%   [K_pos, K_neg, F_pos, F_neg] = two_cx_params_phys_to_model(F_p, PS, v_e, v_p)
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
% 2CXM model is bi-exponential, with end concentration computed as
%   C(t) = [ F_pos.exp(-t.K_pos) + F_neg.exp(-t.K_neg) ] ** Ca(t)
%
% Where
%   K_pos = K_sum + K_root;
%   K_neg = K_sum - K_root;
% 
%   E_pos = (T_neg - Kb) ./ (T_neg + T_pos);
%   F_pos = F_p.E_pos
%   F_neg = F_p.(1 - E_pos)
%
% Derived from
%
%   Kp = (F_p + PS) ./ v_p;
%   Ke = PS ./ v_ecs;
%   Kb = F_p ./ v_p;
%   K_sum = 0.5*(Kp + Ke);
%   K_root = 0.5* sqrt( (Kp + Ke).^2 - 4*Ke .*Kb);
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
% See also: TWO_CX_MODEL, TWO_CX_PARAMS_MODEL_TO_PHYS
%
% Created: 01-Feb-2019
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester

%We can derive the params in a couple of ways, which remain stable under
%different conditions of ve, vp, PS and FP

%The first way is as derived in the Sourbron 2011 MRM paper, which is valid
%except when PS = 0 or FP = 0. The second method is as derived in Luypaert 
%paper 2010 paper. It works when PS or FP = 0, but doesn't like ve or vp = 0
if ~exist('using_Fp', 'var')
    using_Fp = false;
end

method1 = PS > 0 & F_p > 0 & (v_e + v_p) > 0;
method2 = ~method1;

%We're assuming all params have been passed in the same size, not doing any
%error checks here
dims_sz = size(F_p);
K_pos = zeros(dims_sz);
K_neg = zeros(dims_sz);
E_pos = zeros(dims_sz);

%% Method 1: Sourbron 2011
%First derive the secondary parameters from the input Pk parameters
E = PS(method1) ./ (PS(method1) + F_p(method1)); %Extraction fraction
e = v_e(method1) ./ (v_p(method1) + v_e(method1)); %Extractcellular fraction

tau = (E - E.*e + e) ./ (2*E);
tau_root = sqrt(1 - 4*(E.*e.*(1-E).*(1-e)) ./ ((E - E.*e + e).^2) );
tau_pos = tau .* (1 + tau_root);
tau_neg = tau .* (1 - tau_root);

K_pos(method1) = F_p(method1) ./ ((v_p(method1) + v_e(method1)).*tau_neg);
K_neg(method1) = F_p(method1) ./ ((v_p(method1) + v_e(method1)).*tau_pos);

E_pos(method1) = (tau_pos - 1) ./ (tau_pos - tau_neg);

%% Method 2
Kp = (F_p(method2) + PS(method2)) ./ v_p(method2);
Ke = PS(method2) ./ v_e(method2);
Kb = F_p(method2) ./ v_p(method2);

K_sum = 0.5*(Kp + Ke);
K_root = 0.5* sqrt( (Kp + Ke).^2 - 4*Ke .*Kb);
K_pos(method2) = K_sum - K_root;
K_neg(method2) = K_sum + K_root;

E_pos(method2) = (K_neg(method2) - Kb) ./ (K_neg(method2) - K_pos(method2)); 
%%
if using_Fp
    F_pos = F_p;
    F_neg = 1 - E_pos;
else
    F_pos = F_p.*E_pos;
    F_neg = F_p.*(1 - E_pos);
end







