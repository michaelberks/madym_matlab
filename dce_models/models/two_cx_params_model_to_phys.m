function [F_p, PS, v_e, v_p] = ...
    two_cx_params_model_to_phys(F_pos, F_neg, K_pos, K_neg, using_Fp)
%TWO_CX_PARAMS_MODEL_TO_PHYS starting with the derived parameters fitted in
% the 2CXM model, convert to the physiological parameters F_p, PS, ve and vep
%model given input physiological parameters
%   [F_p, PS, v_e, v_p] = two_cx_params_model_to_phys(K_pos, K_neg, F_pos, F_neg)
%
% Inputs:
%
%      K_pos, K_neg - exponents in 2CXM model IRF
%
%      F_pos, F_neg - scalars in 2CXM model IRF
% Outputs:
%      F_p - flow plasma rate
%
%      PS - extraction flow
%
%      v_e - extravascular, extracellular volume
%
%      v_p - plasma volume
%
% Example:
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
% See also: TWO_CX_MODEL, TWO_CX_PARAMS_PHYS_TO_MODEL
%
% Created: 01-Feb-2019
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester
if ~exist('using_Fp', 'var')
    using_Fp = false;
end

%We derive the params based on 2009 Sourbron paper
if ~using_Fp
    F_p = F_pos + F_neg;
    E_neg = F_neg ./ F_p;
else
    F_p = F_pos;
    E_neg = (1 - F_neg);
end

T_B = 1 ./ (K_pos - E_neg .* (K_pos - K_neg));
T_E = 1 ./ (T_B .* K_pos .* K_neg);
T_P_inv = K_pos + K_neg - 1 ./ T_E;

v_p = F_p .* T_B;
PS = F_p .* (T_B .* T_P_inv - 1);
v_e = PS .* T_E; 

apply_tm = K_pos==0 & F_pos==0;
if any(apply_tm(:))
    PS(apply_tm) = NaN;
    v_p(apply_tm) = 0;
    v_e(apply_tm) = F_p(apply_tm) ./ K_neg(apply_tm);
end






