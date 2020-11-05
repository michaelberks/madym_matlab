function [pif, aif] = compute_PIF(aif, irf, dyn_times, prebolus_im, offset)
%COMPUTE_PIF computes hepatic portal vein input function (PIF) as
%convultion of an arterial input function with an impulse response function
%
%   [pif] = compute_PIF(aif, irf, dyn_times)
%
% Inputs:
%      aif - arterial input function, if empty, will use population model
%
%      irf - impulse repsonse function, if empty, will use population model
%
%      dyn_times - only required if generating aif or irf
%
%
% Outputs:
%      pif - hepatic portal vein input function
%
%
% Example:
%
% Notes:
%
% See also: DCE?MRI model selection for investigating disruption of microvascular function in livers with metastatic disease
% A Banerji, J Naish, Y Watson, G Jayson, G Buonaccorsi, G Parker
% 10 October 2011, https://doi.org/10.1002/jmri.22692
%
% Created: 22-May-2018
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 

%If the AIF is empty, generate a population AIF according the Geoff's model
if isempty(aif)
    aif = population_aif(dyn_times, prebolus_im);
end

if ~exist('offset', 'var') || isempty(offset)
    offset = 0;
end

%If the IRF is empty, generate a population IRF according to Anita's model
n_t = length(dyn_times);
if isempty(irf)
    
    irf = zeros(n_t,1);
    for i_t = 1:n_t
        t = dyn_times(i_t) - offset;        
        if t < 0.08% - do nothing, irf(i_t) = 0
        elseif t < 0.17
            irf(i_t) = 24.16*t - 2.01;
        else
            irf(i_t) = 2.83*exp(-10.80*t) + 2.12*exp(-1.82*t);
        end            
    end
    irf = irf / sum(irf);
end

%Convolve the AIF with the IRF to generate the PIF
pif = conv(aif, irf);
pif(n_t+1:end) = [];

% % literal convolution operation
% pif2 = zeros(n_t,1);
% for i_t = 1:n_t
%     k_t = i_t;
%     pif_sum = 0;
%     for j_t = 1:i_t
%         pif_sum = pif_sum + aif(j_t)*irf(k_t);
%     	k_t = k_t - 1;
%     end
%     pif2(i_t) = pif_sum;
% end
% pif - pif2;
