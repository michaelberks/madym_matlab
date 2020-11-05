function [C_t, Cp_t] = ...
    gadoxetate_model(F_p, v_ecs, k_i, k_ef, f_a, aoffset, Ca_t, Cv_t, dyn_times, Hct, voffset)
%GADOXETATE_MODEL compute predicted tissue concentration time-series given
%a set of input PK parameters assuming a Gadoxetate specific DCE model
%   [C_t, Cp_t] = gadoxetate_model(F, ve, vp)
%
% Inputs:
%      F_p - flow plasma rate
%
%      v_ecs - extravascular, extracellular space
%
%      k_i - transport constant to liver tissue
%
%      k_ef - efflux to bile duct transfer constant 
%
%      f_a - the arterial fraction
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
%      Cp_t - concentration input function
%
%
% Example:
%
% Notes:
%
%              |
%        fa.Fp |
%              v
%   fvFp   |--------|  ki   |-------|   kef
%  ------> |  Vecs  |------>|   vi  |------->  
%          |________|       |_______|
%              |
%              |
%           Fp |
%              v
%Fig1. Dual input 2-compartmental uptake and efflux model
%
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
    aoffset = zeros(size(F_p));
end
if ~exist('voffset', 'var') || isempty(voffset)
    voffset = zeros(size(F_p));
end

F_p = double(F_p(:));
v_ecs = double(v_ecs(:));
k_i = double(k_i(:));
k_ef = double(k_ef(:));
f_a = double(f_a(:));
aoffset = double(aoffset(:));
voffset = double(voffset(:));

%Get number of times and number of voxels, noting we allow model params to
%be scalar (so need to check lengths of all inputs) - this makes it easier
%to fix all params but one, and test a range of inputs on the other to
%observe how that parameter alters the model
n_t = length(dyn_times);
n_vox = max([...
    length(F_p)
    length(v_ecs)
    length(k_i)
    length(k_ef)
    length(f_a)
    length(aoffset)
    length(voffset)]);

%Compute derived parameters from input parameters
T_e = v_ecs ./ (F_p + k_i); % extracellular mean transit time
v_i = 1 - v_ecs; % - etsimate of intracellular volume
T_i = v_i ./ k_ef; % intracellular mean transit time
E_i = k_i ./ (F_p + k_i); % the hepatic uptake fraction

f_v = 1 - f_a; %   estimate of hepatic portal venous fraction
k1a = F_p.*f_a;
k1v = F_p.*f_v;

%This can also be precomputed
ETie = E_i ./ (1 - T_e./T_i);

% Let's rewrite the convolution sum, using the exponential "trick" so 
% we can compute everything in one forward loop
Fi_t = 0;
Fe_t = 0;

C_t = zeros(n_vox,n_t); 
Cp_t = zeros(n_vox,n_t);
    
for i_t = 2:n_t
    
    %Get current time, and time change
    t1 = dyn_times(i_t);
    delta_t = t1 - dyn_times(i_t-1);
    
    %Compute (offset) combined arterial and vascular input for this time
    Ca_ti = interp1(dyn_times, Ca_t, t1 - aoffset, 'linear', 'extrap');
    Cv_ti = interp1(dyn_times, Cv_t, t1 - (aoffset + voffset),...
        'linear', 'extrap');   
    Cp_t(:,i_t) = (k1a.*Ca_ti + k1v.*Cv_ti) / (1 - Hct);
    
    %Update the exponentials for the transfer terms in the two compartments
    et_i = exp(-delta_t ./ T_i);
    et_e = exp(-delta_t ./ T_e);
    
    %Use iterative trick to update the convolutions of transfers with the
    %input function
    A_i = delta_t * 0.5 * (Cp_t(:,i_t) + Cp_t(:,i_t-1) .* et_i);
    A_e = delta_t * 0.5 * (Cp_t(:,i_t) + Cp_t(:,i_t-1) .* et_e);
    Fi_t = Fi_t.*et_i + A_i;
    Fe_t = Fe_t.*et_e + A_e;
    
    %Combine the two compartments with the rate constant to get the final
    %concentration at this time point
    C_t(:,i_t) = (ETie.*Fi_t + (1-ETie).*Fe_t);
    
end

% % The old method, we now use the convolution trick, to compute everything in
% % one forward loop
% t = dyn_times;
% ex_i = exp(-t ./ T_i);
% ex_e = exp(-t ./ T_e);
% %Compute the exchange term for the concentration time series
% Cx_t =  F_p.*(ETie .* ex_i + (1 - ETie).*ex_e);
% 
% %Finally convolve the exchange time series with the combined input function
% Cl_t = zeros(n_t, n_vox);
% delta_t = [0; diff(t)];
% for i_vox = 1:n_vox
%     
%     %If needed offset the input functions
%     if aoffset(i_vox) > 0
%         Ca_ti = interp1(t, Ca_t,...
%             t - aoffset(i_vox), 'linear', 0);
%     else
%         Ca_ti = Ca_t;
%     end
%     if aoffset(i_vox)>0 || voffset(i_vox)>0
%         Cv_ti = ...
%             interp1(t, Cv_t,...
%             t - (aoffset(i_vox)+voffset(i_vox)), 'linear', 0);
%     else
%         Cv_ti = Cv_t;
%     end
%     
%     %Make combined input function
%     Cp_t = (f_a.*Ca_ti + f_v.*Cv_ti) / (1 - Hct);
% 
%     Cl_i = conv(Cx_t(:,i_vox), Cp_t(:,i_vox));
%     Cl_t(:,i_vox) = Cl_i(1:n_t) .* delta_t;
% end

%%
%Below is an implementation that loops over time rather voxels, and
%manually computes the convolution - this is how the model is
%implemented in the C++ version of madym. Should produce identical results,
%the implementation above is faster in Matlab, but below is useful for
%debugging vis-a-vis the C++ version
%     Cp_t = zeros(n_t, n_vox);
%     Cx_t = zeros(n_t, n_vox);
%     %Cl_t = zeros(n_t, n_vox);
%     for i_t = 1:n_t
%         ex_i = exp(-t(i_t) ./ T_i);
%         ex_e = exp(-t(i_t) ./ T_e);
%         ETie = E_i ./ (1 - T_e./T_i);
% 
%         %Make combined input function
%         Cp_t(i_t,:) = (f_a*Ca_t(i_t) + f_v*Cv_t(i_t)) / (1 - Hct);
% 
%         %Compute the exchange term for the concentration time series
%         Cx_t(i_t,:) =  F_p.*(ETie .* ex_i + (1 - ETie).*ex_e);
% 
%         %Finally convolve the exchange time series with the combined input function
%         Cl_t_sum = zeros(1, n_vox);
%         for j_t = 1:i_t
%             Cl_t_sum = Cl_t_sum + ...
%                 Cx_t(j_t,:) .*  Cp_t(i_t-j_t+1,:);
%             k_t = k_t - 1;
%         end
%         Cl_t(i_t,:) = Cl_t_sum;
% 
%     end





