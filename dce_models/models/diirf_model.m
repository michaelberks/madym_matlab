function [C_t, Cp_t] = diirf_model(...
    F_pos, F_neg, K_pos, K_neg, ... 
    f_a, aoffset, Ca_t, Cv_t, dyn_times, Hct, voffset, use_fp)
%DIIRF_MODEL compute predicted tissue concentration time-series given
%a set of input functional parameters assuming a generic dual-input, 
% bi-exponential DCE model
%   [C_t, Cp_t] = diirf_model(F, ve, vp)
%
% Inputs:
%      F_pos, F_neg, K_pos, K_neg - functional parameters of bi-exponential 
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
%
% Concentration model equation
%   Cl_t = F_pos.exp(-t.K_pos) + F_neg.exp(-t.K_neg)) * Cp_t
%
% Where
%   Cp_t = (f_a.Ca_t + f_v.Cv_t) / (1 - Hct)
%
%
% See also: GADOXETATE_MODEL, TWO_CX_MODEL
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
    aoffset = zeros(size(F_pos));
end
if ~exist('voffset', 'var') || isempty(voffset)
    voffset = zeros(size(F_pos));
end
if ~exist('use_fp', 'var') || isempty(use_fp)
    use_fp = 0;
end
K_max = 1e9;

F_pos = double(F_pos(:));
F_neg = double(F_neg(:));
K_pos = double(K_pos(:));
K_neg = double(K_neg(:));
f_a = double(f_a(:));
aoffset = double(aoffset(:));
voffset = double(voffset(:));

%Get number of times and number of voxels, noting we allow model params to
%be scalar (so need to check lengths of all inputs) - this makes it easier
%to fix all params but one, and test a range of inputs on the other to
%observe how that parameter alters the model
n_t = length(dyn_times);
n_vox = max([...
    length(F_pos)
    length(F_neg)
    length(K_pos)
    length(K_neg)
    length(f_a)
    length(aoffset)
    length(voffset)]);

if length(F_pos) ~= n_vox
    if length(F_pos) == 1
        F_pos = ones(n_vox,1)*F_pos;
    else
        error('Incorrect size for F_p, should be 1 x %d', n_vox);
    end
end
if length(F_neg) ~= n_vox
    if length(F_neg) == 1
        F_neg = ones(n_vox,1)*F_neg;
    else
        error('Incorrect size for PS, should be 1 x %d', n_vox);
    end
end
if length(K_pos) ~= n_vox
    if length(K_pos) == 1
        K_pos = ones(n_vox,1)*K_pos;
    else
        error('Incorrect size for v_e, should be 1 x %d', n_vox);
    end
end
if length(K_neg) ~= n_vox
    if length(K_neg) == 1
        K_neg = ones(n_vox,1)*K_neg;
    else
        error('Incorrect size for v_p, should be 1 x %d', n_vox);
    end
end
if length(f_a) ~= n_vox
    if length(f_a) == 1
        f_a = ones(n_vox,1)*f_a;
    else
        error('Incorrect size for f_a, should be 1 x %d', n_vox);
    end
end
if length(aoffset) ~= n_vox
    if length(aoffset) == 1
        aoffset = ones(n_vox,1)*aoffset;
    else
        error('Incorrect size for aoffset, should be 1 x %d', n_vox);
    end
end
if length(voffset) ~= n_vox
    if length(voffset) == 1
        voffset = ones(n_vox,1)*voffset;
    else
        error('Incorrect size for voffset, should be 1 x %d', n_vox);
    end
end

C_t = zeros(n_vox,n_t); 
Cp_t = zeros(n_vox,n_t);
f_v = 1 - f_a; %   estimate of hepatic portal venous fraction
Ft_pos = 0;
Ft_neg = 0;

if use_fp
    F_p = F_pos;
    E_pos = F_neg;
    F_pos = F_p .* E_pos;
    F_neg = F_p .* (1 - E_pos);
end

for i_t = 2:n_t
    
    %Get current time, and time change
    t1 = dyn_times(i_t);
    delta_t = t1 - dyn_times(i_t-1);
    
    %Compute (offset) combined arterial and vascular input for this time
    Ca_ti = interp1(dyn_times, Ca_t, t1 - aoffset, 'linear', 'extrap');
    Cv_ti = interp1(dyn_times, Cv_t, t1 - (aoffset + voffset),...
        'linear', 'extrap');   
    Cp_t(:,i_t) = (f_a.*Ca_ti + f_v.*Cv_ti) / (1 - Hct);
    
    %Update the exponentials for the transfer terms in the two compartments        
    et_pos = exp(-delta_t .* K_pos);
    et_neg = exp(-delta_t .* K_neg);
    
    %Use iterative trick to update the convolutions of transfers with the
    %input function. This only works when the exponent is finite, otherwise
    %the exponential is zero, and the iterative solution is not valid. For
    %these voxels, set A_pos/neg to zero
    A_pos = delta_t * 0.5 * (Cp_t(:,i_t) + Cp_t(:,i_t-1).*et_pos);
    A_pos(K_pos > K_max) = 0;    
    
    A_neg = delta_t * 0.5 * (Cp_t(:,i_t) + Cp_t(:,i_t-1).*et_neg);
    A_neg(K_neg > K_max) = 0;
    
    Ft_pos = Ft_pos.*et_pos + A_pos;
    Ft_neg = Ft_neg.*et_neg + A_neg;
    
    %Combine the two compartments with the rate constant to get the final
    %concentration at this time point
    C_t(:,i_t) = F_pos.*Ft_pos + F_neg.*Ft_neg;  
    
end







