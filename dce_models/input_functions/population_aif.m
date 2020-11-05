function [aif] = population_aif(dyn_times, prebolus_im, varargin)
%POPULATION_AIF Compute population AIF using Geoff Parker's analytic method
%   [aif] = population_aif(dyn_times, aif, offset)
%
% Inputs:
%      dyn_times - [n_t,1] times at which to compute AIF values. Assumed to
%      be in minutes. If given in seconds, set option convert_to_minutes to
%      true
%
%      prebolus_im - int, the image at which the bolus is injected, given
%      as an index into dyn_times
%
% Options, given as name/values pairs using the u_packargs interface
%      offset - double, 0, delay time from injection
%
%      dose - double, 0.1, dose of agent in mmol/kg
%
%      convert_to_minutes - bool, false, if true, divides dyn_times by 60
%      to convert from seconds to minutes
%
%      hct - double, 0.42, haematocrit correction applied to scale whole
%      blood to plasma
%
%
% Outputs:
%      aif - [n_t,1] arterial input function values at each time-point
%
%
% Example:
%
% Notes: TODO add reference to Geoff's paper
%
% See also: LOAD_AIF, SAVE_AIF, COMPUTE_PIF
%
% Created: 19-Oct-2017
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 
args = u_packargs(varargin, 0, ...
    'offset', 0, ...
    'dose', 0.1,...
    'convert_to_minutes', false,... %Assume this has already been done
    'hct', 0.42);
clear varargin;

n_times = length(dyn_times);
aif = zeros(n_times, 1);
prebolus_time = dyn_times(prebolus_im);

if args.convert_to_minutes
    to_mins = 60;
else
    to_mins = 1;
end

for i_t = 1:n_times
    %Convert time to minutes (has this already been done in dyn_times?
    t = (dyn_times(i_t) - args.offset - prebolus_time) / to_mins;
	
    %A1/(SD1*sqrt(2*PI)) * exp(-(t-m1)^2/(2*var1))
	%A1 = 0.833, SD1 = 0.055, m1 = 0.171
    gaussian1 = 5.73258 * exp(...
        -1.0 * ...
        (t - 0.17046) * (t - 0.17046) /... 
        (2.0 * 0.0563 * 0.0563) );
    
    %A2/(SD2*sqrt(2*PI)) * exp(-(t-m2)^2/(2*var2))
    %A2 = 0.336, SD2 = 0.134, m2 = 0.364
    gaussian2 = 0.997356 * exp(...
        -1.0 * ...
        (t - 0.365) * (t - 0.365) / ...
        (2.0 * 0.132 * 0.132));
   % alpha*exp(-beta*t) / (1+exp(-s(t-tau)))
   % alpha = 1.064, beta = 0.166, s = 37.772, tau = 0.482
   sigmoid = 1.050 * exp(-0.1685 * t) / (1.0 + exp(-38.078 * (t - 0.483)));
   aif(i_t) = ((args.dose / 0.1) * (gaussian1 + gaussian2 + sigmoid)) /...
       (1.0 - args.hct);
end
