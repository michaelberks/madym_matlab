function [St] = signal_from_T1(T1, S0, FA, TR)
%SIGNAL_FROM_T1 compute MR signal given T1 and flip angle
%   [St] = signal_from_T1(T1, S0, FA, TR)
%
% Inputs:
%      T1 - input T1, all inputs must be scalaras or arrays of the same
%      size
%
%      S0 - signal at baseline
%
%      FA - flip angle
%
%      TR - tissue relaxivity constant
%
%
% Outputs:
%      St - MR signal
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 25-Oct-2017
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 
E1 = exp(-TR ./ T1);
try
    St = S0 .* sind(FA) .* (1 - E1) ./ ...
        (1 - cosd(FA) .* E1);
catch %For versions prior to automatic dimension casting (<R2017)
    St = bsxfun(@times, bsxfun(@times, S0, sind(FA)), (1 - E1)) ./ ...
        (1 - bsxfun(@times, cosd(FA),E1));
end
