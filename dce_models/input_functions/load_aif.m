function [aif] = load_aif(aif_path, h_correction)
%LOAD_AIF load AIF from file
%   [aif] = load_aif(aif_path, h_correction)
%
% Inputs:
%      aif_path - (str) Path to file storing AIF values
%
%      h_correction - (double) haematocrit correction, default 0. If
%      non-zero AIF will be divided by (1 - h_correction)
%
%
% Outputs:
%      aif - 2 x n_t array, 1st column are dynamic times (assumed to be in
%      minutes), 2nd column AIF values
%
%
% Example:
%
% Notes: File format as used by madym, a simple text file with 2 columns of
% data, the first containing times (in minutes), the second containing the
% AIF values
%
% See also:
%
% Created: 12-Apr-2017
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 
if ~exist('h_correction', 'var') || isempty(h_correction)
    h_correction = 0; %0.42%;
end

aif = load(aif_path);
if h_correction
    aif(:,2) = aif(:,2) / (1 - h_correction);
end
