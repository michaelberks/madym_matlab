function [uni_factors, uni_counts, noise_map_scores] = count_DIBEM_factors(...
    indicator_map, roi_masks, noise_map)
%COUNT_DIBEM_FACTORS *Insert a one line summary here*
%   [uni_counts, tse_scores] = count_DIBEM_factors(indicator_map, roi_masks, tse)
%
% Inputs:
%      indicator_map - 3-D array, map of DIBEM factors
%
%      roi_masks - Set of ROIs (stacked on 4th dimension) to aggregate factors
%
%      noise_map - A noise value associated with each voxel, must match
%      size of indicator_map
%
%
% Outputs:
%      uni_factors - list of unique factor values found in ht eindicator
%      map
%
%      uni_counts - 2-d array n_uni x n_rois, each column containing the 
%      number of voxels for each factor
%
%      noise_map_scores - If a noise map is given, the mean score for each
%      factor
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 30-Jun-2020
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 
%Get list of unique factor values
uni_factors = unique(indicator_map(any(roi_masks,4)));
n_uni = length(uni_factors);

%set up containers for counts of each unique value
n_rois = size(roi_masks, 4);
uni_counts = zeros(n_uni, n_rois);
noise_map_scores = zeros(n_uni, 1);

%Loop through ROIs and count number of voxels for each factor
for i_roi = 1:n_rois
    roi_mask = roi_masks(:,:,:,i_roi);
 
    for i_u = 1:n_uni   
        uni_counts(i_u,i_roi) = ...
            nnz(indicator_map(roi_mask)==uni_factors(i_u)); 
    
        if exist('noise_map', 'var') && ~isempty(noise_map)
            noise_map_scores(i_u,1) =...
                mean(noise_map(indicator_map==uni_factors(i_u)));
        end
    end
end