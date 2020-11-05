function [] = display_DIBEM_factor_counts(...
    indicator_map, roi_masks, noise_map, txt_file)
%DISPLAY_DIBEM_FACTOR_COUNTS print summary information on DIBEM factors to
%screen or text file
%   [] = display_DIBEM_factor_counts(indicator_map, roi_masks, tse, txt_file)
%
% Inputs:
%      indicator_map - map of DIBEM factors
%
%      roi_masks - Set of ROI to aggregate factors
%
%      noise_map - If given, displays the mean noise for each factor
%
%      txt_file - If given, output will be written to this path, otherwise
%      it will be printed to the screen (fid = 1)
%
%
% Outputs:
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
if ~exist('noise_map', 'var')
    noise_map = [];
end
if exist('txt_file', 'var') && ~isempty(txt_file)
    fid = fopen(txt_file, 'wt');
else
    fid = 1;
end

%Get unique factors and counts for each, then sort by the count in the
%first ROI
[uni_factors, uni_counts, tse_scores] = count_DIBEM_factors(...
    indicator_map, roi_masks, noise_map);
[uni_counts, sort_idx] = sortrows(uni_counts, 1, 'd');

%Convert counts to pcts
uni_pcts = zeros(size(uni_counts));
for i_roi = 1:size(roi_masks,4)
    n_mask = nnz(roi_masks(:,:,:,i_roi));
    uni_pcts(:,i_roi) = round(100*uni_counts(:,i_roi) / n_mask);
end

%Print pcts for each factor either to screen or text file
fprintf(fid, 'Indicator map factors:\n');
for i_u = 1:length(uni_factors)
    if ~any(uni_pcts(i_u,:))
        continue;
    end

    factor_value = uni_factors(sort_idx(i_u));
    factor_names = DIBEM_indicator_factors(factor_value);
    fprintf(fid, '   % 5d: ', factor_value);
    fprintf(fid, '% 3d  ', uni_pcts(i_u,:));
    fprintf(fid, '%0.2e  ', tse_scores(sort_idx(i_u)));
    for i_f = 1:length(factor_names)
        fprintf(fid, ', %s', factor_names{i_f});
    end
    fprintf(fid, '\n');
end

if fid ~= 1
    fclose(fid);
end