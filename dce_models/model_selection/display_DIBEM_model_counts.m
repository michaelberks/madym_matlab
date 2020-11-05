function [] = display_DIBEM_model_counts(aic_counts, model_names, n_roi_voxels, roi_names, txt_file)
%DISPLAY_DIBEM_MODEL_COUNTS *Insert a one line summary here*
%   [] = display_DIBEM_model_counts(aic_selected, roi_masks, txt_file)
%
% Inputs:
%      aic_selected - *Insert description of input variable here*
%
%      roi_masks - *Insert description of input variable here*
%
%      txt_file - *Insert description of input variable here*
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
% Created: 01-Jul-2020
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester
if exist('txt_file', 'var') && ~isempty(txt_file)
    fid = fopen(txt_file, 'wt');
else
    fid = 1;
end

fprintf('                    ');
for i_model = 1:length(model_names)
    fprintf('%s ', model_names{i_model});
end
n_rois = length(n_roi_voxels);

for i_roi = 1:n_rois
    fprintf('\n%s proportions:  ', roi_names{i_roi});
    fprintf('%2.1f ', ...
        100*aic_counts(:,i_roi) / n_roi_voxels(i_roi));
end
fprintf('\n\n');

if fid ~= 1
    fclose(fid);
end