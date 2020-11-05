function [residuals_dist] = compute_model_residual_distributions(subject_dir, subject_idx, visit, varargin)
%COMPUTE_MODEL_RESIDUAL_DISTRIBUTIONS *Insert a one line summary here*
%   [] = compute_model_residual_distributions(subject_dir, subject_idx, visit, varargin)
%
% Inputs:
%      subject_dir - *Insert description of input variable here*
%
%      subject_idx - *Insert description of input variable here*
%
%      visit - *Insert description of input variable here*
%
%      varargin - *Insert description of input variable here*
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
% Created: 17-Jul-2020
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 
args = u_packargs(varargin, 0, ...
    'dyn_conc_dir', '',...
    'n_vols', 100,...
    'model_conc_dir', [],...
    'roi_mask_path',[],...
    'save_dir', [],...
    'bins', [],...
    'bin_range', [],...
    'n_bins', [],...
    'plot', 0);

fprintf('Making temporal residuals for %d, visit %d\n', ...
    subject_idx, visit);

%Set up paths from input
dyn_conc_dir = [subject_dir args.dyn_conc_dir '\'];
model_conc_dir = [subject_dir args.model_conc_dir '\'];

%Set up ROI masks
hdr = read_analyze_hdr([dyn_conc_dir 'Ct_sig1.hdr']);
if ~isempty(args.roi_mask_path)
    roi_mask = load_img_volume([subject_dir...
        args.roi_mask_path]) > 0;
else
    roi_mask = true(hdr.Dimensions([2 1 3]));
end

%Load the concentration time-series and compute residuals
C_t = get_dyn_vals([dyn_conc_dir 'Ct_sig'], args.n_vols, roi_mask);
C_t_m = get_dyn_vals([model_conc_dir 'Ct_mod'], args.n_vols, roi_mask);

residuals = C_t - C_t_m;

%Set up bins
if ~isempty(args.bins)
    residuals_dist.bins = args.bins;
else
    if isempty(args.bin_range)
        args.bin_range = ...
            prctile(residuals(isfinite(residuals)), [0.5 99.5]);
    end
    residuals_dist.bins = ...
        linspace(args.bin_range(1), args.bin_range(2), args.n_bins);
end     
bin_width = residuals_dist.bins(2)-residuals_dist.bins(1);
bin_min = residuals_dist.bins(1) - bin_width/2;
bin_max = residuals_dist.bins(end) + bin_width/2;

%Compute NaNs, Inf and lo/hi counts
residuals_dist.n_NaN = sum(isnan(residuals), 2);
inf_idx = isinf(residuals);
residuals(inf_idx) = NaN; %Set inf to NaN so they don't count as hi/lo vals
lo_idx = residuals < bin_min;
hi_idx = residuals > bin_max;
residuals(lo_idx | hi_idx) = NaN;

residuals_dist.n_Inf = sum(inf_idx, 2);
residuals_dist.n_lo = sum(lo_idx, 2);
residuals_dist.n_hi = sum(hi_idx, 2);

%Finally, can take the histogram of the data
residuals_dist.counts = hist(residuals, residuals_dist.bins);

if ~isempty(args.save_dir)
    save_dir = [subject_dir args.save_dir '\'];
    save([save_dir 'residuals_dist.mat'], 'residuals_dist');
end