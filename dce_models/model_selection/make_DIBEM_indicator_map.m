function [indicator_map] = make_DIBEM_indicator_map(...
    active_model, cxm_model, selected_irf, ...
    accumulation_map_model, accumulation_map_signal)
%MAKE_DIBEM_INDICATOR_MAP create map where each voxel indicates which
%physiological factors were present in a DIBEM fit
%   [indicator_map] = make_DIBEM_indicator_map(model, aic_preferred, accumulation_map)
%
% Inputs:
%      active/cxm_model - structure containing active and CXM derivations of DIBEM
%      parameters, each of which is a 3-D array
%
%      selected_irf - 3-D array, same size as model parameters specifying
%      which IRF form was preferred at each voxel
%
%      accumulation_map - 3-D logical array, positive where model
%      accumlates CA from 1 to 5 minutes post contrast
%
%
% Outputs:
%      indicator_map - 3-D array, each voxel is an integer with values 0 to
%      2^n_factors. Converting this value to a binary string indicates
%      which factors were present.
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
% factors = {
%     'i2'
%     'i3'
%     'i4'
%     'accum'
%     'accum_signal'
%     'vecs > 1'
%     've + vp > 1'
%     'kef > ki'
%     'active params < 0'
%     'CXM params < 0'
%     'f_a <= 0'
%     'f_a >= 1'};

%Make maps for factors from the model parameters
vecs_map = active_model.v_ecs > 1 | ~isfinite(active_model.v_ecs);
vevp_map = cxm_model.v_e + cxm_model.v_p > 1 | ...
    ~isfinite(cxm_model.v_e) | ~isfinite(cxm_model.v_p);
kefki_map = active_model.k_ef > active_model.k_i;
fa0_map = cxm_model.f_a <= 0;
fa1_map = cxm_model.f_a >= 1;

active_invalid = ...
    active_model.F_p < 0 | ...
    active_model.v_ecs < 0 | ...
    active_model.k_i < 0 | ...
    active_model.k_ef < 0;

cxm_invalid = ...
    cxm_model.F_p < 0 | ...
    cxm_model.PS < 0 | ...
    cxm_model.v_e < 0 | ...
    cxm_model.v_p < 0;

%For each factor, add the factor value to the voxels where this factor is
%present
indicator_map = zeros(size(selected_irf));
indicator_map(selected_irf == 1) = ...
    indicator_map(selected_irf == 1) + 2^0;
indicator_map(selected_irf == 2) = ...
    indicator_map(selected_irf == 2) + 2^1;
indicator_map(selected_irf == 3) = ...
    indicator_map(selected_irf == 3) + 2^2;
indicator_map(accumulation_map_model) = ...
    indicator_map(accumulation_map_model) + 2^3;
indicator_map(accumulation_map_signal) = ...
    indicator_map(accumulation_map_signal) + 2^4;
indicator_map(vecs_map) = indicator_map(vecs_map) + 2^5;
indicator_map(vevp_map) = indicator_map(vevp_map) + 2^6;
indicator_map(kefki_map) = indicator_map(kefki_map) + 2^7;
indicator_map(active_invalid) = indicator_map(active_invalid) + 2^8;
indicator_map(cxm_invalid) = indicator_map(cxm_invalid) + 2^9;
    indicator_map(fa0_map) = indicator_map(fa0_map) + 2^10;
    indicator_map(fa1_map) = indicator_map(fa1_map) + 2^11;