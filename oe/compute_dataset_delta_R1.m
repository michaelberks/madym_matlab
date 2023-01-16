function [delta_R1, R1_baseline, R1_oe] = compute_dataset_delta_R1(...
    visit_dir, oe_limits, varargin)
%BIOCHECC_COMPUTE_DELTA_R1 compute delta R1 in ROI at visits 1 + 2 for a
%given BIOCHECC subject
%   delta_R1 = biochecc_compute_delta_R1(...
%       subject_dir, roi_path, oe_limits, t_border)
%
% Inputs:
%      visit_dir - root folder for the BIOCHECC subject. This is expected
%      to contain folder visit1 and visit2
%
%      roi_path - if a string, defines filepath to an ROI mask image to be
%      loaded. Alternatively this can be the mask itself, in which case it
%      must be a logical array the same size as the OE/DCE data
%
%      oe_limits - 4 element array specifiying the timepoints for:
%       baseline_start, baseline_end, oxygen_start, oxygen_end. Delta R1
%       will be computed as the difference between the mean R1 for 
%       oxygen_start:oxygen_end and baseline_start:baseline_end 
%
%      t_border - int, used to create an expanded cuboid ROI about the
%      tumour (or reference for HV) ROI.
%
% Outputs:
%      delta_R1 - 2 element structure array, with results for visits 1 and
%      2. Each struct contains:
%      - tx1, tx2, ty1, ty2, tz1, tz2: the start and end indices in
%        (x,y,z) for the expanded cuboid ROI
%      - slice: the central slice of the ROI
%      - roi_median: median of delta R1 in the ROI
%      - R1_region: delta R1 for each voxel in the expanded ROI
%      - R1_slice: delta R1 for each voxel in the central slice
%      - to complete...
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: Oct-2022
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester
% Unpack the arguments:
args = u_packargs(varargin, 0, ...
    'oe_dir', 'oxygen',...
    'oe_name', 'oe_',...
    'oe_index_format', '%03d',...
    'T1_path', '',...
    'efficiency_path', '',...
    'average_fun', 'median',...
    'debug', 0);
clear varargin;

%Get average function from user option
switch args.average_fun
    case 'mean'
        avg_fun = @mean;
    case 'median'
        avg_fun = @median;
    otherwise
        error('Average function %s not recognised', average_fun);
end

%Bug out if no data
oe_dir = fullfile(visit_dir, args.oe_dir);
if ~exist(oe_dir, 'dir')
    warning('%s missing oxygen folder', visit_dir);
    return
end

%Load in the dynamic time-series
n_oe = length(dir(fullfile(oe_dir, '*.nii.gz'))) - 1;
oe_vols = get_dyn_vols(...
    fullfile(oe_dir, args.oe_name), n_oe, 0, args.oe_index_format);

%Compute average baseline signal
S_0 = avg_fun(oe_vols(:,:,:,oe_limits(1):oe_limits(2)), 4);

%Compute dynamic R1
xtr_path = fullfile(oe_dir, ...
    sprintf(['%s' args.oe_index_format '.xtr'], args.oe_name, 1));
R1_t = dynamic_R1_oxy(...
    visit_dir, oe_vols, S_0, args.T1_path, args.efficiency_path, xtr_path);


%Compute average at baseline and in enhancing period
R1_baseline = avg_fun(R1_t(:,:,:,oe_limits(1):oe_limits(2)), 4);
R1_oe = avg_fun(R1_t(:,:,:,oe_limits(3):oe_limits(4)), 4);
delta_R1 = R1_oe - R1_baseline;
   
end
%%

function R1_t = dynamic_R1_oxy(...
    visit_dir, dyn_vols, S_0, T1_path, efficiency_path, xtr_path)
    T1_0 = load_img_volume(fullfile(visit_dir, T1_path));
    
    if isempty(efficiency_path)
        eff = [];
    else
        eff = load_img_volume(fullfile(visit_dir, efficiency_path));
    end

    TI = get_oe_TI(xtr_path);
    R1_t = dynamic_R1_IR(dyn_vols, T1_0, S_0, TI, eff);
end

function TI = get_oe_TI(xtr_path)
    [fields, values] = read_xtr_file(xtr_path);
    TI = values(strcmpi('TI', fields));
end



