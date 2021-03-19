function [dyn_vols, dyn_headers] = get_dyn_vols(root_path, num_vols, apply_smoothing, index_fmt)
%GET_DYN_VOLS load dynamic series into 4D data array
%   [times] = get_dyn_vols(root_path, index_fmt, num_vols)
%
% Inputs:
%      root_path - path to each volume to be loaded, so that the full path
%      to the i-th volume is [root_path sprintf(index_fmt, i) '.hdr']
%
%      index_fmt - defines format for indexing of volumes. In most QBi data
%      this is just %01u (equivalent to num2str, with no modifiers)
%
%      num_vols - number of volumes to load
%
%      apply_smoothing - flag to apply tangential smoothing to each volume
%
%
% Outputs:
%      dyn_vols - (Ny, Nx, Nz, Nvols) array containing each dynamic volume
%
%      dyn_headers - header information for each volume
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 29-Mar-2017
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 
if ~exist('apply_smoothing', 'var') || isempty(apply_smoothing)
    apply_smoothing = false;
end

if ~exist('index_fmt', 'var') || isempty(index_fmt)
    index_fmt = '%01u';
end

if nargout > 1
    load_headers = true;
else
    load_headers = false;
end

for i_vol = 1:num_vols
    
    vol_path = [root_path sprintf(index_fmt, i_vol) '.hdr'];
    if load_headers        
        [d, dyn_header] = load_img_volume(vol_path);
    else
        d = load_img_volume(vol_path);
    end
    
    if i_vol == 1
        [n_y, n_x, n_z] = size(d);
        dyn_vols = zeros(n_y, n_x, n_z, num_vols);
        
        if load_headers
            dyn_headers = cell(num_vols,1);
        end
    end
    if apply_smoothing
        d = tangential_smoothing_vol(d);
    end
    dyn_vols(:,:,:,i_vol) = d;
    
    if load_headers
        dyn_headers{i_vol} = dyn_header;
    end

end
