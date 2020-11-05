function [] = save_img_volume(volume, volume_path, vox_sz, scale, flip_y, ...
    to_int, flip_x, swap_axes, int_fmt)
%save_img_volume save 3D image in Analyze 7.5 format, applying default
%rotation and scaling
%   [volume] = load_raw_volume(volume_path)
%
% Inputs:
%      volume - 3D array containing image
%
%      volume_path - filepath to image
%
%      vox_sz - 3-element vector specifying dimensions in mm of voxels
%
%      scale - scaling applied to image values (default 1)
%
%      flip_y - if true, flip the y-axis (applied *after* swap_axes)
%      (default true)
%
%      to_int - if true, convert to integer before saving (default false)
%
%      swap_axes - swap the x,y axis (default true)
%
%      int_fmt - integer format to convert to (default int16, ignored if
%      to_int false)
%
%
% Outputs:
%      volume - 3D array containing image
%
%      header - Analyze header as Matlab struct
%
%
% Example:
%
% Notes:
%
% See also: LOAD_IMG_VOLUME, WRITEANALYZE
%
% Created: 11-Apr-2017
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester
if ~exist('scale', 'var') || isempty(scale)
    scale = 1;
end
if ~exist('flip_y', 'var') || isempty(flip_y)
    flip_y = true;
end
if ~exist('to_int', 'var') || isempty(to_int)
    to_int = false;
end
if ~exist('flip_x', 'var') || isempty(flip_x)
    flip_x = false;
end
if ~exist('swap_axes', 'var') || isempty(swap_axes)
    swap_axes = true;
end
if ~exist('int_fmt', 'var') || isempty(int_fmt)
    int_fmt = 'int16';
end

if flip_y
    try
        volume = flip(volume,1);
    catch
        volume = volume(end:-1:1,:,:);
    end
end

if flip_x
    try
        volume = flip(volume,2);
    catch
        volume = volume(:,end:-1:1,:);
    end
end

volume = single(volume) / scale;
if to_int
    switch int_fmt
        case 'int16'
            volume = int16(volume);
        case 'uint8'
            volume = uint8(volume);
        otherwise
            error('Format %s not supported, must be int16 or uint8', int_fmt);
    end
end

if swap_axes
   volume = permute(volume, [2 1 3]);
end

writeanalyze(volume, volume_path, vox_sz);




