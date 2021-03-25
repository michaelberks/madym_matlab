function [volume, header] = load_img_volume(...
    volume_path, scale, flip_y, swap_axes)
%LOAD_IMG_VOLUME Load image from Analyze 7.5 file, applying default
%rotation and scaling
%   [volume] = load_raw_volume(volume_path)
%
% Inputs:
%      volume_path - filepath to image
%
%      scale - scaling applied to image values (default 1)
%
%      flip_y - if true, flip the y-axis (applied *after* swap_axes)
%      (default true)
%
%      swap_axes - swap the x,y axis (default true)
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
% See also: SAVE_IMG_VOLUME, READ_ANALYZE_IMG, READ_ANALYZE_HDR
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
if ~exist('swap_axes', 'var') || isempty(swap_axes)
    swap_axes = true;
end

if exist('niftiread', 'file')
    try
        header = niftiinfo(volume_path);
    catch
        header = [];
    end
    volume = double(niftiread(volume_path)) / scale;
else
    header = read_analyze_hdr(volume_path);
    volume = double(read_analyze_img([],header)) / scale;
end

if swap_axes
   volume = permute(volume, [2 1 3]);
end

if flip_y
    try
        volume = flip(volume,1);
    catch
        volume = volume(end:-1:1,:,:);
    end
end
