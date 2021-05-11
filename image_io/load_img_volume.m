function [volume, header] = load_img_volume(...
    volume_path, scale, flip_y, swap_axes, use_inbuilt)
%LOAD_IMG_VOLUME Load image from Analyze 7.5 or NIFTI file, applying default
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
%      use_inbuilt - use Matlab's niftiread if it is available (default true)
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
% Notes: We try to use Matlab's inbuilt niftiread functions if the exist,
% however these don't seem to work for all Analyze formats. If they won't
% work (either causing an error or returning empty), we fall back to our
% own read_analayze_img, however this will only work for .hdr/img Analyze
% 7.5 format, and not for later NIFTI nii or nii.gz 
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
if ~exist('use_inbuilt', 'var') || isempty(use_inbuilt)
    use_inbuilt = true;
end

% The steps we take:

if use_inbuilt && exist('niftiread', 'file')
    % 1) If use_inbuilt selected, and niftiread available (>2018b) try loading
    % using niftiread/niftiinfo
    try
        [volume, header] = load_inbuilt(volume_path);
        
        % 2) If either volume or header are empty, try using own
        % read_analyze functions - these will error if nifti format, but we
        % have nothing else to try so just let this error happen
        if isempty(volume) || isempty(header)
            [volume, header] = load_madym(volume_path);
        end
        
    catch err
        
        % 3) If niftiread/info throws an error, again try our own
        % read_analyze, but this time, only try our functions for matching
        % file extensions, otherwise we're better off throwing the original
        % error
        [~,~,ext] = fileparts(volume_path);
        if ismember(lower(ext), {'', '.img', '.hdr'})
            [volume, header] = load_madym(volume_path);
        else
            rethrow(err);            
        end
    end
    
else
    % 4) If inbuilt functions not available, try our own. Again, these will
    % error on .nii and .nii.gz files, but nothing else to try so just let
    % this error happend
    [volume, header] = load_madym(volume_path);
end

%Now apply scaling and axes swaps/flips
volume = double(volume) / scale;

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

%End of main function

%Load volume and header using Matlab's inbuilt fucntion
function [volume, header] = load_inbuilt(volume_path)
    header = niftiinfo(volume_path);
    volume = niftiread(volume_path);
    
%Load volume and header using our own functions
function [volume, header] = load_madym(volume_path)
    header = read_analyze_hdr(volume_path);
    volume = read_analyze_img([],header);
    
