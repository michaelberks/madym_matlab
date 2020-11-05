function [dyn_times, dyn_FA, dyn_TR] = get_dyn_times(root_path, num_vols, index_fmt)
%GET_DYN_TIMES get meta information (scan time, FA, TR) from xtr files for 
% folder of dynamic volumes
%   [times] = get_dyn_times(root_path, index_fmt, num_vols)
%
% Inputs:
%      root_path - folder + filename root where volumes are
%
%      num_vols - number of volumes to load
%
%
%      index_fmt ('%01u') - format that converts indexes into volume suffix
%
%
% Outputs:
%      dyn_times - *Insert description of input variable here*
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

if ~exist('index_fmt', 'var') || isempty(index_fmt)
    index_fmt = '%01u';
end

dyn_times = zeros(num_vols,1);
dyn_TR = zeros(num_vols,1);
dyn_FA = zeros(num_vols,1);

for i_vol = 1:num_vols
    vol_path = [root_path sprintf(index_fmt, i_vol) '.xtr'];
    ff = fopen(vol_path);
    xtr_ss = textscan(ff, '%s');
    fclose(ff);
    xtr_ss = xtr_ss{1};
    
    if strcmpi(xtr_ss{1}, 'voxel')
        [dyn_times(i_vol), dyn_FA(i_vol), dyn_TR(i_vol)] =...
            read_xtr_old(xtr_ss);
    else
        [dyn_times(i_vol), dyn_FA(i_vol), dyn_TR(i_vol)] =...
            read_xtr_new(xtr_ss);
    end
    
end
dyn_times = (dyn_times - dyn_times(1))/60;
% ------------------ End of main function --------------------------------
%% -----------------------------------------------------------------------
function [dyn_t, dyn_FA, dyn_TR] = read_xtr_old(xtr_ss)

    %Find index of time stamp
    t_idx = find(strcmpi(xtr_ss, 'timestamp:'),1);
    
    if isempty(t_idx)
        t_idx = find(strcmpi(xtr_ss, 'timestamp'),1);
    end
    if isempty(t_idx)
        error(['Timestamp not found for volume ' num2str(i_vol)]);
    end
    
    hh = str2double(xtr_ss{t_idx+1});
    mm = NaN;
    ss = NaN;
    if t_idx < length(xtr_ss)-1
        mm = str2double(xtr_ss{t_idx+2});
        ss = str2double(xtr_ss{t_idx+3});
    end
    if isnan(mm) || isnan(ss)
        timeStamp = hh;
        hh = floor(timeStamp / 10000);
		mm = floor((timeStamp - 10000 * hh) / 100);
		ss = (timeStamp...
            - 10000 * hh...
            - 100 * mm);
    end

    dyn_t = 60*60*hh + 60*mm + ss;
    
    %Find index of flip angle
    t_idx = find(strcmpi(xtr_ss, 'angle:'),1);    
    if isempty(t_idx)
        t_idx = find(strcmpi(xtr_ss, 'FlipAngle'),1);
    end
    if isempty(t_idx)
        error(['Flip angle not found for volume ' num2str(i_vol)]);
    end
    dyn_FA = str2double(xtr_ss{t_idx+1});
    
    %Find index of TR
    t_idx = find(strcmpi(xtr_ss, 'TR:'),1);    
    if isempty(t_idx)
        t_idx = find(strcmpi(xtr_ss, 'TR'),1);
    end
    if isempty(t_idx)
        error(['TR not found for volume ' num2str(i_vol)]);
    end
    dyn_TR = str2double(xtr_ss{t_idx+1});
    
%% ------------------------------------------------------------------------    
function [dyn_t, dyn_FA, dyn_TR] = read_xtr_new(xtr_ss)

    %Find index of time stamp
    t_idx = find(strcmpi(xtr_ss, 'timestamp:'),1);
    
    if isempty(t_idx)
        t_idx = find(strcmpi(xtr_ss, 'timestamp'),1);
    end
    if isempty(t_idx)
        error(['Timestamp not found for volume ' num2str(i_vol)]);
    end
    
    hh = str2double(xtr_ss{t_idx+1});
    mm = NaN;
    ss = NaN;
    if t_idx < length(xtr_ss)-1
        mm = str2double(xtr_ss{t_idx+2});
        ss = str2double(xtr_ss{t_idx+3});
    end
    if isnan(mm) || isnan(ss)
        timeStamp = hh;
        hh = floor(timeStamp / 10000);
		mm = floor((timeStamp - 10000 * hh) / 100);
		ss = (timeStamp...
            - 10000 * hh...
            - 100 * mm);
    end
    dyn_t = 60*60*hh + 60*mm + ss;
    
    %Find index of flip angle
    t_idx = find(strcmpi(xtr_ss, 'angle:'),1);    
    if isempty(t_idx)
        t_idx = find(strcmpi(xtr_ss, 'FlipAngle'),1);
    end
    if isempty(t_idx)
        error(['Flip angle not found for volume ' num2str(i_vol)]);
    end
    dyn_FA = str2double(xtr_ss{t_idx+1});
    
    %Find index of TR
    t_idx = find(strcmpi(xtr_ss, 'TR:'),1);    
    if isempty(t_idx)
        t_idx = find(strcmpi(xtr_ss, 'TR'),1);
    end
    if isempty(t_idx)
        error(['TR not found for volume ' num2str(i_vol)]);
    end
    dyn_TR = str2double(xtr_ss{t_idx+1});

    

