function [] = make_postcon_means(...
    output_dir, dynamic_dir, dyn_times, injection_time, ...
    temporal_mean_ranges, temporal_mean_names, do_plot)
%MAKE_POSTCON_MEANS *Insert a one line summary here*
%   [] = georgiou_aif(subject_dir, subject_idx, visit, varargin)
%
%
% Mandatory Arguments:
%
% Optional Arguments:
%
% Outputs:
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 02-Jul-2020
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester
if ~exist('do_plot', 'var')
    do_plot = 0;
end

%Set up path to dynamic data
dyn_header = read_analyze_hdr([dynamic_dir 'dyn_1.hdr']);

%Set times to injection image
t_in_s = 60*(dyn_times - injection_time);

%Loop through each range creating the temporal mean
n_means = size(temporal_mean_ranges, 1);
create_folder(output_dir);

for i_m = 1:n_means
    %Find the image number that best matches the start and end time for
    %this mean
    [~,start_im] = min(abs(t_in_s - temporal_mean_ranges(i_m,1)));
    [~,end_im] = min(abs(t_in_s - temporal_mean_ranges(i_m,2)));
    n_ims = (end_im - start_im + 1);
    for i_t = start_im:end_im
        
        d = load_img_volume([dynamic_dir 'dyn_' num2str(i_t) '.hdr']);
        
        if i_t == start_im
            dyn_sum = d;
        else
            dyn_sum = dyn_sum + d;
        end
    end
    dyn_mean = dyn_sum / n_ims;
    
    %Get  dyn name
    mean_name = temporal_mean_names{i_m};
    if isempty(mean_name)
        sprintf('dyn_mean_postcon%d-%ds', temporal_mean_ranges(i_m,:));
    end
    mean_path = [output_dir '/' mean_name '.hdr'];
    
    save_img_volume(dyn_mean, mean_path,  ...
        dyn_header.PixelDimensions, [], [], 0);
    
    fprintf('Saved temporal mean %s\n', mean_path);
    if do_plot
        figure;
        imgray(dyn_mean(:,:,round(end/2)));
    end
end