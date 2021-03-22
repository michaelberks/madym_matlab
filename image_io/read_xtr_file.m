function [fields, values] = read_xtr_file(xtr_path)
%READ_XTR_FILE read acquisition parameters from xtr file used in Madym
%   [fields, values] = read_xtr_file(xtr_path)
%
% Inputs:
%      xtr_path - path to xtr file
%
%
% Outputs:
%      fields - cell array of parameter names contained in file
%
%      values - values for each field
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 08-Jan-2019
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 

ff = fopen(xtr_path);
xtr_ss = textscan(ff, '%s');
fclose(ff);
xtr_ss = xtr_ss{1};


if strcmpi(xtr_ss{1}, "voxel")
    [fields, values] = read_old_xtr(xtr_ss);
else
    [fields, values] = read_new_xtr(xtr_ss);
end

%%
function [fields, values] = read_old_xtr(xtr_ss)
%More of a pain, we have a fixed format, with an annoying timestamp
fields = cell(0,1);
values = zeros(0,1);
for i_str = 1:length(xtr_ss)
    str_i = xtr_ss{i_str};
    if str_i(end) == ':'
        field = str_i(1:end-1);
        
        [field, value] = parse_old_fields(field, xtr_ss(i_str+1:end));
        fields = [fields; field]; %#ok
        values = [values; value]; %#ok
        
    end
end

%%
function [fields, values] = read_new_xtr(xtr_ss)
%The new format is easier as everything is a name/value pair. Just need to
%convert the values to numbers
fields = xtr_ss(1:2:end);
values = cellfun(@str2double, xtr_ss(2:2:end));

%%
function [field_out, value] = parse_old_fields(field, values)

switch field
    case 'dimensions'
        field_out = ...
            {'VoxelDimensionsX';'VoxelDimensionsY';'VoxelDimensionsZ'}; 
        value = cellfun(@str2double, values(1:3));
    case 'angle'
        field_out = {'FlipAngle'};
        value = str2double(values{1});
    case 'TR'
        field_out = {field};
        value = str2double(values{1});
    case 'timestamp'
        field_out = {field};
        value = str2double(values{4});
        
    otherwise
        warning(['Field name ' field ' not recognised.']);
        field_out = {field};
        value = str2double(values{1});
end
        
