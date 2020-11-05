function [] = write_xtr_file(xtr_path, append, varargin)
%WRITE_XTR_FILE write a text file of name/value pairs, as used in Madym to
%specify additional information not contained in Analyze75 img headers (eg.
%scan flip angle, TR etc)
%   [] = write_xtr_file(xtr_path, varargin)
%
% Inputs:
%      xtr_path - path to write xtr file (typically with extension .xtr)
%
%      append - If file already exists, append or overwrite
%
%      varargin - List of fieldname/value pairs
%
%
% Outputs:
%
% Example: write_xtr_file('temp.xtr', 'TR', 2.4, 'FlipAngle', 20.0, 'TimeStamp', 12345)
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
num_args = length(varargin);
if rem(num_args, 2)
    error('Arguments should be supplied as name/value pairs')
end


if exist(xtr_path, 'file') && append
    new_fields = varargin(1:2:end);
    new_values = varargin(2:2:end);
    [fields, values] = read_xtr_file(xtr_path);
    values = num2cell(values);
    
    %Loop through the new fields, if a match is found in existing fields,
    %update the value, otherwise append the new field and value
    for i_val = 1:length(new_fields)
        match_idx = strcmpi(fields, new_fields{i_val});
        if any(match_idx)
            values(match_idx) = new_values(i_val);
        else
            fields(end+1) = new_fields(i_val); %#ok
            values(end+1) = new_values(i_val); %#ok
        end
    end
else
    fields = varargin(1:2:end);
    values = varargin(2:2:end);
end

%Open a file identifier, either in write or append mode
fid = fopen(xtr_path, 'wt');

%Write out each name/value pair
for i_val = 1:length(fields)
    
    fprintf(fid, '%s%', fields{i_val});
    for j_val = 1:length(values{i_val})
        fprintf(fid, ' %f', values{i_val}(j_val));
    end
    fprintf(fid, '\n');
end

%Close the file ID
fclose(fid);


