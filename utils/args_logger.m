function [info] = args_logger(args_file, delete_file)
%ARGS_LOGGER create log structure of arguments
%   [args] = args_logger()
%
% Inputs:
%       args_file - path to saved args or the args themselves in a
%       structure
%
%       delete_file - flag to delete args_file after reading, default true
% Outputs:
%      info - structure logging the calling function, args, username, computername
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 14-Dec-2022
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester
info = struct();
stack = dbstack();
info.calling_function = which(stack(2).name);
info.datetime = datetime('now');
info.username = get_username();
info.computername = get_computername();

if ischar(args_file)
    info.args = load(args_file);
    for field = fieldnames(info.args)'
        try
            if strcmpi (info.args.(field{1}), args_file)
                info.args = rmfield(info.args, field{1});
            end
        catch
        end
    end
    if ~exist('delete_file', 'var') || delete_file
        delete(args_file);
    end
else
    info.args = args_file;
end







