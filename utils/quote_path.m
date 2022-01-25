function [quoted_path_str] = quote_path(path_str)
%QUOTE_PATH Enclose path in double quotes, helps deal with paths containing
%spaces. 
%   [quoted_path_str] = quote_path(path_str)
%
% Inputs:
%      path_str - filepath to be quote enclosed
%
% Outputs:
%      quoted_path_str - filepath enclosed in double quotes
%
% Example:
%      [quoted_path_str] = quote_path('C:\a path with spaces')
%      returns '"C:\a path with spaces"'
% Notes:
%       Any existing double-quotes are first removed from the path to
%       ensure they don't interfere with the new enclosing quotes
%
path_str(path_str == '"') = [];
quoted_path_str = ['"' path_str '"'];
  
