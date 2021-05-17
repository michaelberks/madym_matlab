function [status, result] = run_madym(model, output_dir, varargin)
%RUN_MADYM historical wrapper function to call C++ tool Madym. Fits
%   tracer-kinetic models to DCE time-series stored in Analyze format images,
%
% Notes:
%   This is a wrapper to the RUN_MADYM_DCE wrapper, to maintain back
%   compatability.
% 
%
% See also: RUN_MADYM_DCE
%
% Created: 17-May-2021
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk  
% Copyright: (C) University of Manchester 
if ~nargin
    run_madym_DCE();
else
    varargin = [...
        {'model'}, {model}, ...
        {'output_dir'}, {output_dir},...
        varargin];
        
    [status, result] = run_madym_DCE(varargin);
end