function [] = run_madym_tests(test_level)
%RUN_MADYM_TESTS run set of tests tests on C++ madym tools.
%   [] = run_madym_tests()
%
% Inputs:
%       test_level (int, >= 1) - set level of tests required.
%           1 = Basic tests only (currently just this implemented)
%
% Outputs:
%
% Example:
%
% Notes:
%
% See also: INSTALL_MADYM
%
% Created: 01-May-2019
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester
if ~exist('test_level', 'var') || isempty(test_level)
    test_level = 1;
end

fprintf('\n*****************************************************\n');
fprintf('Running madym tests\n');
fprintf('*****************************************************\n');

%Calling functions with no input runs basic tests
run_madym_lite();
run_madym();
run_madym_T1();
run_madym_DWI_lite();

if test_level > 1
    %Apply extended tests
    
end

fprintf('\n*****************************************************\n');
fprintf('All tests passed\n');
fprintf('*****************************************************\n');


