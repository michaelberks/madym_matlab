function [cmd] = add_option(option_type, cmd, option_str, option)
% 
%  
%   See also STRUCT, VARARGIN.

switch option_type
    
    case 'bool'
        cmd = option_bool(cmd, option_str, option);
        
    case 'int'
        cmd = option_int(cmd, option_str, option);
        
    case 'float'
        cmd = option_float(cmd, option_str, option);
        
    case 'string'
        cmd = option_string(cmd, option_str, option);
        
    case 'int_list'
        cmd = option_int_list(cmd, option_str, option);
        
    case 'float_list'
        cmd = option_float_list(cmd, option_str, option);
        
    case 'string_list'
        cmd = option_string_list(cmd, option_str, option);
        
    otherwise
        error('Option type %s not recognised', option_type);
end

function cmd = option_bool(cmd, option_str, option)
    if isfinite(option)
        cmd = sprintf('%s %s %d', cmd, option_str, option);
    end
    
%--------------------------------------------------------------------------
function cmd = option_int(cmd, option_str, option)
    if isfinite(option)
        cmd = sprintf('%s %s %d', cmd, option_str, option);
    end
    
%--------------------------------------------------------------------------
function cmd = option_float(cmd, option_str, option)
    if isfinite(option)
        cmd = sprintf('%s %s %d', cmd, option_str, option);
    end
    
%--------------------------------------------------------------------------
function cmd = option_string(cmd, option_str, option)
    if isfinite(option)
        cmd = sprintf('%s %s %s', cmd, option_str, option);
    end
    
%--------------------------------------------------------------------------
function cmd = option_int_list(cmd, option_str, option)
    if ~isempty(option)
        options = sprintf('%d', option(1));
        for i_t = 2:length(option)
            options = sprintf('%s,%d', options, option(i_t));
        end
    
        cmd = sprintf('%s %s %s', cmd, option_str, options);
    end
    
%--------------------------------------------------------------------------
function cmd = option_float_list(cmd, option_str, option)
    if ~isempty(option)
        options = sprintf('%d', option(1));
        for i_t = 2:length(option)
            options = sprintf('%s,%d', options, option(i_t));
        end
    
        cmd = sprintf('%s %s %s', cmd, option_str, options);
    end
    
%--------------------------------------------------------------------------
function cmd = option_string_list(cmd, option_str, option)
    if ~isempty(option)
        options = sprintf('%s', option{1});
        for i_t = 2:length(option)
            options = sprintf('%s,%s', options, option{i_t});
        end
    
        cmd = sprintf('%s %s %s', cmd, option_str, options);
    end
  
