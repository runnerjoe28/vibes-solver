% Essentially a library for solving forcing functions

% NEED TO IMPLEMENT LATER
function response_func = handle_forcing(cos_matrix, sin_matrix, imp_matrix, other_func)
    
% First verify input format
cos_matrix = check_matrix_form(cos_matrix, "cos");
sin_matrix = check_matrix_form(sin_matrix, "sin");
imp_matrix = check_matrix_form(imp_matrix, "impulse");

% Create variables to track number of forcings and their types
% THIS IS A TODO

response_func = other_func;

end

% Checks for properly formatted sin/cos/impulse matrices
% If it is not formatted, sets matrix to [0] and prints message
%   matrix - cos, sin, or impulse matrix to be verified
%   type - "cos", "sin" or "impulse" to format error message
function new_matrix = check_matrix_form(matrix, type)
    
    % Get the dimensions of the matrix
    [size_x,size_y] = size(matrix);
    
    % Checks based on size
    if (size_y == 2)    % Case 1: Properly formated with values
        new_matrix = matrix;
        return;
    elseif (size_x == 1 && size_y == 1 && matrix(1,1) == 0)     % Case 2: Properly formated with [0]
        new_matrix = matrix;
        return;
    else    % Case 3: Ill formated matrix
        fprintf("The %s matrix is ill-formed, setting it to [0], restart program by ctrl+C if desired\n", type);
        new_matrix = [0];
    end
        
end