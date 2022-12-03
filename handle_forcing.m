% Essentially a library for solving forcing functions

% NEED TO IMPLEMENT LATER
function response_func = handle_forcing(m, c, k, cos_matrix, sin_matrix, imp_matrix, other_func)
clc;
% First verify input format
cos_matrix = check_matrix_form(cos_matrix, "cos");
sin_matrix = check_matrix_form(sin_matrix, "sin");
imp_matrix = check_matrix_form(imp_matrix, "impulse");

% Create variables to track number of forcings and their types
num_cos_forcings = count_matrix(cos_matrix);
num_sin_forcings = count_matrix(sin_matrix);
num_imp_forcings = count_matrix(imp_matrix);
has_extra_forcing = 1;  % 1 if other func is non-zero, 0 is 
if (other_func == 0)
    has_extra_forcing = 0;
end

% Ask for initial conditions (will be parameters for further handling
% functions
x0 = input("Initial displacement (m): ");
v0 = input("Initial velocity (m/s):   ");
clc;

response_func = 0;

% Handle variations between easy, hardcoded examples and general solution
% examples
if (has_extra_forcing == 0)
    if (num_cos_forcings == 1 && num_sin_forcings == 0)
        % Handle single cosine case
        fprintf("Single cos\n");
        return;
    elseif (num_sin_forcings == 1 && num_cos_forcings == 0)
        % Handle single sine case
        fprintf("Single sin\n");
        return;
    end
end

% Remaining possibilities are non-standard routes, employ other methods
    

response_func = other_func;

end

% Count elements of sin/cos/impulse matrices (should be pre-verified)
%   matrix - cose, sin, or impulse matrix
%   count - number of rows in matrix
function count = count_matrix(matrix)
% Set up count variable    
count = 0;

% Perform size checks and counts
[num_elems num_cols] = size(matrix);
if (num_cols == 1)  % 1 column => [0] matrix/no forcings
    return;
else
    count = num_elems;
    return;
end

end

% Checks for properly formatted sin/cos/impulse matrices
% If it is not formatted, sets matrix to [0] and prints message
% If it is formatted, sort the matrix based if it is an impulse matrix
%   matrix - cos, sin, or impulse matrix to be verified
%   type - "cos", "sin" or "impulse" to format error message
function new_matrix = check_matrix_form(matrix, type)
    
    % Get the dimensions of the matrix
    [size_x,size_y] = size(matrix);
    
    % Checks based on size
    if (size_y == 2)    % Case 1: Properly formated with values
        new_matrix = matrix;
        
        % Sort impulse matrix
        if (strcmp("impulse", type) == 1)
            new_matrix = sortrows(new_matrix, 2);
        end
        
        return;
    elseif (size_x == 1 && size_y == 1 && matrix(1,1) == 0)     % Case 2: Properly formated with [0]
        new_matrix = matrix;
        return;
    else    % Case 3: Ill formated matrix
        fprintf("The %s matrix is ill-formed, setting it to [0], restart program by ctrl+C if desired\n", type);
        new_matrix = [0];
    end
    
end