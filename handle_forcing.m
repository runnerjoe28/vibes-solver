% Essentially a library for solving forcing functions

% Create the response for a forced single DOF function
%   m - mass of the system
%   c - damping coefficient
%   k - spring constant
%   cos_matrix - N x 2 matrix with amplitude,frequency pairs
%   sin_matrix - N x 2 matrix with amplitude,frequency pairs
%   imp_matrix - N x 2 matrix with amplitude,time pairs
%   other_func - symbolic function other other forcings
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
has_extra_forcing = 1;  % 1 if other func is non-zero, 0 if the other func is 0
if (other_func == 0)
    has_extra_forcing = 0;
end

% Ask for initial conditions (will be parameters for further handling
% functions
x0 = input("Initial displacement (m): ");
v0 = input("Initial velocity (m/s):   ");
clc;

response_func = 0;

% Get the matrix of impulse reactions
imp_responses = find_impulse_response(imp_matrix, m, c, k);
disp(imp_responses)

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
[num_elems, num_cols] = size(matrix);
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
        new_matrix = 0;
    end
    
end

% Computes a response of all impulses
% If there are no impulses (input is 0) then the response is zero
%   impulse: N X 2 properly formatted impulse matrix
%   imp_resp: response function from the impulses
function imp_resp = find_impulse_response(impulses, m, c, k)

% Set up variables and get number of impulses
imp_resp = 0;
[num_impulses, num_cols] = size(impulses);

% Handle no impulse case
if (num_cols == 1)
    return;
end

% Compute key parameters
zeta = c/(2*sqrt(k*m));
wn = sqrt(k/m);
wd = wn*sqrt(1 - zeta^2);

% Set up return statements
syms t;

% Loop through impulses, construct response
for i = 1:num_impulses

    % Get impulse information here
    impulse_mag = impulses(i,1);
    impulse_time = impulses(i,2);
    fprintf("Impulse of %i at t=%i\n", impulse_mag, impulse_time);
    cur_impulse_function = 0;
    
    % Set initial position and initial velocity from impulse
    x0 = 0;
    v0 = impulse_mag/m;
    
    % Handle all 4 types of system
    if (c == 0) % No damping

        % Compute response
        amp = sqrt(wn^2 * x0^2 + v0^2)/wn;
        phi = atan(wn*x0/v0);
        cur_impulse_function = amp*sin(wn*(t - impulse_time) + phi);

    elseif (zeta < 1) % Under damped

        % Compute response
        amp = sqrt( (v0 + zeta*wn*x0)^2 + (x0*wd)^2 )/wd;
        phi = atan( (x0*wd) / (v0 + zeta*wn*x0) );
        cur_impulse_function = amp * exp(-zeta*wn*(t - impulse_time)) * sin(wd*(t - impulse_time) + phi);

    elseif (zeta == 1) % Critically damped

        % Compute response
        a1 = x0;
        a2 = v0 + wn*x0;
        cur_impulse_function = (a1 + a2*(t - impulse_time)) * exp(-wn*(t - impulse_time));

    else % Over damped

        % Compute response
        zeta_term = sqrt(zeta^2 - 1);
        a1 = (-v0 + (-zeta+zeta_term)*wn*x0)/(2*wn*zeta_term);
        a2 = (v0 + (zeta+zeta_term)*wn*x0)/(2*wn*zeta_term);
        cur_impulse_function = exp(-zeta*wn*(t - impulse_time))*(a1*exp(-wn*zeta_term*(t - impulse_time)) + a2*exp(wn*zeta_term*(t - impulse_time)));

    end
    
    % Add current function to total response
    imp_resp = imp_resp + heaviside(t - impulse_time)*cur_impulse_function;
    
end

end


