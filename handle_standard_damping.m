% Handles the solving of standard damped systems
function success = handle_standard_damping()

success = 1; % Indicates failed, change before returning

% Get system parameters (mass, damping coefficient, and spring constant)
fprintf("Enter system parameters (enter 0 for damping if no damping is present)\n");
m = input("Mass (kg):                  ");
c = input("Damping coefficient (kg/s): ");
k = input("Spring constant (N/m):      ");

% Ask if there are forcing functions and handle the NO case
is_forcing = input("\nIs there a forcing function [Y/N]: ", 's');
if (is_forcing == 'N')
    success = solve_EOM(m, c, k);
    return;
end

end


% Solves basic no/damped system without forcing functions
function success = solve_EOM(m, c, k)

% Set up conditional for determing damping classification
success = 1; % Indicates failed, change before returning
c_cr = c/(2*sqrt(k*m));

if (c ==0) % No damping
    fprintf("No damping\n");
elseif (c_cr < 1) % Under damped
    fprintf("Under damped\n");
elseif (c_cr == 1) % Critically damped
    fprintf("Critically damped\n");
else % Over damped
    fprintf("Over damped\n");
end

end

