% Handles the solving of standard damped systems
function success = handle_standard_damping()

success = 1; % Indicates failed, change before returning

% Get system parameters (mass, damping coefficient, and spring constant)
fprintf("Enter system parameters (enter 0 for damping if no damping is present)\n");
m = input("Mass (kg):                  ");
c = input("Damping coefficient (kg/s): ");
k = input("Spring constant (N/m):      ");

% Create response_func for general plotting later
syms response_func(t)
response_func = 1;

% Ask if there are forcing functions and handle the NO case
is_forcing = input("\nIs there a forcing function [Y/N]: ", 's');
if (is_forcing == 'N')
    response_func = solve_EOM(m, c, k);
    success = 0;
else
    % handle later
end

% Plot the response function
fplot(response_func);
upper_limit = 6*pi*sqrt(k/m)/(1+c);
xlim([0 upper_limit]);

end


% Solves basic no/damped system without forcing functions
function response_func = solve_EOM(m, c, k)

% Set up conditional for determing damping classification + key variables
zeta = c/(2*sqrt(k*m));
wn = sqrt(k/m);
syms t;

% Ask for initial conditions
clc;
fprintf("Enter initial conditions\n");
x0 = input("x0: ");
v0 = input("v0: ");
clc;

% Handle all 4 types of system
if (c == 0) % No damping
    
    % Compute response
    amp = sqrt(wn^2 * x0^2 + v0^2)/wn;
    phi = atan(wn*x0/v0);
    response_func = amp*sin(wn*t + phi);
    
    % Outpute EOM in variable and numeric form
    fprintf("EOM for undamped, unforced system: \n");
    fprintf("Variable form: mx''(t) + kx(t) = 0\n");
    fprintf("Numeric form:  (%d)*x''(t) + (%d)*x(t) = 0\n", m, k);

    % Output solved x(t) expression in variable and numeric form
    fprintf("\nSolved x(t) expression\n");
    fprintf("Variable form: A*sin(wn*t + Φ)\n");
    fprintf("Numeric form:  %d*sin(%d*t + %d)\n", amp, wn, phi);

    % Display critical system parameters/values
    fprintf("\nKey system parameters for undamped, unforced response\n");
    fprintf("Mass (kg):                 %d\n", m);
    fprintf("Spring constant (N/m):     %d\n", k);
    fprintf("Natural frequency (rad/s): %d\n", wn);
    fprintf("Intial displacement (m):   %d\n", x0);
    fprintf("Intial velocity (m/s):     %d\n", v0);
    fprintf("Response amplitude (m):    %d\n", amp);
    fprintf("Response offset (rad):     %d\n\n", phi);
    
elseif (zeta < 1) % Under damped
    
    % Compute response
    wd = wn*sqrt(1 - zeta^2);
    amp = sqrt( (v0 + zeta*wn*x0)^2 + (x0*wd)^2 )/wd;
    phi = atan( (x0*wd) / (v0 + zeta*wn*x0) );
    response_func = amp * exp(-zeta*wn*t) * sin(wd*t + phi);
    
    % Outpute EOM in variable and numeric form
    fprintf("EOM for underdamped, unforced system: \n");
    fprintf("Variable form: mx''(t) + cx'(t) kx(t) = 0\n");
    fprintf("Numeric form:  (%d)*x''(t) + (%d)*x'(t) + (%d)*x(t) = 0\n", m, c, k);

    % Output solved x(t) expression in variable and numeric form
    fprintf("\nSolved x(t) expression\n");
    fprintf("Variable form: A * exp(-ζ*wn*t) * sin(wn*t + Φ)\n");
    fprintf("Numeric form:  %d * exp(%d*t) * sin(%d*t + %d)\n", amp, -zeta*wn, wn, phi);

    % Display critical system parameters/values
    fprintf("\nKey system parameters for underdamped, unforced response\n");
    fprintf("Mass (kg):                        %d\n", m);
    fprintf("Damping coefficient (kg/s):       %d\n", c);
    fprintf("Spring constant (N/m):            %d\n", k);
    fprintf("Natural frequency (rad/s):        %d\n", wn);
    fprintf("Damped natural frequency (rad/s): %d\n", wd);
    fprintf("Damping ratio:                    %d\n", zeta);
    fprintf("Intial displacement (m):          %d\n", x0);
    fprintf("Intial velocity (m/s):            %d\n", v0);
    
elseif (zeta == 1) % Critically damped

    % Compute response
    a1 = x0;
    a2 = v0 + wn*x0;
    response_func = (a1 + a2*t) * exp(-wn*t);
    
    % Outpute EOM in variable and numeric form
    fprintf("EOM for critcally damped, unforced system: \n");
    fprintf("Variable form: mx''(t) + cx'(t) kx(t) = 0\n");
    fprintf("Numeric form:  (%d)*x''(t) + (%d)*x'(t) + (%d)*x(t) = 0\n", m, c, k);

    % Output solved x(t) expression in variable and numeric form
    fprintf("\nSolved x(t) expression\n");
    fprintf("Variable form: (a1 + a2*t) * exp(-wn*t)\n");
    fprintf("Numeric form:  (%d + %d*t) * exp(%d*t)\n", a1, a2, -wn);

    % Display critical system parameters/(values
    fprintf("\nKey system parameters for critcally damped, unforced response\n");
    fprintf("Mass (kg):                  %d\n", m);
    fprintf("Damping coefficient (kg/s): %d\n", c);
    fprintf("Spring constant (N/m):      %d\n", k);
    fprintf("Natural frequency (rad/s):  %d\n", wn);
    fprintf("Intial displacement (m):    %d\n", x0);
    fprintf("Intial velocity (m/s):      %d\n", v0);
    
else % Over damped
    
    % Compute response
    zeta_term = sqrt(zeta^2 - 1);
    a1 = (-v0 + (-zeta+zeta_term)*wn*x0)/(2*wn*zeta_term);
    a2 = (v0 + (zeta+zeta_term)*wn*x0)/(2*wn*zeta_term);
    response_func = exp(-zeta*wn*t)*(a1*exp(-wn*zeta_term*t) + a2*exp(wn*zeta_term*t));
    
    % Outpute EOM in variable and numeric form
    fprintf("EOM for overdamped, unforced system: \n");
    fprintf("Variable form: mx''(t) + cx'(t) kx(t) = 0\n");
    fprintf("Numeric form:  (%d)*x''(t) + (%d)*x'(t) + (%d)*x(t) = 0\n", m, c, k);

    % Output solved x(t) expression in variable and numeric form
    fprintf("\nSolved x(t) expression\n");
    fprintf("Variable form: exp(-zeta*wn*t) * ( a1*exp(-wn*sqrt(zeta^2 - 1)*t) + a2*exp(wn*sqrt(zeta^2 - 1)*t) )\n");
    fprintf("Numeric form:  exp(%d*t) * ( %d*exp(%d*t) + %d*exp(%d*t) )\n", -zeta*wn, a1, -wn*zeta_term, a2, wn*zeta_term);

    % Display critical system parameters/(values
    fprintf("\nKey system parameters for overdamped, unforced response\n");
    fprintf("Mass (kg):                  %d\n", m);
    fprintf("Damping coefficient (kg/s): %d\n", c);
    fprintf("Spring constant (N/m):      %d\n", k);
    fprintf("Natural frequency (rad/s):  %d\n", wn);
    fprintf("Intial displacement (m):    %d\n", x0);
    fprintf("Intial velocity (m/s):      %d\n", v0);
    
end

end

