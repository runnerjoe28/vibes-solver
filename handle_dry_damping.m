% Handles the solving of dry damped systems
function success = handle_dry_damping()

success = 1; % Indicates failed, change before returning

% Ask user for mass, spring constant, normal force, and both friction coefficients
fprintf("Enter system parameters\n");
m =      input("Mass (kg):                       ");
k =      input("Spring Constant (N/m):           ");
normal = input("Normal Force (N):                ");
uk =     input("Coefficient of Kinetic Friction: ");
us =     input("Coefficient of Static Friction:  ");

% Clear command window and ask user for initial conditions
clc;
x0 = input("Enter inital displacement (m): ");

% Output EOM in variable and numeric form
clc;
fprintf("EOM for Dry Damped System: \n");
fprintf("Variable form: mx''(t) + kx(t) = -µN*sgn(x'(t))\n");
fprintf("Numeric form:  (%.3e)*x''(t) + (%.3e)*x(t) = -(%.3e)*(%.3e)*sgn(x'(t))\n", m, k, uk, normal);

% Output solved x(t) expression in variable and numeric form
friction_factor = uk*normal/k;
wn = sqrt(k/m);
fprintf("\nSolved x(t) expression\n");
fprintf("Variable form: (x_0 - (2n-1)*μN/k)*cos(ω_n t) - (-1)^n * μN/k\n");
fprintf("Numeric form:  (%.3e - (2n-1)*%.3e) * cos(%.3e * t) - (-1)^n * %.3e\n", x0, friction_factor, wn, friction_factor);

% Compute values correlating to stopping
crit_stop_distance = us*normal/k;
stop_cycles = round((x0*k - us*normal)/(2*uk*normal) - 0.5);
if (stop_cycles < 0)
    stop_cycles = 0;
end
stop_time = stop_cycles*pi()/wn;

% Print summary of new computed values
fprintf("\nCritical computed values\n");
fprintf("Natural frequency (rad/s): %.3e\n", wn);
fprintf("Stopping distance (m):     %.3e\n", crit_stop_distance);
fprintf("Stopping cycles:           %u\n", stop_cycles);
fprintf("Stopping time (s):         %.3e\n", stop_time);

% Create symbolic function to use for response (needs n and t)
syms response(n,t)
response(n,t) = (x0 - (2*n-1)*friction_factor)*cos(wn*t) - (-1)^n * friction_factor;

% Generate graph at each piece of the piecewise response function
hold on;
for n = 1:(stop_cycles + 1)
    t1 = (n-1)*pi()/wn;
    t2 = n*pi()/wn;
    t_plot = linspace(t1, t2, 10);
    x_plot = response(n, t_plot);
    if (n >= stop_cycles + 1)
        plot(t_plot, x_plot, 'b--');
    else
        plot(t_plot, x_plot, 'b');
    end
end

% Generate threshold lines indicating termination of vibration
yline(crit_stop_distance, 'r', "+ Threshold");
yline(-crit_stop_distance, 'r', "- Threshold");

% Set title and axes
title("Dry damping response with terminator lines");
xlabel("Time (s)");
ylabel("Displacement (m)");
hold off;

success = 0;

end