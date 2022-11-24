% Performs a fourier transform on the specified function
function success = handle_fourier()

% Set up variables for use in computation
success = 1; % Indicates failed, change before returning
syms t

% Ask for critical information from user
fprintf("Enter problem parameters\n");
func =     input("f(t) whose Fourier series is to be found:    ");
i_bounds = input("Bounds of one period of oscillation [a , b]: ");
num_h =    input("Number of harmonics:                         ");

% Bounds of integration
lower_bound = i_bounds(1);
upper_bound = i_bounds(2);

% Set up key parts of fourier function
T = upper_bound - lower_bound;
a0 = (2/T)*(int(func, lower_bound, upper_bound));
fourier_func = a0/2; % Fourier series function of key interest

% Set up an and bn vectors
an = zeros(1, num_h);
bn = zeros(1, num_h);

% Calculate the first num_h harmonics
for n=1:num_h
    
    % Compute a_n and b_n for this series
    omegan = 2*pi*n/T;
	an(n)=(2/T) * (int(func * cos(omegan * t), lower_bound, upper_bound));
	bn(n)=(2/T) * (int(func * sin(omegan * t), lower_bound, upper_bound));
    
    % Update fourier function
	fourier_func = fourier_func + an(n)*cos(omegan * t) + bn(n)*sin(omegan * t);
			
end

% Print all associated information with the transform
clc;
disp("Information associated with the fourier series transform")
disp(strcat("Original function: ", char(func)));
disp(strcat("Fourier series:    ", char(fourier_func)));
fprintf("Value of a0:       %.3e\n", a0);
fprintf("Values of a_n:     ");
fprintf("%d ", an);
fprintf("\nValues of b_n:     ");
fprintf("%d ", bn);
fprintf("\n");

% Create a graphical representation of the original and fourier series
hold on;
ezplot(func, [lower_bound,upper_bound]);
ezplot(fourier_func, [lower_bound,upper_bound]);
title([num_h, " Harmonics"]);
legend("Original Function", "Fourier Series");
hold off;
end

