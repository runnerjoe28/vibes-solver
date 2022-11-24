% Clear workspace and command window
clc;
clear;

% Initiate options for vibrations solving, and display to user in prompt
% There are 5 options for users
%       1. Standard damping (none, under, critically, and over dampled)
%       2. Dry damping
%       3. Based Excitation
%       4. Rotating Unbalance
%       5. Fourier Transform
fprintf("Problem types available to solve\n");
fprintf("\t1 - Standard damping (none, under, critically, and over dampled)\n");
fprintf("\t2 - Dry damping\n");
fprintf("\t3 - Based Excitation\n");
fprintf("\t4 - Rotating Unbalance\n");
fprintf("\t5 - Fourier Transform\n");
problem_type_num = input("\nPlease enter number of problem type: ");
fprintf("%i\n", problem_type_num);

% Handle each problem type (delegate to other .m files)
clc;
switch problem_type_num
    case 1
        handle_standard_damping();
    case 2
        handle_dry_damping();
    case 3
    case 4
    case 5
        handle_fourier();
    otherwise
        fprintf("INVALID INPUT, TERMINATING PROGRAM\n");
end


