
function response_func = handle_rotating_imbalance(x, w, m, X, w2)
clear;
clc;
fprintf("Enter system parameters (enter 0 for damping if no damping is present)\n");
m = input("Mass (kg):                  ");
x = input("Spring deflection m:      ");
w = input("device operating at (rad/s):      ");
X = input("Amplitude of oscillations (m):      ");
w2 = input("device will be operating at (rad/s):    ");
if m == 0
    wn = sqrt(9.81/x);
    eqtn = X/((w^2)/(1-(w^2)/wn^2));
    X2 = eqtn*((w2^2)/(1-((w2^2)/(wn^2))))
else
    X1 = input("Orignial Amplitude of oscillations (m):      ");
    wn = sqrt(9.81/x)
    eqtn = X1/((w^2)/(1-(w^2)/wn^2))
    r = w2/wn;
    z = (sqrt(-w^4*X^2+2*w^2*X^2*wn^2+eqtn^2*w^4*wn^4-X^2*wn^4))/(2*w*X*wn)
    c = z*(2*m*wn)
    X2 = eqtn * ((w2^2) / (sqrt((1 - r)^2 + (2 * z * r)^2)))    
end
end
    
    


   



