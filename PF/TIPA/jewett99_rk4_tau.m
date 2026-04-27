function [t, u, y] = jewett99_rk4_tau(time, input, initial_conditions, tau)
% Description:
% Fixed-step implementation of the ODEs from the paper below. I've found
% that the ODE solvers from MATLAB are all variable step, and are strictly
% slower than this implementation.
%
% Takes tau as an argument for use in the particle filtering algorithm.
%
% Paper:
% Revised limit cycle oscillator model of human circadian pacemaker -
% Jewett (1999).

% Process L Parameters.
a0 = 0.16;
p = 0.6;
beta = 0.013;
G = 19.875;
I0 = 9500;

% Process P Parameters.
mu = 0.13;
k = 0.55;
q = 1/3;
kc = 0.4;
% tau = 24.2;

t = time;
dt = t(2) - t(1);
y = zeros(length(time), length(initial_conditions));
y(1, :) = initial_conditions;
for idx = 2:length(t)
    y(idx, :) = jfk_forward(y(idx-1, :), input(idx), dt);
end

u = G*(a0*(input/I0).^p).*(1-y(:,3));

%%
function xk_1 = jfk_forward(xk, Ik, h)  
x_ = xk(1);
xc = xk(2);
n = xk(3);

a = a0*(Ik/I0)^p;

k1 = f(x_, xc, n, a);    
y1 = xk + (k1*h/2);

k2 = f(y1(1), y1(2), y1(3), a);
y2 = xk + (k2*h/2);

k3 = f(y2(1), y2(2), y2(3), a);
y3 = xk + (k3*h);

k4 = f(y3(1), y3(2), y3(3), a);

xk_1 = xk + (1/6)*( k1 + (2*k2) + (2*k3) + k4 )*h;   
end

function xdot = f(x, xc, n, a)
    uk = G*a*(1 - n);
    B = (1 - 0.4*x)*(1 - kc*xc)*uk;

    xdot = [(pi/12)*( xc + mu*( (x/3) + (4*x^3/3) - (256*x^7/105) ) + B ),...
        (pi/12)*( q*xc*B - x*( (24/0.99729/tau)^2 + k*B ) ),...
        60*( (a*(1-n)) - beta*n )];
end
end