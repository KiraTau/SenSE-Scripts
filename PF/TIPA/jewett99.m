function [t, u, y] = jewett99(time, input, initial_conditions)
% Paper:
% Revised limit cycle oscillator model of human circadian pacemaker -
% Jewett (1999).

% Process L Parameters.
a0 = 0.16;
p = 0.6;
beta = 0.013;
G = 19.875;

% I have no clue why, but Serkh (2014) and Jiawei Yin (2020) both use these
% parameters. Jiawei cites Serkh, who cites Jewett (1999), but her paper
% uses the parameters above. I have these here for clarity I guess.
% Running simulations with both obviously yields different qualitative
% behavior.
% a0 = 0.05;
% p = 0.5;
% beta = 0.0075;
% G = 33.75;

I0 = 9500;

% Process P Parameters.
mu = 0.13;
k = 0.55;
q = 1/3;
kc = 0.4;
tau = 24.2;

% options = odeset("InitialStep", time(2)-time(1));
[t, y] = ode23(@(t,y) dyn(t,y,time,input),...
               time,...
               initial_conditions(:));

% Calculate u after the fact.
if length(input) == length(t)
    u = G*(a0*(input/I0).^p).*(1-y(:,3));
else
    u = G*(a0*(input/I0).^p).*(1-y([1,end],3));
end

function xdot = dyn(t, x, t_rad, I_vec)
    I = interp1(t_rad, I_vec, t);

    % States.
    x_ = x(1);
    xc = x(2);
    n = x(3);

    xdot = zeros(3,1);

    % Process L.
    alpha = (a0*(I/I0)^p);
    xdot(3) = 60*((alpha*(1-n)) - (beta*n));
    B_hat = G*alpha*(1-n);

    % Process P.
    B = B_hat*(1 - 0.4*x_)*(1 - kc*xc);
    xdot(1) = (pi/12)*(xc + mu*((x_/3) + (4*x_^3/3) - (256*x_^7/105)) + B);
    xdot(2) = (pi/12)*(q*xc*B - x_*( (24/0.99729/tau)^2 + k*B )); 
end
end