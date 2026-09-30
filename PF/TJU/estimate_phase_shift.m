function [theta, phase_shift] = estimate_phase_shift(day1, day2, omg, t, x1, x2)
% Rensselaer Polytechnic Institute - Julius Lab
% SenSE Project
% Author - Chukwuemeka Osaretin Ike
%
% Description:
% Estimates the phase shift given xHat from the KF/OBF.

    % Calculate the phase based on x1 and x2
    theta = wrapTo2Pi(-atan2(x2, omg*x1));
    % theta = mod(-atan2(x2, omg*x1) + pi/2, 2*pi)-pi;
    t = t - t(1);

    % Calculate the phase shift between day 1 and 2
    day1Unwrap = unwrap(theta((t >= ((day1-1)*24)) & (t <= day1*24)));
    day2Unwrap = unwrap(theta((t >= ((day2-1)*24)) & (t <= day2*24)));
    l1 = length(day1Unwrap);
    l2 = length(day2Unwrap);

    if l1 >= l2
        phase_shift = (1/omg)*mean(day1Unwrap(1:l2) - day2Unwrap);
    else
        phase_shift = (1/omg)*mean(day1Unwrap - day2Unwrap(1:l1));
    end
    phase_shift = mod(12+phase_shift, 24) - 12;
end