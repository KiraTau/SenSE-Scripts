function L = activity_to_light(a)
% Rensselaer Polytechnic Institute - Julius Lab
% SenSE Project
% Author - Chukwuemeka Osaretin Ike
%
% Description:
% Converts an activity signal to light exposure using the method from:
% "Predicting circadian phase across populations: a comparison of
% mathematical models and wearable devices" - Huang (2021).
%
% Parameters:
%   a - vector of step counts.
% Returns:
%   L - vector of estimated light exposure.

    L = zeros(size(a));
    m = max(a)/2;
    for i = 1:length(a)
        if a(i) <= 0
            L(i) = 0;
        elseif a(i) > 0 && a(i) < 0.1*m
            L(i) = 100;
        elseif a(i) >= 0.1*m && a(i) < 0.25*m
            L(i) = 200;
        elseif a(i) >= 0.25*m && a(i) < 0.4*m
            L(i) = 500;
        else
            L(i) = 2000;
        end
    end
end