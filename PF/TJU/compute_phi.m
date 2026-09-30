function phi = compute_phi(t, theta, omg)
% Rensselaer Polytechnic Institute - Julius Lab
% SenSE Project
% Author - Chukwuemeka Osaretin Ike
%
% Description:
% Computes phi given t, theta, and omega.
% Returns phi in [0, 2*pi]. Expects theta in radians in same interval.

theta = reshape(theta, size(t));
phi = unwrap(theta) - omg*t;
phi = wrapTo2Pi(phi);
end