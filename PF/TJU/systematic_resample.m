function [x_new, w_new] = systematic_resample(x, w, Ns)
% Rensselaer Polytechnic Institute - Julius Lab
% SenSE Project
% Author - Chukwuemeka Osaretin Ike
%
% Description:
% Systematic resampling for a particle filter.
    Q = cumsum(w);

    % Random starting point.
    u0 = rand/Ns;

    x_new = zeros(size(x));
    w_new = (1/Ns)*ones(size(w));
    m = 1;
    for idx_resamp = 1:Ns
        u_i = u0 + (idx_resamp-1)/Ns;
        while u_i > Q(m)
            m = m+1;
        end
        x_new(idx_resamp,:) = x(m,:);
    end
end