function run_particle_filter_delta_tau(subject, run, Q, R, I, t)
% Rensselaer Polytechnic Institute - Julius Lab
% SenSE Project
% Original Author - Chukwuemeka Osaretin Ike
% Modified by Zidi Tao
%
% LCO Model + PF + Delta + Tau.
%
% Description:
% Implements the LCO model state estimation with a PF. Augments the "state"
% with delta bias in phi and tau in the model.


%% Load subject data.
load(strcat("Data\A", num2str(subject), "_KF_", num2str(run),".mat"), "phi")
load(strcat("Data\A", num2str(subject), "_MTN.mat"), "init")

y = wrapTo2Pi(movmean(unwrap(phi), 1440));

% Time to start incorporating measurements. Greater than t(1) because
% there's about 2 days of transients in KF output.
measure_start = 48;
measure_idx = find(t >= measure_start, 1, "first");

omg = 2*pi/24;

%% Particle filtering.

% Start timing.
pf_start = datetime;

% PF Hyperparameters.
Ns = 500;            % Number of particles.
N_thresh = Ns/4;    % Threshold for resampling.

% Expected model covariance parameters.
% Q = [5e-3, 5e-3, 1e-3, 1e-2, 1e-2];  % Roughening when we resample.
% Q = [1e-3, 1e-3, 0, 1e-4];      % Direct roughening.
% R = 0.85;

% Initial state distribution.
% x, xc, n, Delta, tau.
numStates = 5;
x_max = [1.5; 1.5; 1; 2*pi; 24.6];
x_min = [-1.5; -1.5; 0; 0; 23.8];

xHat = zeros(Ns, numStates, length(y));

% Uniformly distributed initial states.
xHat(:, :, 1) = (x_min + (x_max - x_min).*rand(numStates, Ns))';

% % Normally distributed x and xc, and uniform in n and delta.
% stdev = 0.5;
% xHat(:,1:2,1) = normrnd(repmat(init(1:2), [Ns, 1]), stdev);
% % xHat(:,3:4,1) = (x_min(3:4) + (x_max(3:4) - x_min(3:4)).*rand(2, Ns))';

% Normally distributed tau based on Duffy (24h09m +- 12m).
xHat(:, 5, 1) = normrnd(repmat(24.15, [Ns, 1]), 0.2);

% Set the Huang initial condition.
xHat(1, 1:3, 1) = init;
xHat(1, 4, 1) = pi;
% xHat(1, 5, 1) = 24.2;

% Filter estimate (weighted mean of all particles).
x_mean = zeros(length(y), numStates);
x_mean(1, :) = mean(xHat(:, :, 1));

% Output estimates.
yHat = zeros(Ns, length(y));
theta = atan2(x_mean(1,1), xHat(1,2));
phi = compute_phi(t(1), theta, omg);
yHat(:,1) = phi + x_mean(1,4);

% Particle weights.
weights = zeros(Ns, length(y));
weights(:,1) = 1/Ns;

% Start the filter.
epoch_start = datetime;
for k = 2:length(y)
    % Start timing the iteration.
    step_start = datetime;
    for i = 1:Ns
        % Propagate the particle forward.
        [~,~,tx] = jewett99_rk4_tau(t(k-1:k), I(k-1:k), squeeze(xHat(i, 1:3, k-1))', xHat(i, 5, k-1));

        % Add some jitter to the states.
        xHat(i, 1:3, k) = tx(end,:);
        xHat(i, 4, k) = wrapTo2Pi(xHat(i, 4, k-1));
        xHat(i, 5, k) = xHat(i, 5, k-1);

        % Compute output due to current state.
        theta = atan2(xHat(i, 1, k), xHat(i, 2, k));
        phi = compute_phi(t(k), theta, omg);
        yHat(i, k) = wrapTo2Pi(phi + xHat(i, 4, k));

        % Compute weight update if measurement available.
        if k > measure_idx
            % Difference on a circle.
            circle_diff = wrapToPi(yHat(i,k) - y(k));
            weights(i, k) = weights(i, k-1)*normpdf(circle_diff, 0, R);
        else
            weights(i, k) = weights(i, k-1);
        end
    end

    % Normalize the weights.
    weights(:,k) = weights(:,k)/sum(weights(:,k));

    % Compute the filter weighted mean.
    x_mean(k, :) = weights(:, k)'*squeeze(xHat(:, :, k));

    % Resample.
    Neff = 1/sum(weights(:,k).^2);
    if Neff < N_thresh
        % % Resample-move.
        % [xHat(:,:,k), weights(:,k)] = metropolis_hastings(t(k), omg, xHat(:,:,k), weights(:,k), y(k), Ns, R);

        % Systematic resample alone.
        [xHat(:,:,k), weights(:,k)] = systematic_resample(xHat(:,:,k), weights(:,k), Ns);

        % Add jitter to resampled particles.
        xHat(:,:,k) = xHat(:,:,k) + normrnd(zeros(Ns,numStates), repmat(Q(1:5), [Ns, 1]));

        % Maintain tau within the selected bounds.
        xHat(:,5,k) = max(x_min(5), min(x_max(5), xHat(:,5,k)));
    end

    % Progress tracking.
    if k < 10
        fprintf("Step %d runtime: %3.4f s\n", k, seconds(datetime-step_start))
    elseif mod(k, 1500) == 0
        fprintf("Step %d runtime: %3.4f s\n", k, seconds(datetime-epoch_start))
        epoch_start = datetime;
    end
end

theta = atan2(x_mean(:,1), x_mean(:,2));
phi = compute_phi(t, theta, omg);
y_mean = phi + x_mean(:,4);

pf_end = datetime;
fprintf("Total runtime: %s\n", pf_end - pf_start)


%%
save(strcat("Results\A", num2str(subject), "_PF_Out_Delta_Tau_", num2str(run), ".mat"), ...
    "xHat", "x_mean", "yHat", "y_mean", "weights", "Q", "R", "measure_start")
end