function run_particle_filter_delta_tau_ks_phi_compressed(subject, run, test, Q, R, I, t, option)
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
if option == 2
    option_info = "Act_";
elseif option == 3
    option_info = "BPM_";
elseif option == 4
    option_info = "AVGIBI_";
elseif option == 5
    option_info = "MinIBI_";
elseif option == 6
    option_info = "MaxIBI_";
elseif option == 7
    option_info = "SD_";
elseif option == 8
    option_info = "RMSSD_";
elseif option == 9
    option_info = "LF_";
elseif option == 10
    option_info = "HF_";
elseif option == 11
    option_info = "LFHF_";
elseif option == 12
    option_info = "Wrist_Act_";
end

load(strcat("A", num2str(subject), "_KS_",option_info, num2str(run),".mat"), "phi")
load(strcat("MatData\NSF0", num2str(subject), "_MTN.mat"), "init")

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

% seed = test + 1;
% rng(seed, 'twister');

xHat = zeros(Ns, numStates, length(y));

% Uniformly distributed initial states.
xHat(:, :, 1) = (x_min + (x_max - x_min).*rand(numStates, Ns))';

% % Normally distributed x and xc, and uniform in n and delta.
% stdev = 0.5;
% xHat(:,1:2,1) = normrnd(repmat(init(1:2), [Ns, 1]), stdev);
% % xHat(:,3:4,1) = (x_min(3:4) + (x_max(3:4) - x_min(3:4)).*rand(2, Ns))';

% Normally distributed tau based on Duffy (24h09m +- 12m).
xHat(:, 5, 1) = normrnd(repmat(24.15, [Ns, 1]), 0.2);

% Set the Huang initial condition for the first particle.
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

    % Compute the filter weighted mean for linear states (x, xc, n, tau)
    x_mean(k, [1, 2, 3, 5]) = weights(:, k)' * squeeze(xHat(:, [1, 2, 3, 5], k));
    
    % Projects the angles off the 1D circular boundary into 2D Cartesian space
    % Compute the circular weighted mean for angular state (Delta)
    sin_val = weights(:, k)' * sin(xHat(:, 4, k));
    cos_val = weights(:, k)' * cos(xHat(:, 4, k));
    x_mean(k, 4) = wrapTo2Pi(atan2(sin_val, cos_val));

    % Resample.
    Neff = 1/sum(weights(:,k).^2);
    if Neff < N_thresh
        % Systematic resample alone.
        [xHat(:,:,k), weights(:,k)] = systematic_resample(xHat(:,:,k), weights(:,k), Ns);

        % Add jitter to resampled particles.
        xHat(:,:,k) = xHat(:,:,k) + normrnd(zeros(Ns,numStates), repmat(Q(1:5), [Ns, 1]));

        % Maintain tau within the selected bounds.
        xHat(:,5,k) = max(x_min(5), min(x_max(5), xHat(:,5,k)));

        % fprintf("Resample at iter: %d\n", k);
    end

    % Progress tracking.
    if k < 10
        % fprintf("Step %d runtime: %3.4f s\n", k, seconds(datetime-step_start));
    elseif mod(k, 1500) == 0
        % fprintf("Step %d runtime: %3.4f s\n", k, seconds(datetime-epoch_start));
        epoch_start = datetime;
    end
end

theta = atan2(x_mean(:,1), x_mean(:,2));
phi = compute_phi(t, theta, omg);
y_mean = wrapTo2Pi(phi + x_mean(:,4));

pf_end = datetime;
fprintf("Total runtime: %s\n", pf_end - pf_start)


%% DLMO and Standard Deviation Calculation
% Based on the paper: DLMO is approximately CBT_min - 7 hours.
% CBT_min occurs at the local minimum of the circadian state x (State 1).

% --- Day 7 DLMO ---
% Isolate the time window for Day 7 (hours 144 to 168)
day7_idx = find(t >= 144 & t < 168);
if ~isempty(day7_idx)
    [~, min_loc_7] = min(x_mean(day7_idx, 1)); % Find CBT_min index
    idx_dlmo_7 = day7_idx(min_loc_7);         % Exact array index
    
    cbt_min_time_7 = t(idx_dlmo_7);
    dlmo_hour_7 = mod(cbt_min_time_7 - 7, 24); % Calculate DLMO and map to 0-24h
    
    % Get standard deviation of all 5 states across the 500 particles at this exact moment
    std_dlmo_7 = std(xHat(:, :, idx_dlmo_7), 0, 1); 
else
    dlmo_hour_7 = NaN;
    std_dlmo_7 = nan(1, numStates);
end

% --- Day 14 DLMO ---
% Isolate the time window for Day 14 (hours 312 to 336)
day14_idx = find(t >= 312 & t < 336);
if ~isempty(day14_idx)
    [~, min_loc_14] = min(x_mean(day14_idx, 1)); % Find CBT_min index
    idx_dlmo_14 = day14_idx(min_loc_14);         % Exact array index
    
    cbt_min_time_14 = t(idx_dlmo_14);
    dlmo_hour_14 = mod(cbt_min_time_14 - 7, 24); % Calculate DLMO and map to 0-24h
    
    % Get standard deviation of all 5 states across the 500 particles at this exact moment
    std_dlmo_14 = std(xHat(:, :, idx_dlmo_14), 0, 1);
else
    dlmo_hour_14 = NaN;
    std_dlmo_14 = nan(1, numStates);
end

estimate_dlmo = [dlmo_hour_7, dlmo_hour_14];

dlmo_hour = [];
t_days = t / 24; % Use a temporary variable so the original 't' is preserved in hours for saving
for j = 1:round(t_days(end))
    day_idx = find(round(t_days) == j);
    if ~isempty(day_idx)
        day_t = t_days(day_idx);
        day_x = x_mean(day_idx, 1);
        [~, min_idx] = min(day_x);
        dlmo_hour = [dlmo_hour, day_t(min_idx) * 24];
    end
end
dlmo_hour = mod(dlmo_hour - 7, 24);


%% Compress and Downsample for storage
downsample_step = 10; % Change this to 5 or 10 as needed

% Downsample the massive particle arrays by taking every Nth time step
xHat = xHat(:, :, 1:downsample_step:end);
weights = weights(:, 1:downsample_step:end);
% Optional: You can also downsample yHat if it is taking up too much space
% yHat = yHat(:, 1:downsample_step:end);

% Compress to single precision to cut the remaining storage in half
xHat = single(xHat);
weights = single(weights);

%% Save all data to .mat file
save(strcat("E:\Particle Filter Data Recreate\A", num2str(subject), "_PF_Out_Delta_Tau_",option_info, num2str(run+test), "_phi.mat"), ...
    "xHat", "x_mean", "yHat", "y_mean", "weights", "Q", "R", "measure_start")

save(strcat('E:\Particle Filter Data Recreate\A',num2str(subject),'min_CBT_time',option_info, num2str(run+test), "_phi.mat"), "dlmo_hour", "estimate_dlmo", "std_dlmo_7", "std_dlmo_14")
end