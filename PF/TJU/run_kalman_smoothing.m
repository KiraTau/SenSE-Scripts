function run_kalman_smoothing(t_, y_, run, subject, option, seed)
% Rensselaer Polytechnic Institute - Julius Lab
% SenSE Project
% Author - Zidi Tao
%
% Description:
% Kalman Smoothing optimization and simulation. Structured as a callable 
% function to find optimal parameters, simulate the smoother, and save 
% the resulting states and parameters to a .mat file.

% Optional random seed for reproducibility (default: 12345)
if nargin < 6 || isempty(seed)
    seed = 12345;
end
rng(seed, 'twister');

%% Parameters.

% Optimization Hyperparameters.
iterations = 150;           % Maximum iterations.
mu = 100;                   % Number of every generation.
lambda = 50;                % Number of children.
rho = 2;                    % Number of parents to create one child.

% Optimization bounds (Adjusted for normalized y data).
LB = -5;                % Q Lower bound - 10^LB.
lEnd = 0;               % Q lower end.
rEnd = 18;              % Twice the number of orders of magnitudes.
rLA = 1e2;              % R Lower bound.
rLB = 1e8;              % R upper bound. 


% Filter params.
order = 1;
fprintf("Smoother Order %d\n", order)
numStates = (2*order + 1);        % State length based on order
qSize = numStates^2;
omg = 2*pi/24;

% Load and orient data.
y = y_';
t = t_';

% Handle missing data and normalize
y(y == 0) = NaN;
mu_y = mean(y, 'omitnan');
sigma_y = std(y, 'omitnan');
y_norm = (y - mu_y) / sigma_y;

%% Optimize the Smoother.
dt = t(2) - t(1); % Get the sampling time for the subject.
[A, ~, C, ~] = KS.createKalmanStateSpace(order, numStates, dt);

startT = datetime;                  % Start timing the optimization.

% Run optimization using the normalized data
[Cost, avgCost, Q_pop, R_pop] = KS.optimizeParameters(...
            iterations, mu, lambda, rho, order, lEnd, rEnd, LB, rLA, rLB,...
            A, C, t, y_norm);
            
runTime = datetime - startT;  % Get the optimization runtime.

%% Use the optimized matrices for smoothing and phase shift estimation.
[minCost, idx] = min(Cost);

Q = reshape(Q_pop(idx, 1:qSize), numStates, numStates);
Q = Q*Q';
R = R_pop(idx);

% Simulate Kalman Smoother using the optimal Q, R and normalized signal
[x_smooth, y_smooth, ~] = KS.simulateKalmanSmoother(A, C, Q, R, y_norm);

% Reconstruct the signal filling in gaps with the smoothed estimate
nan_idx = isnan(y_norm);
y_rec = y_norm;
y_rec(nan_idx) = y_smooth(nan_idx);

% Compute spectra and final cost
[rectifiedSpectrum, f, ~] = Utils.computeSpectrum(t, y_rec);
filteredSpectrum = Utils.computeSpectrum(t, y_smooth);
estimCost = Utils.computeCost(rectifiedSpectrum, filteredSpectrum, f, order);

day1 = 7;           % First day to calculate phase shift.
day2 = 14;          % Second day to calculate phase shift.

[theta, estimPhaseShift] = estimate_phase_shift(day1, day2, omg, t, ...
                                                    x_smooth(1,:), x_smooth(2, :));
phi = compute_phi(t, theta, omg);

fprintf("Optimization Runtime: %s\n", runTime)
fprintf("Cost: %3.3f\n", estimCost)
fprintf("Phase shift: %3.4f hr\n", estimPhaseShift)

%% Plots.

% % Plot original vs smoothed signal
% figure()
% plot(t, y_norm)
% hold on
% plot(t, y_smooth)
% legend('Normalized Original Signal','Smoothed Signal')
% title('Kalman Smoother Estimation')
% xlabel('Time (h)')
% grid on; axis padded;
% hold off
% 
% % Plot states.
% figure();
% for i = 1:numStates
%     subplot(numStates, 1, i)
%     plot(t, x_smooth(i,:))
%     title(strcat("x\_smooth_{", num2str(i), "}"))
%     xlabel("Time (h)")
%     xticks(0:24:t(end)); grid on; axis padded;
% end
% 
% % Plot theta and phi evolution.
% figure()
% subplot(2,1,1)
% plot(t, theta)
% title("\theta")
% xticks(0:24:t(end)); grid on; axis padded;
% 
% subplot(2,1,2)
% plot(t, phi)
% title("\phi")
% xticks(0:24:t(end)); grid on; axis padded;
% xlabel("Time (h)")
% 
% % Plot spectra.
% figure()
% hold on
% plot(f(1:30),rectifiedSpectrum(1:30),'b')
% plot(f(1:30),filteredSpectrum(1:30),'y')
% xregion(f(1), f(3), 'FaceColor', 'g', 'FaceAlpha', 0.2);
% xregion(f(12), f(18), 'FaceColor', 'g', 'FaceAlpha', 0.2);
% legend({'Rectified Spectrum','Filtered Rectified Spectrum $\hat{Y}(w)$'}, 'Interpreter', 'latex');
% title(strcat('Frequency Of Subject ', num2str(subject)))
% hold off
% 
% % Plot average costs of each output.
% figure()
% semilogy(avgCost)
% xlabel("Iteration (k)")
% title("Average Cost (Genetic Algorithm)")
% grid on;

%% Save the values.
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
else
    option_info = "Unknown_";
end


% Also saving mu_y and sigma_y so the signal can be denormalized later if needed.
save(strcat("A", num2str(subject), "_KS_",option_info, num2str(run),".mat"),...
        "x_smooth", "y_smooth", "theta", "phi", "Q", "R", "estimCost", "mu_y", "sigma_y","minCost","y_norm")

end