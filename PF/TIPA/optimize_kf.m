function optimize_kf(t_,y_,run,subject)
% Rensselaer Polytechnic Institute - Julius Lab
% SenSE Project
% Original Author - Chukwuemeka Osaretin Ike
% Modified by Zidi Tao
%
% Description:
% 

%% Parameters.

% Optimization Hyperparameters.
iterations = 150;            % Maximum iterations.
mu = 100;                   % Number of every generation.
lambda = 50;                % Number of children.
rho = 2;                    % Number of parents to create one child.

% Optimization bounds.
LB = -5;                % Q Lower bound - 10^LB.
lEnd = 0;               % Q lower end.
rEnd = 18;              % Twice the number of orders of magnitudes.
rLA = 1e2;              % R Lower bound.
rLB = 1e8;              % R upper bound. 

% Filter params.
order = 1;
fprintf("Filter Order %d\n", order)
numStates = (2*order + 1);        % State length based on order
qSize = numStates^2;
omg = 2*pi/24;

y = y_';
t = t_';


% Take the FFT of original signal for cost calculation
numOutputs = size(y,1);
titles = ["Actigraphy", "Light"];

[originalSpectrum, f, ~] = Utils.computeSpectrum(t, y);

%% Optimize the filter.
dt = t(2) - t(1); % Get the sampling time for the subject.
[A, ~, C, ~] = KF.createKalmanStateSpace(order, numStates, dt);

startT = datetime;                  % Start timing the optimization.
[Cost, avgCost, Q_pop, R_pop] = KFSS.optimizeFilter(...
            iterations, mu, lambda, rho, order, lEnd, rEnd, LB, rLA, rLB,...
            A, C, t, y, f, originalSpectrum);
runTime = datetime - startT;  % Get the optimization runtime.

%% Use the optimized matrices in estimating the phase shift.
[~, idx] = min(Cost);

Q = reshape(Q_pop(idx, 1:qSize), numStates, numStates);
Q = Q*Q';
R = R_pop(idx);

% Compute P from the ARE.
[P, ~, ~] = idare(A', C', Q, R, [], []);
L = A*P*C'/(C*P*C' + R);


% Reload the data.
y = y_';
t = t_';

day1 = 7;           % First day to calculate phase shift.
% Second day to calculate phase shift.
day2 = 14;


[xHat, yHat] = KFSS.simulateFilter(A, C, L, y);

[originalSpectrum, f, ~] = Utils.computeSpectrum(t, y);
filteredSpectrum = Utils.computeSpectrum(t, yHat);
estimCost = Utils.computeCost(originalSpectrum,...
    filteredSpectrum, f, order);

[theta, estimPhaseShift] = estimate_phase_shift(day1, day2, omg, t, ...
                                                    xHat(1,:), xHat(2, :));
phi = compute_phi(t, theta, omg);

% [dailyPhases, thetaDReasonable] = get_average_daily_phase(t, theta, omg);

fprintf("Optimization Runtime: %s\n", runTime)
fprintf("Cost: %3.3f\n", estimCost)
fprintf("Phase shift: %3.4f hr\n", estimPhaseShift)

%% Plots.

% Plot states.
figure(1);
for i = 1:numStates
    subplot(numStates, 1, i)
    plot(t, xHat(i,:))
    title(strcat("x_", num2str(i)))
    xlabel("Time (h)")
    xticks(0:24:t(end)); grid on; axis padded;
end


% Plot theta and thetaDot evolution.
figure(2)
subplot(2,1,1)
plot(t-t(1), theta)
title("\theta")
xticks(0:24:t(end)); grid on; axis padded;

subplot(2,1,2)
plot(t, phi)
title("\phi")
xticks(0:24:t(end)); grid on; axis padded;
xlabel("Time (h)")



% Plot spectra.
figure(3)
for d = 1:numOutputs
    subplot(numOutputs, 1, d)
    Utils.plotSpectra(f(d,:), originalSpectrum(d,:), filteredSpectrum(d,:), titles(d))
    hold on;
    xline(0, 'LineWidth', 2); 
    xline(1/24, 'LineWidth', 2); 
    xline(2/24, 'LineWidth', 2); 
    xline(3/24, 'LineWidth', 2); 
    hold off
    xlabel("Frequency (Hz)")
    ylabel("Magnitude")
    xlim([0 5/24])
end

% Plot average costs of each output.
figure(4)
for d = 1:numOutputs
    subplot(numOutputs, 1, d)
    semilogy(avgCost(:,d))
    xlabel("Iteration (k)")
    title(strcat("Average Cost for ", titles(d)))
end
%% Save the values.
save(strcat("Data\A", num2str(subject), "_KF_", num2str(run),".mat"),...
        "xHat", "yHat", "theta", "phi", "Q", "R")

end