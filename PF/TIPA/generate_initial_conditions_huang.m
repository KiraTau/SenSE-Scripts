function generate_initial_conditions_huang(subject)
% Rensselaer Polytechnic Institute - Julius Lab
% SenSE Project
% Author - Chukwuemeka Osaretin Ike
%
% Description:
% Adds initial conditions to the mat file. The initial conditions are
% gotten based on the method from Huang (2021).
% We generate a trajectory with 16 hours of light and 8 of darkness around
% average habitual wake of 7am and sleep of 11pm.
% After generating that trajectory and ensuring it converges, we look at
% the limit cycle on the last day and take the average x, xc, n for the
% hour that subject's data starts.
% These conditions are used for both the converted and measured light
% cases.

%% 
T = 24;                 % Light input period - 24 hours.
dt = 1/60;              % Sampling time - 1 minute.
omg = 2*pi/24;          % Fundamental frequency in rad/h.

numDays = 30;           % Length of data in days.
numSamplesPerDay = 1440;
numHours = T*numDays;   % Length of data in hours.
offset = 0.5;           % Offset to shift the square waves.

% Create the time and light vectors.
t = (0:dt:numHours)';  % Time vector.

I_max = 800;           % Max light value.

% Duty cycle is flipped because we start at midnight with 0 light.
I = I_max*(1 - (offset*square(t*omg, 100/3) + offset))';

init_cond = [0.1, 0.1, 0.5];        % x, xc, n.

% Generate the trajectories.
[t_out, ~, y_out] = jewett99(t, I, init_cond);
cbt_mins = getCBTMins(t_out, y_out(:,1), y_out(:,2), "jewett99");

%% Plots.
figure(1); clf;

subplot(2,1,1)
plot(t, I)
xticks(0:24:t(end))
xlabel("Time (h)")

subplot(2,1,2)
plot(t_out, y_out(:,1), t_out, y_out(:,2))
xline(cbt_mins, "--")
legend("x", "xc", "CBT_{min}")
xticks(0:24:t(end))
xlabel("Time (h)")
title("State Trajectories")


%% Set initial conditions.
% Find the average x, xc, n for the hour that the subject's data
% started at.
data_file = strcat("Data\A", num2str(subject), "_MTN.mat");
load(data_file, "MotionData")
start_hour = floor(MotionData(1,1)*24);

last_day_t = t_out(end-numSamplesPerDay+1:end);
last_day_y = y_out(end-numSamplesPerDay+1:end, :);
last_day_hour = mod(last_day_t, 24);
init = mean(last_day_y(floor(last_day_hour) == start_hour, :), 1);

%% Save the initial conditions into the mat file.
save(data_file, "init", "-append")
end