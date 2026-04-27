%This is the main function for the whole package
%create the folder
folderName = fullfile(pwd, 'Results');

% Check if the folder exists at this exact absolute path
if ~exist(folderName, 'dir')
    mkdir(folderName);
    fprintf('Folder "%s" created.\n', folderName);
else
    fprintf('Folder "%s" already exists.\n', folderName);
end

%Load the DLMO data, DLMO considered as 7 hours before minimum core body
%temperature
load("DLMO.mat")

%load subject 8
test_subject = 8;
run = 1;
load(strcat('Data/A',num2str(test_subject),'_MTN'))
t = MotionData(:,1)*24;%Time in days, convert to hours
y = MotionData(:,2);%actigraph data
I = MotionData(:,3);%light data

if ~any(I)%if there are not light, generate light with activity
    I = activity_to_light(y);
end

generate_initial_conditions_huang(test_subject)%generate the initial condition and add that to the file
optimize_kf(t,y,run,test_subject)%This process compute the parameters of the kalman filter and system model
% Estimate the circadian phase
% The process is stochastic so you can run it multiple times to obtain the best estimate

load(strcat("Data\A", num2str(test_subject), "_KF_", num2str(run),".mat"));
%Particle filter is stochastic, so you can run multiple times to find the
%best estimate
run_particle_filter_delta_tau(test_subject, run, Q, R, I, t)

%size of this file is very big, loading can be slow
load(strcat("Results\A", num2str(test_subject), "_PF_Out_Delta_Tau_", num2str(run), ".mat"))
figure()
subplot(2,1,1)
plot(t,x_mean(:,1))
title('Estimated x')
xticks(0:24:t(end)); 
grid on; axis padded;

subplot(2,1,2)
plot(t,x_mean(:,2))
title("Estimated  x_c")
xticks(0:24:t(end)); 
grid on; axis padded;
xlabel("Time (h)")

load(strcat('Data\A',num2str(test_subject),'_MTN'))
t = MotionData(:,1);%in days
min_x = [];
%Assume that minimum core body temperature occurs during night time
for j = 1:round(t(end))
    day_t = t(find(round(t)==j));
    day_x = x_mean(find(round(t)==j),1);
    [xmin,idx] = min(day_x);
    min_x = [min_x day_t(idx)*24];
end
%assumes dlmo = core body temperature min - 7h
dlmo_hour = mod(min_x-7,24);
save(strcat('Results/A',num2str(test_subject),'min_CBT_time'),"dlmo_hour")

disp(strcat('Day 7 Estimated DLMO ',num2str(dlmo_hour(7))))
disp(strcat('Day 7 Salivary DLMO ',num2str(DLMO(1))))

disp(strcat('Day 14 Estimated DLMO ',num2str(dlmo_hour(14))))
disp(strcat('Day 14 Salivary DLMO ',num2str(DLMO(2))))