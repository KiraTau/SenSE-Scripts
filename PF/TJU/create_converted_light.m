function create_converted_light(subject)
% Rensselaer Polytechnic Institute - Julius Lab
% SenSE Project
% Author - Chukwuemeka Osaretin Ike
%
% Description:
% Converts actigraphy to light and saves the new data in a mat file.

% subject = 10;

[t, actigraphy, I, steps] = Utils.loadRealData(subject);

step_light = activity_to_light(steps);
act_light = activity_to_light(actigraphy);
conv_light = step_light;

%%
figure(1);
subplot(3,1,1)
plot(t, I)
title("Measured Light")
xticks(0:24:t(end))
xlabel("Time (h)")

subplot(3,1,2)
plot(t, step_light)
title("Light from Steps")
xticks(0:24:t(end))
xlabel("Time (h)")

subplot(3,1,3)
plot(t, act_light)
title("Light from Actigraphy")
xticks(0:24:t(end))
xlabel("Time (h)")

%%
figure(2)
subplot(2,1,1)
plot(t, actigraphy)
title("Actigraphy")
xticks(0:24:t(end))
xlabel("Time (h)")

subplot(2,1,2)
plot(t, steps)
title("Steps")
xticks(0:24:t(end))
xlabel("Time (h)")

%% Find the hour that the subject's data started at.
warning off
% source_dir = "G:\My Drive\Career\Projects\SenSE\Code\Data\UNM Originals\";
source_dir = "Data\UNM Originals\";
filename = dir(strcat( ...
    source_dir, "A", ...
    num2str(subject), '_*.csv')).name;
subject_data = readtable(strcat(source_dir, filename));

start_time = hours(subject_data.Time(1));
start_hour = floor(start_time);

save(strcat("Data\A", num2str(subject), ".mat"), "conv_light",...
          "start_hour", "start_time","t", "actigraphy", "I", "steps", "-append")
% save(strcat("Data\A", num2str(subject), ".mat"), "conv_light",...
%           "start_hour", "start_time", "-append")
% save(strcat("A", num2str(subject), ".mat"), ...
%     "t", "actigraphy", "I", "steps", "conv_light")
end