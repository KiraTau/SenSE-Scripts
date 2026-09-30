% TJU theta experiments: wrist-activity light and saved KS run 11.
subjects = [7 8 1 5];
options = 2:11;
%2Actigraphy	3BPM	4Ave IBI (ms)	5Min IBI (ms)	6Max IBI (ms)	7SD
%8RMSSD	9LF (ms2, 0.04-0.15 Hz)	10HF (ms2, 0.15-0.40 Hz)	11LF/HF	12
%Wrist actigraphy
run = 11;
tests = 20:29;
Q = [5e-3 5e-3 1e-3 1e-2 1e-2];
R = 0.85;
rerun_kalman = false; % Reuse saved KS run 11; set true to rerun the smoother.

for subject = subjects
    H = load(fullfile('MatData', "NSF0" + subject + "_HRV.mat"));
    W = load(fullfile('MatData', "NSF0" + subject + "_MTN.mat"));
    [heart, wrist, common_start_date, days_duration] = match_time( ...
        H.HRV_person, W.MotionData, H.start_day, W.start_day);

    switch subject
        case {7, 8}
            t = wrist(:,1);
            activity = wrist(:,2);

        case 1
            % Preserve the existing padded KS layout, including its endpoint.
            fill_time = (heart(end,1):1/60:14.25*24)';
            heart = [heart; fill_time zeros(numel(fill_time),11)];
            wrist = [wrist; fill_time zeros(numel(fill_time),2)];
            t = heart(:,1);
            activity = wrist(:,2);

        case 5
            % Raw wrist time is in days; matched heart time is in hours.
            minute_counts = W.MotionData(:,1)*24*60;
            wrist_clean = W.MotionData(abs(minute_counts-round(minute_counts))<1e-5,:);
            first_idx = find(abs(wrist_clean(:,1)-wrist(1,1)/24)<1e-6,1);
            fill_time = (heart(end,1)+1/60:1/60:wrist_clean(end,1)*24)';
            heart = [heart; fill_time zeros(numel(fill_time),11)];
            t = heart(:,1);
            activity = wrist_clean(first_idx:end,2);
    end

    wrist_light = activity_to_light(activity);
    for option = options
        y = heart(:,option);
        I = wrist_light;
        if subject ~= 1 && ~any(I)
            I = activity_to_light(y);
        end

        generate_initial_conditions_huang(subject)
        if rerun_kalman
            run_kalman_smoothing(t,y,run,subject,option,run)
        end

        for test = tests
            run_particle_filter_delta_tau_ks_theta_compressed( ...
                subject,run,test,Q,R,I,t,option)
        end
    end
end
