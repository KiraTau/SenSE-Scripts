% TJU experiments: subject-specific preparation, shared signal/trial loops.
subjects = [7 8 1 5];
options = 2:11;
run = 11;
tests = 20:29;
Q = [5e-3 5e-3 1e-3 1e-2 1e-2];
R = 0.85;
rerun_kalman = false; % Reuse saved run 11. Set true to regenerate KS and init.

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
            pf_function = @run_particle_filter_delta_tau_ks_phi_compressed;
    end

    wrist_light = activity_to_light(activity);
    for option = options
        y = heart(:,option);
        I = wrist_light;
        if subject ~= 1 && ~any(I)
            I = activity_to_light(y);
        end

        if rerun_kalman
            generate_initial_conditions_huang(subject)
            run_kalman_smoothing(t,y,run,subject,option,run)
        end

        for test = tests
            run_particle_filter_delta_tau_ks_phi_compressed(...
                subject,run,test,Q,R,I,t,option)
        end
    end
end
