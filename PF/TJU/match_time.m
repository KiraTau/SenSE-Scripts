function [matched_heart, matched_wrist,common_start_date,days_duration] = match_time(heart, wrist, start_day_1, start_day_2)
    % Convert everything to Absolute Hours (Day 1 at 00:00 = Hour 0)
    % This handles the different start days perfectly.
    base_date = min(datetime(start_day_1), datetime(start_day_2));
    
    % Get absolute hours for heart
    h_abs_start = 24 * (datenum(datetime(start_day_1)) - datenum(base_date));
    h_times = h_abs_start + (heart(:,1) * 24); % Convert heart days to hours
    
    % Get absolute hours for wrist
    % Convert wrist days to hours
    wrist_hrs = wrist(:,1) * 24;
    
    % Identify the "on-the-minute" marks
    % We multiply by 60 so that every integer (1, 2, 3...) represents a minute.
    minute_counts = wrist_hrs * 60;
    
    % Create a mask for rows that are extremely close to an integer minute
    % We use a tolerance of 1e-5 to catch those :59 second rows
    is_minute_mark = abs(minute_counts - round(minute_counts)) < 1e-5;
    wrist_clean = wrist(is_minute_mark, :);

    w_abs_start = 24 * (datenum(datetime(start_day_2)) - datenum(base_date));
    w_times = w_abs_start + (wrist_clean(:,1) * 24); % Convert wrist days to hours


    % Find the Shared Time Window (Intersection)
    common_start = max(min(h_times), min(w_times));
    common_end   = min(max(h_times), max(w_times));


    idx_h_start = find(abs(h_times - common_start) < 1e-7);
    idx_h_end   = find(abs(h_times - common_end) < 1e-7);
    
    % Find the start and end indices for Wrist
    idx_w_start = find(abs(w_times - common_start) < 1e-7);
    idx_w_end   = find(abs(w_times - common_end) < 1e-7);
    
    % Crop the original arrays using these indices
    matched_heart = [h_times(idx_h_start:idx_h_end)-h_abs_start,heart(idx_h_start:idx_h_end, 2:end)];
    matched_wrist = [w_times(idx_w_start:idx_w_end)-w_abs_start,wrist_clean(idx_w_start:idx_w_end, 2:end)];
    
    common_start_date = max(datetime(start_day_1), datetime(start_day_2));

    days_duration = (matched_heart(end,1) - matched_heart(1,1))/24;
end