function cbt_mins = getCBTMins(t, x, xc, model)
%% Rensselaer Polytechnic Institute - Julius Lab
% SenSE Project
% Author - Chukwuemeka Osaretin Ike
%
% Description:
% Calculates CBTmin based on the model and phi_c.

T = 24;

switch model
    case "jewett99"
    %% Jewett (1999).
    % If time vector doesn't start at 0, remove everything before 24h.
    % Assumes that if someone put a watch on, it was after they woke up and
    % thus, CBTmin has passed.
    if t(1) ~= 0
        x(t < T) = [];
        t(t < T) = [];
    end
    new_start = t(1);
    t = t - t(1);

    % [~, idx] = findpeaks(-x);
    % cbt_mins = t(idx);

    % Gather the minimum values in every 24 hour period.
    num_days = ceil(t(end)/T);
    cbt_mins = zeros(num_days, 1);
    for i = 1:num_days
        % [~, idx] = min(x( (t >= (i-1)*T) & (t < i*T) ));
        % cbt_mins(i) = t( idx + length(t(t < (i-1)*T)) );
        t_idx = (t >= (i-1)*T) & (t < i*T);
        x_day = x(t_idx);
        [~, idx] = findpeaks(-x_day, "MinPeakDistance", min(120, length(x_day)-2));
        if isempty(idx)
            cbt_mins(i) = NaN;
        elseif length(idx) == 1
            cbt_mins(i) = t( idx + length(t(t < (i-1)*T)) );
        else
            [~, idx2] = min(x_day(idx));
            cbt_mins(i) = t( idx(idx2) + length(t(t < (i-1)*T)) );
        end
    end

    % Add reference from Jewett (1999).
    phi_c = 0.8;

    % Add new_start to cast cbt_mins back in original time vector.
    cbt_mins = cbt_mins + new_start + phi_c ;
    if new_start ~= 0
        cbt_mins = [NaN; cbt_mins];
    end

    % % Post-process to add NaNs for days with no detected CBTmin.
    % cbt_mins = postprocess(t(end)+new_start, T, cbt_mins);

    case "hilaire07"
    %% St. Hilaire (2007).
    % Time of CBTmin = Time at which [atan(xc/x) = -170.7] + phi_c.

    angles = atan2(xc,x);

    % Find values within a 0.5 deg tolerance of -170.7 degrees.
    target = -170.7/180*pi;
    tol = 1.5/180*pi;
    mins = t(angles <= target+tol & angles > target-tol);

    % The mins above tends to contain > 1 value for each CBTmin.
    % Average values that are really close together (within 1 hour).
    cbt_mins = [];
    cluster_start = 1;
    threshold = 1;
    for i = 2:length(mins)
        if abs(mins(i) - mins(i-1)) > threshold
            cbt_mins = [cbt_mins; mean(mins(cluster_start:i-1))];
            cluster_start = i;
        end
    end
    if cluster_start <= length(mins)
        cbt_mins = [cbt_mins; mean(mins(cluster_start:end))];
    end

    % Add reference from St. Hilaire (2007).
    phi_c = 0.97;
    cbt_mins = cbt_mins + phi_c;

    % Post-process to add NaNs for days with no detected CBTmin.
    cbt_mins = postprocess(t(end), T, cbt_mins);
end

    function c_mins = postprocess(t_end, T, cbt_mins)
        % Add NaNs on days without detected CBTmins.
        numDays = ceil(t_end/T);
        c_mins = zeros(numDays, 1);
        j = 1;
        for idx1 = 1:numDays
            % if cbtmin is not within the current day window, NaN.
            if ~((cbt_mins(j) >= (idx1-1)*T) && (cbt_mins(j) < idx1*T))
                c_mins(idx1) = NaN;
            else
                c_mins(idx1) = cbt_mins(j);
                j = j+1;
            end
        end
    end

end