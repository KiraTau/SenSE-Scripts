% Rensselaer Polytechnic Institute - Julius Lab
% SenSE Project
% Author - Chukwuemeka Osaretin Ike

% Provides general utilities for the optimization and evaluation of 
% both filters.
classdef Utils
    methods (Static, Access='public')

        function [t, y, I, steps] = loadRealData(subject)
            % Loads the data for the specified subject

%             filename = strcat('..\..\Data\A', num2str(subject), '.mat');
            filename = strcat('Data\A', num2str(subject), '.mat');
            load(filename, 'Data')
            t = Data(:,1);
            steps = Data(:,5);
            I = Data(:,6);
            y = Data(:,11);
        end

        function [t, y, phaseShift] = loadSynthData(subject)
            % Loads the data for the specified subject
            
            filename = strcat('Data\SynthData', num2str(subject), '.mat');
            load(filename, 't', 'y', 'phaseShift')
        end
        
        function L = loadOptimalGains(filename)
            load(filename, 'optimGain')
            L = optimGain;
        end
        
        function [P, f, len] = computeSpectrum(t, y)
            T = t(2) - t(1);             % Sampling period
            Fs = 1/T;                    % Sampling frequency
            len = length(t);             % Length of signal
            Y = fft(y);                  % True output 
            P2 = abs(Y/len);             % Compute the two-sided spectrum P2
            P = P2(1:floor(len/2)+1);    % Compute the single-sided spectrum P
            P(2:end-1) = 2*P(2:end-1);   % Double everything except the DC term
            f = Fs*(0:floor(len/2))/len; % Freq based on the sampling rate
        end

        % Computes the FFT of the ANF output, then computes the cost 
        function cost = computeCost(originalSpectrum, filteredSpectrum, f, order)
            % Get the indices of f closest to each harmonic.
            [~, N1] = min(abs(f - 1/24)); 
            [~, N2] = min(abs(f - 2/24));
            [~, N3] = min(abs(f - 3/24));
            [~, N4] = min(abs(f - 4/24));
            [~, N5] = min(abs(f - 5/24));
            [~, n1] = min(abs(f - 0.0309)); 
            NN = N1 - n1;
            harmonicIdxs = [N1, N2, N3, N4, N5];

            % J_harm is the square error within the band around each specified 
            % harmonic and the DC term.
            J_harmo = trapz((filteredSpectrum(1:NN) - originalSpectrum(1:NN)).^2); % DC term.
            
            % J_noise is the square of the signal outside the bands around each
            % harmonic and beyond the last one.
            J_noise = trapz(filteredSpectrum(NN+1:N1-NN-1).^2); % DC term to 1st Harmonic.
            
            for i = 1:order
                idx = harmonicIdxs(i);
                J_harmo = J_harmo + trapz(...
                    (filteredSpectrum(idx-NN:idx+NN) -...
                        originalSpectrum(idx-NN:idx+NN)).^2 );
                if i ~= order
                    idx2 = harmonicIdxs(i+1);
                    J_noise = J_noise + trapz(filteredSpectrum(idx+NN+1:idx2-NN-1).^2);
                end
            end
            J_noise = J_noise + trapz(filteredSpectrum(idx+NN+1:end).^2);

            cost = J_harmo + J_noise;
        end

        function theta = calculateTheta(xHat, omg)
            % Get the first harmonic and its derivative.
            x1 = xHat(1,:);     x2 = xHat(2,:);
            theta = mod(-atan2(x2, omg*x1) + pi/2, 2*pi)-pi;
        end

        function theta_d_hours = calculateContinuousPhase(t, theta, omg)
            theta = reshape(theta, size(t));
            theta_unwrapped = unwrap(theta);
            theta_detrended = theta_unwrapped - (omg*t);
            theta_d_hours = theta_detrended/omg;
        end

        function [theta, phaseShift] = estimatePhaseShift(xHat, omg, day1, day2)
            % Estimate the phase shift in hours given the OBF output and state
            % evolution
            
            % Isolate the first Harmonic and its derivative
            x1 = xHat(1,:);     x2 = xHat(2,:);
            
            % Calculate the phase based on x1 and x2
            theta = mod(-atan2(x2, omg*x1) + pi/2, 2*pi)-pi;
            
            % Calculate the phase shift between day 1 and 2
            day1Series = (( (day1-1)*24*60)+1):(day1*24*60);
            day2Series = (( (day2-1)*24*60)+1):(day2*24*60);
            phaseShift = (1/omg)*mean(unwrap(theta(day1Series)) - unwrap(theta(day2Series)) );
            phaseShift = mod(12+phaseShift, 24) - 12;
        end
        
        function [markers, phaseShift] = estimateCBTPhaseShift(t, x, day1, day2)
            T = 24; % period of the signal.
            
            % Gather the minimum values in every 24 hour period.
            intervals = t(1):T:t(end);
            markers = zeros(length(intervals)-1,1);
            for i = 2:length(intervals)
                [~, idx] = min(x(t >= intervals(i-1) & t < intervals(i)));
                markers(i-1) = t(((i-2)*1440)+idx);
            end
            markers = mod(markers,24);
            phaseShift = markers(day2)-markers(day1);
        end

        function [min_values, min_indices] = find_daily_min(x, window_size)
            % Calculate the number of groups.
            num_groups = floor(length(x)/window_size);

            % Initialize a vector to store the values.
            min_values = zeros(num_groups, 1);
            min_indices = zeros(num_groups, 1);

            % Loop through the data in chunks of window_size.
            for i = 1:num_groups
                start_idx = (i-1)*window_size+1;
                end_idx = i*window_size;
                [min_values(i), min_indices(i)] = min(x(start_idx:end_idx));
                min_indices(i) = min_indices(i)+start_idx-1;
            end
        end

        function plotAvgCostEvolution(avgCost)
            % Plot the average cost evolution.
            figure(1); clf
            semilogy(avgCost)
            title('Evolution of Average Population Cost')
            xlabel('Iteration')
            ylabel('Cost')
            grid on; axis padded
        end
        
        function plotFilterOutputs(t, y, yHat, xHat, theta, day1, day2, label)
            % Plot the outputs.
            figure(2); clf
            subplot(2,1,1)
            plot(t, y, 'b-.')
            hold on
            plot(t, yHat, 'r-.', 'LineWidth', 2)
            hold off
            legend('y', 'yHat')
            title(strcat("Filter Output vs ", label))
            xlabel('t (h)')
            ylabel(strcat(label, " Magnitude"))
            xticks(t(1):72:t(end)); grid on; axis padded;
            
            % Plot continuous phase offset.
            subplot(2,1,2)
            theta = reshape(theta, size(t));
            thetaU = unwrap(theta);
            thetaD = thetaU - ((2*pi/24)*t);
            thetaDHours = thetaD*24/2/pi;
            plot(t,...
                thetaDHours - mean( thetaDHours(t >= (day1-1)*24 & t <= day1*24) )...
            ); 
            hold on; 
            xline(t(t==day1*24), 'LineWidth', 2); 
            xline(t(t==(day1-1)*24), 'LineWidth', 2);
            xline(t(t==day2*24), 'LineWidth', 2); 
            xline(t(t==(day2-1)*24), 'LineWidth', 2); 
            hold off
            title(strcat("Phase Offset Relative to Mean Phase on Day ", num2str(day1)))
            xlabel('Time (h)'); ylabel('Phase Offset (h)'); 
            xticks(t(1):72:t(end)); grid on; axis padded;

            % Plot the state estimates.
            figure(3); clf
            subplot(3,1,1); plot(t, xHat(1,:), 'k-'); legend('Harmo_1')
            title('State Estimates')
            xticks(t(1):72:t(end)); grid on; axis padded;
            subplot(3,1,2); plot(t, xHat(2,:), 'k-'); legend('Harmo_1 Derivative')
            xticks(t(1):72:t(end)); grid on; axis padded;
            subplot(3,1,3); plot(t, xHat(end,:), 'k-'); legend('Bias Term')
            xlabel('t (h)'); 
            xticks(t(1):72:t(end)); grid on; axis padded;
        end
        
        function plotSpectra(f, originalSpectrum, filteredSpectrum, label)
            plot(f, originalSpectrum);
            hold on;
            plot(f, filteredSpectrum);
            hold off;
            legend('Original Spectrum', 'Filtered Spectrum')
            title(strcat("Original and Filtered Spectra of ", label))
            xlim([0 1])
        end

        % Functions for analysis.
        % .
        function plotPhaseShiftBoxPlot(results, timeShift, subject, dlmoSubject, filter)
            boxplot(results(results(:,1) == subject, 4))
            if ~isnan(timeShift(dlmoSubject,1))
                yline(timeShift(dlmoSubject,1), '-.r', 'LineWidth', 1.5, 'DisplayName', 'Lewy');
            end
            if ~isnan(timeShift(dlmoSubject,2))
                yline(timeShift(dlmoSubject,2), '-.k', 'LineWidth', 1.5, 'DisplayName', 'Eastman');
            end
            if ~isnan(timeShift(dlmoSubject,3))
                yline(timeShift(dlmoSubject,3), '-.c', 'LineWidth', 1.5, 'DisplayName', 'Voultsios');
            end
            hold off; legend; ylim([-5 5]);
            title(strcat(filter, " Phase Shift Estimate for Subject ", num2str(subject)));
            ylabel('Phase Shift (h)');
        end

        function plotPhaseShiftHistogram(results, timeShift, subject, dlmoSubject, filter)
            histogram(results(results(:,1) == subject, 4), 24, 'DisplayName', 'OBF Estimates')
            if ~isnan(timeShift(dlmoSubject,1))
                xline(timeShift(dlmoSubject,1), '-.r', 'LineWidth', 1.5, 'DisplayName', 'Lewy'); 
            end
            if ~isnan(timeShift(dlmoSubject,2))
                xline(timeShift(dlmoSubject,2), '-.k', 'LineWidth', 1.5, 'DisplayName', 'Eastman'); 
            end
            if ~isnan(timeShift(dlmoSubject,3))
                xline(timeShift(dlmoSubject,3), '-.c', 'LineWidth', 1.5, 'DisplayName', 'Voultsios');
            end
            hold off; legend;
            title(strcat(filter, " Phase Shift Estimates for Subject ", num2str(subject))); 
            xlabel('Phase Shift (h)'); ylabel('Count')
            xlim([-3 3])
        end

        function plotCosts(results, subject, filter)
            histogram(results(results(:,1) == subject, 6), 24)
            title(strcat(filter, " Optimal Cost for Subject ", num2str(subject)));
            ylabel('Count'); xlabel('Optimal Costs')
            xlim([9e4 2e6])
        end

        function plotRuntimes(results, subject, filter)
            boxplot(results(results(:,1) == subject, 5))
            title(strcat(filter, " Runtimes for Subject ", num2str(subject)));
            ylabel('Runtime (s)');
            ylim([0 100])
        end
        
        
    end
end