% Rensselaer Polytechnic Institute - Julius Lab
% SenSE Project
% Author - Chukwuemeka Osaretin Ike
% 
% Steady-State Kalman Filter Class.
%
% Description:
% Contains functions for the optimization and simulation of the
% Steady-State Kalman Filter.

classdef KFSS
    methods (Static, Access='public')
        function [xHat, yHat] = simulateFilter(A, C, L, y)
            seqLength = length(y);
            xHat = zeros(size(A,1), seqLength);
            xHat(end, 1) = 1000;

            % Run the simulation. If y=0, run without correction.
            for i = 2:seqLength
                k = i-1;
                if y(k) == 0
                    xHat(:,i) = A*xHat(:,k);
                else
                    xHat(:,i) = (A - L*C)*xHat(:,k) + L*y(k);
                end
            end
            yHat = C*xHat;
        end

        function [Cost, avgCost, Q_pop, R_pop] = optimizeFilter(...
                iterations, mu, lambda, rho, order, lEnd, rEnd, LB, rLA, rLB,...
                A, C, t, y, f, originalSpectrum)

            % Arrays to hold the costs.
            Cost = zeros(mu,1);    
            avgCost = zeros(iterations, 1);

            stateLength = size(A, 1);
            qSize = stateLength^2;
            [Q_pop, R_pop] = KF.initializePopulation(stateLength, mu,...
                lEnd, rEnd, LB, rLA, rLB);

            for member = 1:mu
                Q = reshape(Q_pop(member, 1:qSize), stateLength, stateLength);
                Q = Q*Q';
                R = R_pop(member);

                % Compute P from the ARE.
                [P, ~, ~] = idare(A', C', Q, R, [], []);
                if isempty(P)
                    reshape(Q_pop(member, 1:qSize), stateLength, stateLength)
                    Q
                    R
                    Cost(member) = NaN;
                    continue;
                end
                L = A*P*C'/(C*P*C' + R);
                % simStart = clock;
                [~, yHat] = KFSS.simulateFilter(A, C, L, y);
                % simEnd = clock;
                % fprintf("Simulation Runtime: %3.3f\n", etime(simEnd, simStart))
                filteredSpectrum = Utils.computeSpectrum(t, yHat);
                Cost(member) = Utils.computeCost(originalSpectrum, filteredSpectrum, f, order);
            end

            for iteration = 1:iterations
                % All possible combinations of [1:50] in pairs then take 1st 50.
                Combination = nchoosek([1:50], rho);
                Labels = randperm(size(Combination, 1) );

                % Add 50 members to the gene pool, simulate the dynamics
                % each time and collect the costs.
                for j = 1:lambda
                    m = mu + j;

                    % Make offspring using the mean then get Q and R.
                    Q_pop(m,:) = mean(Q_pop(Combination(Labels(j),:), :));
                    R_pop(m) = mean(R_pop(Combination(Labels(j), :), :));

                    Q = reshape(Q_pop(m, 1:qSize), stateLength, stateLength);
                    Q = Q*Q';
                    R = R_pop(m);

                    % Compute P from the ARE.
                    [P, ~, ~] = idare(A', C', Q, R, [], []);
                    if isempty(P)
                        reshape(Q_pop(m, 1:qSize), stateLength, stateLength)
                        Q
                        R
                        Cost(m) = NaN;
                        continue;
                    end
                    L = A*P*C'/(C*P*C' + R);

                    [~, yHat] = KFSS.simulateFilter(A, C, L, y);
                    filteredSpectrum = Utils.computeSpectrum(t, yHat);
                    Cost(m) = Utils.computeCost(originalSpectrum, filteredSpectrum, f, order);
                end

                % Remove the lambda highest costs - maxk only removes real values,
                % so we have to prioritize removing NaN values ourselves.
                nan_idx = find(isnan(Cost));
                [~,n] = maxk(Cost, lambda);
                remove_idx = [nan_idx; n];
                Cost(remove_idx(1:lambda)) = [];
                Q_pop(remove_idx(1:lambda), :) = [];
                R_pop(remove_idx(1:lambda), :) = [];

                % Collect the average cost of the iteration.
                avgCost(iteration) = mean(Cost, 'omitnan');
            end
        end
    end
end