% Rensselaer Polytechnic Institute - Julius Lab
% SenSE Project
% Author - Chukwuemeka Osaretin Ike
% 
% Kalman Filter Class.
%
% Description:
% Contains functions for the optimization and simulation of the
% Kalman Filter.

classdef KF
    methods (Static, Access='public')
        function [A, B, C, D] = createKalmanStateSpace(order, stateLength, dt)
            w = 2*pi/24;
            A = zeros(stateLength); 
            C = zeros(1, stateLength);
            for k = 1:order
                i = k*2;
                A(i-1:i, i-1:i) = [0, 1; -(k*w)^2, 0];
                C(i) = 2/(k*w);
            end
            C(end) = 1;
            B = zeros(stateLength, 1);
            D = 0;
            contSystem = ss(A, B, C, D);

            % Discretize the system with dt sample time.
            discSystem = c2d(contSystem, dt, 'impulse');
            A = discSystem.A;
            B = discSystem.B;
            C = discSystem.C;
            D = discSystem.D;
        end

        function [Q_pop, R_pop] = initializePopulation(stateLength, mu,...
                                            lEnd, rEnd, LB, rLA, rLB)
            % Create the initial Q population.
            qSize = stateLength^2;
            mid = (rEnd + lEnd)/2;
            N = (rEnd-lEnd)*rand(mu, qSize) + lEnd;
            Q_pop = zeros(size(N));
            Q_pop(N > mid+.5) = 10.^(N(N > mid+.5) - mid + LB);
            Q_pop(N < mid-.5) = -10.^(mid - N(N < mid-.5) + LB);

            % Create the initial R population.
            R_pop = (rLA + (rLB-rLA)*rand(mu, 1));
            
        end
        
        function [xHat, yHat] = simulateKalmanFilter(A, C, Q, R, y)
            % Runs the Kalman Filter given the appropriate settings.
            seqLength = length(y);
            stateLength = size(A, 1);
            xHat = zeros(stateLength, seqLength);
            xHat(end, 1) = 1500;
            yTilde = zeros(1, seqLength);
            x_k = zeros(stateLength, 1);
            x_k(end, 1) = 1500;
            I = eye(size(A));
            P = 0.01*eye(size(Q));
%             Ks = xHat;

            for i = 2:seqLength
        %         k = i-1;

                % Prediction step.
                xHat(:,i) = A*x_k;
                PHat = A*P*A' + Q;

                % Update step.
                yTilde(i) = y(i) - C*xHat(:,i);
                S = C*PHat*C' + R;
                K = PHat*C'/S;
%                 Ks(:,i) = K;
                x_k = xHat(:,i) + K*yTilde(i);
                P = (I - K*C)*PHat;

            end
            yHat = C*xHat;
        end
        
        function [Cost, avgCost, Q_pop, R_pop] = optimizeFilter(...
            iterations, mu, lambda, rho, order, lEnd, rEnd, LB, rLA, rLB,...
            A, C, t, y, f, originalSpectrum)

            % Arrays to hold the costs.
            Cost = zeros(mu, 1);    
            avgCost = zeros(iterations, 1);

            stateLength = size(A, 1);
            [Q_pop, R_pop] = KF.initializePopulation(stateLength, mu, lEnd,...
                                                      rEnd, LB, rLA, rLB);

            for member = 1:mu
                Q = reshape(Q_pop(member, :), stateLength, stateLength);
                Q = Q*Q';
                R = R_pop(member);

                [~, yHat] = KF.simulateKalmanFilter(A, C, Q, R, y);
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

                    Q = reshape(Q_pop(m, :), stateLength, stateLength);
                    Q = Q*Q';
                    R = R_pop(m);

                    [~, yHat] = KF.simulateKalmanFilter(A, C, Q, R, y);
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
                R_pop(remove_idx(1:lambda)) = [];

                % Collect the average cost of the iteration.
                avgCost(iteration) = mean(Cost);
            end
        end
    end
    
end