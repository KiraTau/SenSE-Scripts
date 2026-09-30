% Rensselaer Polytechnic Institute - Julius Lab
% SenSE Project
% Author - Zidi Tao
% 
% Kalman Smoothing Class.
%
% Description:
% Contains functions for the optimization and simulation of the
% Kalman Smoothing.

classdef KS
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

        function [x_smooth,y_smooth,x_post] = simulateKalmanSmoother(A,C,Q,R,y)
            %Forward Kalman filter
            N = length(y);
            x_priori = zeros(3, N);  % x_{k|k-1}
            P_priori = zeros(3, 3, N); % P_{k|k-1}
            x_post   = zeros(3, N);  % x_{k|k}
            P_post   = zeros(3, 3, N); % P_{k|ki
            
            
            
            
            for k = 1:N
                % --- Predict Step ---
                if k == 1
                    x_priori(:,k) = [0;0;0];
                    P_priori(:,:,k) = 1e3*eye(3);
                else
                    x_priori(:,k) = A * x_post(:,k-1);
                    P_priori(:,:,k) = A * P_post(:,:,k-1) * A' + Q;
                end
                
                % --- Update posterior estimates
                if isnan(y(k))
                    %No measurement
                    x_post(:,k) = x_priori(:,k);
                    P_post(:,:,k) = P_priori(:,:,k);
                else
                    % Measurement Available
                    y_tilde = y(k) - C * x_priori(:,k);
                    K = P_priori(:,:,k) * C' / (C * P_priori(:,:,k) * C' + R); % Kalman Gain
                    
                    x_post(:,k) = x_priori(:,k) + K * y_tilde;
                    P_post(:,:,k) = (eye(3) - K * C) * P_priori(:,:,k);
                end
            end
            
            % Backward Pass: RTS Smoother
            x_smooth = zeros(3, N);
            P_smooth = zeros(3, 3, N);
            
            x_smooth(:,N) = x_post(:,N);
            P_smooth(:,:,N) = P_post(:,:,N);
            
            
            for k = N-1:-1:1
                % Smoother gain
                Ks = P_post(:,:,k) * A' / P_priori(:,:,k+1);
                
                % Smoothed state and covariance update
                x_smooth(:,k) = x_post(:,k) + Ks * (x_smooth(:,k+1) - x_priori(:,k+1));
                P_smooth(:,:,k) = P_post(:,:,k) + Ks * (P_smooth(:,:,k+1) - P_priori(:,:,k+1)) * Ks';
            end
            
            %Find the smoothed output
            y_smooth = C * x_smooth;
        end

        function [Cost, avgCost, Q_pop, R_pop] = optimizeParameters(...
                iterations, mu, lambda, rho, order, lEnd, rEnd, LB, rLA, rLB,...
                A, C, t, y)

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

                [~,y_smooth,~] = KS.simulateKalmanSmoother(A,C,Q,R,y);

                nan_idx = isnan(y);
                y_rec = y;
                y_rec(nan_idx) = y_smooth(nan_idx);

                [rectifiedSpectrum, f, ~] = Utils.computeSpectrum(t, y_rec);
                filteredSpectrum = Utils.computeSpectrum(t, y_smooth);
                Cost(member) = Utils.computeCost(rectifiedSpectrum, filteredSpectrum, f, order);
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

                    [~,y_smooth,~] = KS.simulateKalmanSmoother(A,C,Q,R,y);

                    nan_idx = isnan(y);
                    y_rec = y;
                    y_rec(nan_idx) = y_smooth(nan_idx);

                    [rectifiedSpectrum, f, ~] = Utils.computeSpectrum(t, y_rec);
                    filteredSpectrum = Utils.computeSpectrum(t, y_smooth);
                    Cost(m) = Utils.computeCost(rectifiedSpectrum, filteredSpectrum, f, order);
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