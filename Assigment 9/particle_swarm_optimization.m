clc
clear
close all

% Start Timing
tic;

% Parameters
N = 25;  % Number of particles
n = 3;   % Number of dimensions (independent variables)
phi1_max = 2.0;
phi2_max = 2.0;
vmax = 3;
sigma = 3; % Neighborhood size
w = 0.4;  % Inertia weighting
max_iter = 30; % Maximum number of iterations

% Gain constraint
GMh_min = 0;
GMh_max = 2;

% Parameters initial limits
initial_limit = 15;

% Initialize the swarm
x = unifrnd(-initial_limit, initial_limit, N, n);

% Gain Gmh between 0 and 2
x(:,3) = (x(:,3)+initial_limit)/initial_limit;

v = zeros(N, n);
b = x;

% Evaluate initial population inverting the function to maximize result
fx = arrayfun(@(i) -simFox(x(i, :)), 1:N)';
fb = fx;

% Prepare for dynamic plot
best_velocity = zeros(1, max_iter);
best_position = zeros(max_iter, n);
figure(1);
plotHandle = plot(nan, nan);
hold on;
grid on;
xlabel('Iteration');
ylabel('Best Velocity');
title('Best Velocity Over Iterations');


figure(2);
scatterHandle = scatter3(nan,nan,nan);
hold on;
grid on;
xlabel('PEAhip');
ylabel('AEAhip');
zlabel('GMh');
title('Particles Position Over Iterations');
trails = zeros(N, n, max_iter);
trails(:,:,1) = x;
xlim([-initial_limit initial_limit]);
ylim([-initial_limit initial_limit]);
zlim([0 2]);

% PSO loop
for iter = 1:max_iter
    for i = 1:N
        % Find the neighborhood
        distances = sum((x - x(i, :)).^2, 2);
        [~, sorted_indices] = sort(distances);
        H_i = sorted_indices(2:sigma+1); % Exclude the particle itself

        % Find the best position in the neighborhood
        [~, best_neighbor_idx] = min(fx(H_i));
        h_i = x(H_i(best_neighbor_idx), :);

        % Generate random vectors
        phi1 = unifrnd(0, phi1_max, 1, n);
        phi2 = unifrnd(0, phi2_max, 1, n);

        % Update velocity
        v(i, :) = w * v(i, :) + phi1 .* (b(i, :) - x(i, :)) + phi2 .* (h_i - x(i, :));

        % Apply velocity constraint
        if norm(v(i, :)) > vmax
            v(i, :) = (v(i, :) / norm(v(i, :))) * vmax;
        end

        % Update position
        x(i, :) = x(i, :) + v(i, :);

        % Apply position constraint
        x(i, 3) = max(min(x(i, 3), GMh_max), GMh_min);

        % Evaluate new position inverting the function to maximization
        fx(i) = -simFox(x(i, :));

        % Update the best position
        if fx(i) < fb(i)
            b(i, :) = x(i, :);
            fb(i) = fx(i);
        
        end
        % Store trail for scatter
        trails(i,:,iter) = x(i, :);
    end
    % Store the best velocity (negated to get the actual max velocity)
    [best_velocity_neg,best_velocity_index] = min(fb);
    best_velocity(iter) = -best_velocity_neg;
    best_position(iter,:) = b(best_velocity_index,:);

    % Update the plot dynamically
    set(scatterHandle, 'XData', x(:,1), 'YData', x(:,2), 'ZData', x(:,3));
    set(plotHandle, 'XData', 1:iter, 'YData', best_velocity(1:iter));
    drawnow;

    % Display the best objective value found so far
    disp(['Iteration ', num2str(iter), ': Best Objective = ', num2str(best_velocity(iter)), ' Position =',num2str(best_position(iter,:))]);
end
% Plot trails after optimization
figure(2);
hold on;
for i = 1:N
    plot3(squeeze(trails(i,1,:)), squeeze(trails(i,2,:)), squeeze(trails(i,3,:)));
end
hold off;

% Stop Timing
runtime = toc
