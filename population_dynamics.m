% Population Dynamics of Two Interacting Species
% Using 4th Order Runge-Kutta Method
function population_dynamics
    % Part (a): Non-dimensionalization
    % Defining the dimensional parameters
    % Parameters (example values; these can be manipulated for getting different plots)
    r1 = 2; r2 = 1;      % Growth rates
    K1 = 100; K2 = 150;  % Carrying capacities
    b1 = 0.01; b2 = 0.005; % Interaction coefficients
    
    % Non-dimensional parameters
    alpha = r2/r1;       % Ratio of growth rates
    beta1 = b1*K2/r1;    % Scaled interaction coefficient 1
    beta2 = b2*K1/r2;    % Scaled interaction coefficient 2
    kappa = K2/K1;       % Ratio of carrying capacities
    
    % Time settings
    tspan = [0 50];
    dt = 0.01;
    t = tspan(1):dt:tspan(2);
    
    % Part (b) and (c): Phase portraits for different parameter regimes
    % Case 1: Stable coexistence
    params1 = struct('alpha', 1.0, 'beta1', 0.4, 'beta2', 0.4, 'kappa', 1.0);
    plotPhasePortrait(t, params1, 'Case 1: Stable Coexistence', 1);
    
    % Case 2: Competitive exclusion (Species 1 wins)
    params2 = struct('alpha', 0.5, 'beta1', 0.3, 'beta2', 1.5, 'kappa', 1.0);
    plotPhasePortrait(t, params2, 'Case 2: Competitive Exclusion (Species 1 wins)', 2);
    
    % Case 3: Competitive exclusion (Species 2 wins)
    params3 = struct('alpha', 1.2, 'beta1', 1.5, 'beta2', 0.3, 'kappa', 1.0);
    plotPhasePortrait(t, params3, 'Case 3: Competitive Exclusion (Species 2 wins)', 3);
    
    % Case 4: Bistability
    params4 = struct('alpha', 1.0, 'beta1', 1.2, 'beta2', 1.2, 'kappa', 1.0);
    plotPhasePortrait(t, params4, 'Case 4: Bistability', 4);
end

function plotPhasePortrait(t, params, titleStr, figNum)
    % Initial conditions for different trajectories
    N0_sets = [
        [0.2, 0.3];
        [0.8, 0.9];
        [0.4, 0.6];
        [0.9, 0.3];
        [0.3, 0.9]
    ];
    
    figure('Position', [100, 100, 1200, 400]);
    
    % Phase portrait
    subplot(1, 2, 1)
    hold on
    
    % Plotting the trajectories
    for i = 1:size(N0_sets, 1)
        [t, N] = rk4_solver(t, N0_sets(i,:), params);
        plot(N(:,1), N(:,2), 'LineWidth', 1.5)
    end
    
    % Adding the nullclines
    x = linspace(0, 1.5, 100);
    y = linspace(0, 1.5, 100);
    [X, Y] = meshgrid(x, y);
    
    % First nullcline: dx/dt = 0
    contour(X, Y, X.*(1 - X) - params.beta1*X.*Y, [0 0], '--r', 'LineWidth', 1.5)
    % Second nullcline: dy/dt = 0
    contour(X, Y, params.alpha*Y.*(1 - Y/params.kappa) - params.beta2*X.*Y, [0 0], '--b', 'LineWidth', 1.5)
    
    xlabel('x (Scaled Species 1)')
    ylabel('y (Scaled Species 2)')
    title([titleStr ' - Phase Portrait'])
    grid on
    legend('Trajectory 1', 'Trajectory 2', 'Trajectory 3', 'Trajectory 4', 'Trajectory 5', ...
        'dx/dt = 0', 'dy/dt = 0', 'Location', 'best')
    axis([0 1.5 0 1.5])
    
    % Time series
    subplot(1, 2, 2)
    [t, N] = rk4_solver(t, N0_sets(1,:), params);
    plot(t, N(:,1), 'b-', t, N(:,2), 'r-', 'LineWidth', 1.5)
    xlabel('Scaled Time')
    ylabel('Scaled Population')
    title([titleStr ' - Time Evolution'])
    legend('Species 1', 'Species 2')
    grid on
end

function [t, N] = rk4_solver(t, N0, params)
    % Initialize solution arrays
    N = zeros(length(t), 2);
    N(1,:) = N0;
    dt = t(2) - t(1);
    
    % Runge-Kutta 4th Order Method
    for i = 1:length(t)-1
        k1 = dt * f(t(i), N(i,:), params);
        k2 = dt * f(t(i) + dt/2, N(i,:) + k1/2, params);
        k3 = dt * f(t(i) + dt/2, N(i,:) + k2/2, params);
        k4 = dt * f(t(i) + dt, N(i,:) + k3, params);
        N(i+1,:) = N(i,:) + (k1 + 2*k2 + 2*k3 + k4)/6;
    end
end

function dN = f(t, N, params)
    % Non-dimensionalized system of ODEs
    dN = zeros(1,2);
    dN(1) = N(1)*(1 - N(1)) - params.beta1*N(1)*N(2);
    dN(2) = params.alpha*N(2)*(1 - N(2)/params.kappa) - params.beta2*N(1)*N(2);
end