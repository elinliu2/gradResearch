addpath(genpath('.'))

% Sim Constants
SC.T = 1e-3;
SC.t_final = 10;

% Model
% pendulum = PendulumModel(0.2, 0.2);
% plot_all_the_plots(SC, pendulum);

[drone, initialState, ref, dist] = initializeDrone(SC);

zbar = findClosestzbar(drone, SC, initialState, ref, dist);
[z, p] = findOptimalqGradAscent(drone, SC, initialState, ref, dist, zbar);
drone = assignP(drone, p);

% dist.Fx(1000:2000) = -5.9*ones(1, 1001);
% dist.Fy(1000:2000) = -5.9*ones(1, 1001);
% 
% [time, stateProgression, drone, stable] = trajectory(drone, SC, initialState, false, ref, dist);
% splotting.trajectoryPlot3d(time, stateProgression(:, 1), stateProgression(:, 2), stateProgression(:, 3), 1, "State Progression");




