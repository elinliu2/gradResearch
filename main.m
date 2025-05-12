addpath(genpath('.'))

% Sim Constants
SC.T = 1e-3;
SC.t_final = 5;

% Model
% pendulum = PendulumModel(0.2, 0.2);
% plot_all_the_plots(SC, pendulum);

[drone, initialState, ref] = initializeDrone();
[time, stateProgression, drone, stable] = trajectory(drone, SC, initialState, false, ref);

% testTrajectorySensitivities(drone, SC, initialState, stateProgression, ref);
% splotting.trajectoryPlot3d(time, stateProgression(:, 1), stateProgression(:, 2), stateProgression(:, 3), 1, "State Progression");

zbar = findClosestzbar(drone, SC, initialState, ref);
z = findOptimalqLinear(drone, SC, initialState, ref, zbar);


