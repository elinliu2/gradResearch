addpath(genpath('.'))

% Sim Constants
SC.T = 1e-3;
SC.t_final = 10;

% Model
% pendulum = PendulumModel(0.2, 0.2);
% plot_all_the_plots(SC, pendulum);

dist.Fx = zeros(1, floor(SC.t_final/SC.T) + 1);
% dist.Fx(1000:2000) = ones(1, 1001);
dist.Fy = zeros(1, floor(SC.t_final/SC.T) + 1);
% dist.Fy(1000:4000) = ones(1, 3001);
dist.Fz = zeros(1, floor(SC.t_final/SC.T) + 1);

dist.Tx = zeros(1, floor(SC.t_final/SC.T) + 1);
dist.Tx(1000:2000) = -60*ones(1, 1001);
dist.Ty = zeros(1, floor(SC.t_final/SC.T) + 1);
dist.Tx(1000:2000) = -51*ones(1, 1001);
dist.Tz = zeros(1, floor(SC.t_final/SC.T) + 1);

[drone, initialState, ref] = initializeDrone();
[time, stateProgression, drone, stable] = trajectory(drone, SC, initialState, false, ref, dist);
% testTrajectorySensitivities(drone, SC, initialState, stateProgression, ref);
splotting.trajectoryPlot3d(time, stateProgression(:, 1), stateProgression(:, 2), stateProgression(:, 3), 1, "State Progression");

% zbar = findClosestzbar(drone, SC, initialState, ref);
% [z, p] = findOptimalqGeometric(drone, SC, initialState, ref, zbar);


