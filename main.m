addpath(genpath('.'))
warning('off', 'all')

% Sim Constants
SC.T = 1e-3;
SC.t_final = 60;

% Model
% pendulum = PendulumModel(0.2, 0.2);
% plot_all_the_plots(SC, pendulum);

[drone, z0, ref, dist] = initializeDrone(SC);

tic;
disp("find closest zbar")
zbar = findClosestzbar(drone, SC, z0, ref, dist, z0);
toc;

tic;
disp("find optimal p")
[z, p, p_iteration] = findOptimalqGradAscent(drone, SC, z0, ref, dist, zbar);
drone = assignP(drone, p);
toc;

t = 1;
dist.Fx(1:1+t) = -6*ones(1, 1+t);
dist.Fy(1:1+t) = 6*ones(1, 1+t);

[time, stateProgression, drone, stable] = trajectory(drone, SC, z0, false, ref, dist);
% trajSens = droneTrajectorySensitivitiesHybrid.calculateTrajectorySensitivitiesParams(drone, SC, stateProgression, ref);
% secondTrajSensParam = droneTrajectorySensitivitiesHybrid.calculateSecondOrderTrajectorySensitivitiesParams(drone, SC, stateProgression, trajSens, ref);
% G = norms.infimumOfReciprocalOfOneNorm(trajSens);
splotting.trajectoryPlot3d(time, stateProgression(:, 1), stateProgression(:, 2), stateProgression(:, 3), 1, "State Progression");


figure(2)
for i = 1:18
    subplot(3, 6, i);
    p_iterate = p_iteration(i, :);
    plot(p_iterate)
end


