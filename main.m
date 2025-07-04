addpath(genpath('.'))

% Sim Constants
SC.T = 2e-3;
SC.t_final = 10;

% Model
% pendulum = PendulumModel(0.2, 0.2);
% plot_all_the_plots(SC, pendulum);

[drone, initialState, ref, dist] = initializeDrone(SC);

% tic;
% disp("find closest zbar")
% zbar = findClosestzbar(drone, SC, initialState, ref, dist);
% toc;
% tic;
% disp("find optimal p")
% [z, p] = findOptimalqGradAscent(drone, SC, initialState, ref, dist, zbar);
% drone = assignP(drone, p);
% toc;
% distance = distToRoA(drone, SC, initialState, zbar, ref, dist);
% 
% dist.Fx(1000:2000) = -6*ones(1, 1001);
% dist.Fy(1000:2000) = -6*ones(1, 1001);
% 
[time, stateProgression, drone, stable] = trajectory(drone, SC, initialState, false, ref, dist);
% tic;
% trajSens = droneTrajectorySensitivitiesHybrid.calculateTrajectorySensitivitiesParams(drone, SC, stateProgression, ref);
% secondTrajSensParam = droneTrajectorySensitivitiesHybrid.calculateSecondOrderTrajectorySensitivitiesParams(drone, SC, stateProgression, trajSens, ref);
% 
% G = norms.infimumOfReciprocalOfOneNorm(trajSens);
% tp = trajectorySensitivities.findtp(trajSens);
% DG = trajectorySensitivities.calculate_dG(G, trajSens, secondTrajSensParam, tp, drone);   
% DG_new = droneTrajectorySensitivitiesHybrid.calculate_dG(drone, SC, stateProgression, ref);
splotting.trajectoryPlot3d(time, stateProgression(1, :), stateProgression(2, :), stateProgression(3, :), 1, "State Progression");




