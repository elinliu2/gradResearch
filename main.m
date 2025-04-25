addpath(genpath('.'))

tjs = trajectorySensitivities;

% Sim Constants
SC.T = 1e-3;
SC.t_final = 10;
SC.start = -pi/2;
SC.stop = pi/2;
SC.step = 0.1;
SC.ref = 0;

% Model
% pendulum = PendulumModel(0.2, 0.2);
% plot_all_the_plots(SC, pendulum);


[drone, initialState, ref] = initializeDrone();
[time, stateProgression, drone] = trajectory(drone, SC, initialState, false, ref);

trajSens = droneTrajectorySensitivitiesHybrid.calculateTrajectorySensitivitiesParams(drone, SC, stateProgression, ref);
trajSensTest =    trajectorySensitivitiesTest.calculateTrajectorySensitivitiesParams(drone, SC, initialState, size(stateProgression, 1), ref);
disp("traj sens diff")
disp(norms.supremum(abs(trajSens-trajSensTest)));

% secondTrajSensTestParam = tjs.calculateSecondOrderTrajectorySensitivitiesTestParams(drone, SC, initialState, size(stateProgression, 1), ref);
% secondTrajSensParam = droneTrajectorySensitivitiesHybrid.calculateSecondOrderTrajectorySensitivitiesParams(drone, SC, stateProgression, trajSens, trajSensParam, secondTrajSensTestParam, ref);
% disp("2nd traj sens param diff 1-19 2nd timestamp")
% disp(norms.secondOrderMax(abs(secondTrajSensParam(1:19,:,:,2)-secondTrajSensTestParam(1:19,:,:,2))));
% disp("2nd traj sens param diff 20-23 2nd timestamp")
% disp(norms.secondOrderMax(abs(secondTrajSensParam(20:23,:,:,2)-secondTrajSensTestParam(20:23,:,:,2))));
% disp("2nd traj sens param max diff 2nd timestamp")
% disp(max(abs(secondTrajSensParam(:,:,:,2)-secondTrajSensTestParam(:,:,:,2)), [], "all"));
% disp("2nd traj sens param max diff")
% disp(max(abs(secondTrajSensParam(1:19,:,:,3)-secondTrajSensTestParam(1:19,:,:,3)), [], "all"));

% splotting.trajectoryPlot3d(time, stateProgression(:, 1), stateProgression(:, 2), stateProgression(:, 3), 1, "State Progression");




