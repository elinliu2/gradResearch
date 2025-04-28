addpath(genpath('.'))

tjs = trajectorySensitivities;

% Sim Constants
SC.T = 1e-3;
SC.t_final = 5;
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
disp("traj sens param diff")
disp(norms.supremum(abs(trajSens-trajSensTest)));
disp("traj sens param individual diff")
disp(max(abs(trajSens-trajSensTest), [], "all"));

secondTrajSensParam = droneTrajectorySensitivitiesHybrid.calculateSecondOrderTrajectorySensitivitiesParams(drone, SC, stateProgression, trajSens, ref);
secondTrajSensParamTest =    trajectorySensitivitiesTest.calculateSecondOrderTrajectorySensitivitiesParams(drone, SC, initialState, size(stateProgression, 1), ref);
disp("2nd traj sens param diff")
disp(norms.secondOrderMax(abs(secondTrajSensParam-secondTrajSensParamTest)));
disp("2nd traj sens param individual diff")
disp(max(abs(secondTrajSensParam-secondTrajSensParamTest), [], "all"));

G = norms.infimumOfReciprocalOfOneNorm(trajSens);
tp = trajectorySensitivities.findtp(trajSens);
DG = trajectorySensitivities.calculate_dG(G, trajSens, secondTrajSensParam, tp, drone);
DG_test = trajectorySensitivitiesTest.calculate_dG(drone, SC, initialState, ref);
disp("DG diff")
disp(sum(abs(DG-DG_test), "all"));
% splotting.trajectoryPlot3d(time, stateProgression(:, 1), stateProgression(:, 2), stateProgression(:, 3), 1, "State Progression");




