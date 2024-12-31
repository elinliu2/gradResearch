addpath('./utils')  
addpath('./pendulum_pd')
addpath('./drone_pid')

% Sim Constants
SC.T = 1e-2;
SC.t_final = 15;
SC.start = -1;
SC.stop = 1;
SC.step = 0.5;
SC.ref = 0;

% % Model
% pendulum = PendulumModel(0.2, 0.2);
% plot_all_the_plots(SC, pendulum);

% initialState = [-2, -0.5];
% [time, stateProgression, pendulum] = trajectory(pendulum, SC, initialState, false);
% splotting.StatePlot(time, stateProgression(:,1), stateProgression(:,2), 1, "State Progression");
% splotting.uplup(time, stateProgression(:,1), 2, "x1 over time", "x1");
% splotting.uplup(time, stateProgression(:,2), 3, "x2 over time", "x2");
kp1 = 5;
ki1 = 0.5;
kd1 = 5;

kp2 = 5;
ki2 = 0.5;
kd2 = 5;

kp3 = 5;
ki3 = 0.5;
kd3 = 5;

kp4 = 10;
ki4 = 0.5;
kd4 = 5;

drone = DroneModel(kp1, ki1, kd1, kp2, ki2, kd2, kp3, ki3, kd3, kp4, ki4, kd4);
initialState = zeros(15, 1);
initialState(3) = 100; % Drone starts 100 m in the air
[time, stateProgression, drone] = trajectory(drone, SC, initialState, false);
splotting.trajectoryPlot3d(time, stateProgression(:,1), stateProgression(:,2), stateProgression(:,3), 1, "State Progression");


