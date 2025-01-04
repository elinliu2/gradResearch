addpath('./utils')  
addpath('./pendulum_pd')
addpath('./drone_pid')

% Sim Constants
SC.T = 1e-3;
SC.t_final = 60;
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

a = 5; 
b = 0.01;
c = 10;

kp1 = a;
ki1 = b;
kd1 = c;

kp2 = a;
ki2 = b;
kd2 = c;

kp3 = a;
ki3 = b;
kd3 = c;

kp4 = 5;
ki4 = 0.05;
kd4 = 5;

drone = DroneModel(kp1, ki1, kd1, kp2, ki2, kd2, kp3, ki3, kd3, kp4, ki4, kd4);
initialState = zeros(16, 1);
initialState(3) = 100; % Drone starts 100 m in the air

ref = [0.01; 0.01; 101];
[time, stateProgression, drone] = trajectory(drone, SC, initialState, ref, false);
splotting.trajectoryPlot3d(time, stateProgression(:,1), stateProgression(:,2), stateProgression(:,3), 1, "State Progression");
