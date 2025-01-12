addpath('./utils')  
addpath('./pendulum_pd')
addpath('./drone_pid')

% Sim Constants
SC.T = 1e-3;
SC.t_final = pi;
SC.start = -2*pi;
SC.stop = 2*pi;
SC.step = 0.1;
SC.ref = 0;

% Model
pendulum = PendulumModel(0.2, 0.2);
plot_all_the_plots(SC, pendulum);

% initialState = [-2, -0.5];
% [time, stateProgression, pendulum] = trajectory(pendulum, SC, initialState, false);
% splotting.StatePlot(time, stateProgression(:,1), stateProgression(:,2), 1, "State Progression");
% splotting.uplup(time, stateProgression(:,1), 2, "x1 over time", "x1");
% splotting.uplup(time, stateProgression(:,2), 3, "x2 over time", "x2");

% kp1 = 5;
% ki1 = 0.5;
% kd1 = 10;
% 
% kp2 = 5;
% ki2 = 1;
% kd2 = 5;
% 
% kp3 = 5;
% ki3 = 0.5;
% kd3 = 10;
% 
% kp4 = 20;
% ki4 = 3;
% kd4 = 20;

% drone = DroneModel(kp1, ki1, kd1, kp2, ki2, kd2, kp3, ki3, kd3, kp4, ki4, kd4);
% initialState = zeros(16, 1);
% initialState(3) = 100; % Drone starts 100 m in the air
% initialState(8) = pi/6;
% 
% [time, stateProgression, drone] = trajectory(drone, SC, initialState, false);
% splotting.trajectoryPlot3d(time, stateProgression(:,1), stateProgression(:,2), stateProgression(:,3), 1, "State Progression");
% plotDroneXYZ(time, stateProgression, 2);
% plotDroneAngles(time, stateProgression, 3);

% plot_all_the_plots(SC, drone);



