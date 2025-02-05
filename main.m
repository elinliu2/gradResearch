addpath('./utils')  
addpath('./pendulum_pd')
addpath('./drone_pid')
addpath('./drone_integral_backstepping/')

% Sim Constants
SC.T = 1e-3;
SC.t_final = 200;
SC.start = -pi/2;
SC.stop = pi/2;
SC.step = 0.1;
SC.ref = 0;

% Model
% pendulum = PendulumModel(0.2, 0.2);
% plot_all_the_plots(SC, pendulum);

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
% 
% drone = DroneModel(kp1, ki1, kd1, kp2, ki2, kd2, kp3, ki3, kd3, kp4, ki4, kd4);
% initialState = zeros(16, 1);
% initialState(3) = 100; % Drone starts 100 m in the air
% 
% [time, stateProgression, drone] = trajectory(drone, SC, initialState, false);
% splotting.trajectoryPlot3d(time, stateProgression(:,1), stateProgression(:,2), stateProgression(:,3), 1, "State Progression");
% % splotting.plotDroneXYZ(time, stateProgression, 2);
% splotting.plotDroneAngles(time, stateProgression, 3);
c1 = 10.5;
c2 = 2;
c3 = 10;
c4 = 2;
c5 = 2;
c6 = 2;
c7 = 4;
c8 = 2;
c9 = 1;
c10 = 0.5;
c11 = 1;
c12 = 0.5;
lambda1 = 1e-1;
lambda2 = 1e-1;
lambda3 = 1e-1;
lambda4 = 1e-1;
lambda5 = 1e-1;
lambda6 = 1e-1;
drone = DroneModelIB(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, lambda1, lambda2, lambda3, lambda4, lambda5, lambda6);
initialState = zeros(18, 1);
for i = 1:12
    initialState(i) = 0.9;
end
initialState(3) = 100; % Drone starts 100 m in the air

[time, stateProgression, drone] = trajectory(drone, SC, initialState, false);
splotting.trajectoryPlot3d(time, stateProgression(:,1), stateProgression(:,2), stateProgression(:,3), 1, "State Progression");

% plot_all_the_plots(SC, drone);



