addpath('./utils')  
addpath('./pendulum_pd')
addpath('./drone_pid')
addpath('./drone_integral_backstepping/')
addpath('./drone_integral_backstepping_syms/')

tjs = trajectorySensitivities;

% Sim Constants
SC.T = 1e-4;
SC.t_final = 20;
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
m = 1.4; % kg
g = 9.8;
drone = DroneModelIB_syms(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, lambda1, lambda2, lambda3, lambda4, lambda5, lambda6, m, g);

ref.x = 1;
ref.x_dot = 0;
ref.x_ddot = 0;
ref.x_dddot = 0;
ref.x_ddddot = 0;
ref.y = 1;
ref.y_dot = 0;
ref.y_ddot = 0;
ref.y_dddot = 0;
ref.y_ddddot = 0;
ref.z = 100;
ref.z_dot = 0;
ref.z_ddot = 0;
ref.z_dddot = 0;
ref.z_ddddot = 0;
ref.psi = 0;
ref.psi_dot = 0;
ref.psi_ddot = 0;
ref.psi_dddot = 0;
ref.psi_ddddot = 0;

initialState = zeros(23, 1);
for i = 1:12
    initialState(i) = 0;
end
x1 = initialState(1);
x2 = initialState(2);
x3 = 100; % Drone starts 100 m in the air
initialState(3) = x3; 
x4 = initialState(4);
x5 = initialState(5);
x6 = initialState(6);
x7 = initialState(7);
x8 = initialState(8);
x9 = initialState(9);
x10 = initialState(10);
x11 = initialState(11);
x12 = initialState(12);
x16 = initialState(16);
x17 = initialState(17);
x18 = initialState(18);

e7 = ref.z - x3;
e8 = c7*e7 + ref.z_dot + lambda4*x16 - x6;
U_z = m/(cos(x7)*cos(x8)*4)*(g + (1-c7^2+lambda4)*e7 + (c7+c8)*e8 - c7*lambda4*x16);
initialState(19) = U_z;

e9 = ref.x - x1;
e10 = c9*e9 + ref.x_dot + lambda5*x17 - x4;
theta_ref = m/U_z*((1-c9^2  + lambda5)*e9 + (c9 +c10)*e10 + ref.x_ddot - c9*lambda5*x17);
initialState(21) = theta_ref;

e11 = ref.y - x2;
e12 = c11*e11 + ref.y_dot + lambda6*x18 - x5;
phi_ref =   m/U_z*( -(1-c11^2 + lambda6)*e11 - (c11+c12)*e12 - ref.y_ddot + c11*lambda6*x18);
initialState(22) = phi_ref;  

z_ddot = 4*U_z*cos(x7)*cos(x8)/m-g;
e8_dot = c7*(ref.z_dot - x6)+ref.z_ddot+lambda4*e7-z_ddot;
ez     = (1-c7^2+lambda4)*e7+(c7+c8)*e8-c7*lambda4*x16+g;
ez_dot = ((1-c7^2+lambda4)*(ref.z_dot-x6)+(c7+c8)*e8_dot-c7*lambda4*e7);
U_zdot = m/4*( cos(x7)*cos(x8)*ez_dot + ez*(-sin(x7)*cos(x8)*x10-cos(x7)*sin(x8)*x11) )/(cos(x7)*cos(x8))^2;
initialState(20) = U_zdot; 

x_ddot = 4*U_z*(sin(x7)*cos(x8)*cos(x9) + sin(x8)*sin(x9))/(m);
e10_dot= c9*(ref.x_dot - x4) - ref.x_ddot + lambda5*e9 - x_ddot;
ex_dot = (1-c9^2+lambda5)*(ref.x_dot-x4)+(c9+c10)*e10_dot+ref.x_dddot-c9*lambda5*e9;
ex     = (1-c9^2+lambda5)*e9 + (c9+c10)*e10 + ref.x_ddot - c9*lambda5*x17;

y_ddot = 4*U_z*(sin(x7)*cos(x8)*sin(x9) - sin(x8)*cos(x9))/(m);
e12_dot= c11*(ref.y_dot - x5) - ref.y_ddot + lambda6*e11 - y_ddot;
ey_dot = -(1-c11^2+lambda6)*(ref.y_dot-x5)-(c11+c12)*e12_dot-ref.y_dddot+c11*lambda6*e11;
ey     = -(1-c11^2+lambda6)*e11 - (c11+c12)*e12 - ref.y_ddot + c11*lambda6*x18;

theta_ref_dot = m * (U_z*ex_dot - ex*U_zdot)/U_z^2;
initialState(23) = theta_ref_dot;

phi_ref_dot = m * (U_z*ey_dot - ey*U_zdot)/U_z^2;
initialState(24) = phi_ref_dot;

[time, stateProgression, drone] = trajectory(drone, SC, initialState, false, ref);
% trajSens = tjs.calculateTrajectorySensitivities(drone, SC, stateProgression, ref);
% trajSensTest = tjs.calculateTrajectorySensitivitiesTest(drone, SC, initialState, size(stateProgression, 1), ref);
% 
% trajSensParam = tjs.calculateTrajectorySensitivitiesParameters(drone, SC, stateProgression, ref);
% trajSensTestParam = tjs.calculateTrajectorySensitivitiesTestParams(drone, SC, initialState, size(stateProgression, 1), ref);
% 
% secondTrajSensParam = tjs.calculateSecondOrderTrajectorySensitivitiesParams(drone, SC, stateProgression, trajSens, trajSensParam, ref);
% secondTrajSensTestParam = tjs.calculateSecondOrderTrajectorySensitivitiesTestParams(drone, SC, initialState, size(stateProgression, 1), ref);
                        
% splotting.trajectoryPlot3d(time, stateProgression(:,1), stateProgression(:,2), stateProgression(:,3), 1, "State Progression");
d2fdx2 = drone.d2fdx2(drone, initialState, ref);
                
% paramNorm = zeros(1, 50000);
% for x = 1:50000
%     paramNorm(x) = norm(trajSensParam(:,:,x)-trajSensTestParam(:,:,x), 1);
% end
% plot(1:50000, paramNorm)

