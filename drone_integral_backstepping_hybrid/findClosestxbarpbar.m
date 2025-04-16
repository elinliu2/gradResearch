addpath('./../utils/')  
tjs = trajectorySensitivities;

% Sim Constants
SC.T = 1e-3;
SC.t_final = 10;

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
drone = droneModelIBHybrid(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, lambda1, lambda2, lambda3, lambda4, lambda5, lambda6, m, g);

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

x0 = zeros(drone.numStates, 1);
for i = 1:12
    x0(i) = 0;
end
x3 = 100; % Drone starts 100 m in the air
x0(3) = x3; 
x6 = x0(6);
x7 = x0(7);
x8 = x0(8);
x16 = x0(16);

e7 = ref.z - x3;
e8 = c7*e7 + ref.z_dot + lambda4*x16 - x6;
U_z = m/(cos(x7)*cos(x8)*4)*(g + (1-c7^2+lambda4)*e7 + (c7+c8)*e8 - c7*lambda4*x16);
x0(19) = U_z;

xk = x0;
lambdak = 1/2;

n = drone.numStates;
p = drone.numParams-drone.numStates;

for i = 1:2
    [time, stateProgression, drone] = trajectory(drone, SC, xk, false, ref);
    
    trajSens = droneTrajectorySensitivitiesHybrid.calculateTrajectorySensitivities(drone, SC, stateProgression, ref);
    trajSensParam = droneTrajectorySensitivitiesHybrid.calculateTrajectorySensitivitiesParameters(drone, SC, stateProgression, ref);
    secondTrajSensParam = droneTrajectorySensitivitiesHybrid.calculateSecondOrderTrajectorySensitivitiesParams(drone, SC, stateProgression, trajSens, trajSensParam, ref);
    
    G = norms.infimumOfReciprocalOfElementOneNorm(trajSens);
    tp = tjs.findtp(trajSens);
    dGp = tjs.calculate_dG(G, trajSens, secondTrajSensParam, tp, drone)*1e15;
    
    A = [eye(n)                       transpose(dGp*[eye(n); zeros(p, n)]); 
         dGp*[eye(n); zeros(p, n)]    0];

    kp1 = A\[-x0; -dGp*[xk; zeros(p, 1)] + G];
    xk = kp1(1:23);
    lambdak = kp1(24);
    SC.t_final = 20e-3;
end

splotting.trajectoryPlot3d(time, stateProgression(:, 1), stateProgression(:, 2), stateProgression(:, 3), 1, "State Progression");

