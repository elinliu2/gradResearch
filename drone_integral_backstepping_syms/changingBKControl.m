addpath('./../utils/')
tjs = trajectorySensitivities;

% Sim Constants
SC.T = 1;
SC.t_final = 6;
SC.start = -pi/2;
SC.stop = pi/2;
SC.step = 0.1;
SC.ref = 0;

% Model
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

model = DroneModelIB_syms(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, lambda1, lambda2, lambda3, lambda4, lambda5, lambda6, m, g);

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

initialState = zeros(19, 1);
for i = 1:12
    initialState(i) = 0;
end

x3 = 100; % Drone starts 100 m in the air
initialState(3) = x3; 
x6 = initialState(6);
x7 = initialState(7);
x8 = initialState(8);
x16 = initialState(16);

e7 = ref.z - x3;
e8 = c7*e7 + ref.z_dot + lambda4*x16 - x6;
U_z = m/(cos(x7)*cos(x8)*4)*(g + (1-c7^2+lambda4)*e7 + (c7+c8)*e8 - c7*lambda4*x16);
initialState(19) = U_z;

numIterations = 150;
param_tracker = zeros(model.numParams, numIterations);
G_tracker = zeros(1, numIterations);

for counter=1:numIterations
    param_tracker(:, counter) = [
        model.c1; model.c2; model.c3; model.c4; model.c5; model.c6; 
        model.c7; model.c8; model.c9; model.c10; model.c11; model.c12; 
        model.lambda1; model.lambda2; model.lambda3; 
        model.lambda4; model.lambda5; model.lambda6];

    [~, stateProgression, ~] = trajectory(model, SC, initialState, false, ref);
    trajSens = tjs.calculateTrajectorySensitivities(model, SC, stateProgression, ref);
    tp = tjs.findtp(trajSens);

    G = norms.infimumOfReciprocalOfElementOneNorm(trajSens);
    trajSensParam = tjs.calculateTrajectorySensitivitiesParameters(model, SC, stateProgression, ref);
    secondTrajSensParam = tjs.calculateSecondOrderTrajectorySensitivitiesParams(model, SC, stateProgression, trajSens, trajSensParam, ref);

    
    dGp = tjs.calculate_dG(G, trajSens, secondTrajSensParam, tp, model);

    scale = 0.01;
    normdGp = sqrt(dGp(1)^2 + dGp(2)^2);

    model.c1 = model.c1 + dGp(1)/normdGp*scale;
    model.c2 = model.c2 + dGp(2)/normdGp*scale;
    model.c3 = model.c3 + dGp(3)/normdGp*scale;
    model.c4 = model.c4 + dGp(4)/normdGp*scale;
    model.c5 = model.c5 + dGp(5)/normdGp*scale;
    model.c6 = model.c6 + dGp(6)/normdGp*scale;
    model.c7 = model.c7 + dGp(7)/normdGp*scale;
    model.c8 = model.c8 + dGp(8)/normdGp*scale;
    model.c9 = model.c9 + dGp(9)/normdGp*scale;
    model.c10 = model.c10 + dGp(10)/normdGp*scale;
    model.c11 = model.c11 + dGp(11)/normdGp*scale;
    model.c12 = model.c12 + dGp(12)/normdGp*scale;

    model.lambda1 = model.lambda1 + + dGp(13)/normdGp*scale;
    model.lambda2 = model.lambda2 + + dGp(14)/normdGp*scale;
    model.lambda3 = model.lambda3 + + dGp(15)/normdGp*scale;
    model.lambda4 = model.lambda4 + + dGp(16)/normdGp*scale;
    model.lambda5 = model.lambda5 + + dGp(17)/normdGp*scale;
    model.lambda6 = model.lambda6 + + dGp(18)/normdGp*scale;

    G_tracker(counter) = G;
end

% for i = 1:18
%     splotting.uplup(1:counter, param_tracker(1:counter), i, "param over iteration", string(i));
% end

splotting.uplup(1:counter, G_tracker(1:counter), 19, "G over iteration", "G");