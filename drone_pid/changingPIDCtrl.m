addpath('./../utils/')
tjs = trajectorySensitivities;

% Sim Constants
SC.T = 1e-3;
SC.t_final = pi;
SC.start = -pi/2;
SC.stop = pi/2;
SC.step = 0.1;
SC.ref = 0;

% Model
kp1 = 5;
ki1 = 0.5;
kd1 = 10;

kp2 = 5;
ki2 = 1;
kd2 = 5;

kp3 = 5;
ki3 = 0.5;
kd3 = 10;

kp4 = 20;
ki4 = 3;
kd4 = 20;

drone = DroneModel(kp1, ki1, kd1, kp2, ki2, kd2, kp3, ki3, kd3, kp4, ki4, kd4);

initialState = zeros(16, 1);
initialState(3) = 100; % Drone starts 100 m in the air
initialState(8) = pi/6;

numIterations = 500;
parameter_tracker = zeros(drone.numParams, numIterations);
G_tracker = zeros(1, numIterations);
for counter=1:numIterations
    parameters = [drone.kp1, drone.ki1, drone.kd1, drone.kp2, drone.ki2, drone.kd2, drone.kp3, drone.ki3, drone.kd3, drone.kp4, drone.ki4, drone.kd4];
    for index=1:drone.numParams
        parameter_tracker(index, counter) = parameters(index);
    end

    [~, stateProgression, ~] = trajectory(drone, SC, initialState, false);
    trajSens = tjs.calculateTrajectorySensitivities(drone, SC, stateProgression);
    G = norms.infimumOfReciprocalOfElementOneNorm(trajSens);
    trajSensParam = tjs.calculateTrajectorySensitivitiesParameters(drone, SC, stateProgression);
    secondTrajSensParam = tjs.calculateSecondOrderTrajectorySensitivitiesParams(drone, SC, stateProgression, trajSens, trajSensParam);

    tp = tjs.findtp(trajSens);
    dGp = tjs.calculate_dG(G, trajSens, secondTrajSensParam, tp, drone);

    scale = 0.1;
    normdGp = sqrt(dGp(1)^2 + dGp(2)^2);

    for index=1:drone.numParams
        parameters(index) = parameters(index) + dGp(index)/normdGp*scale;
    end

    drone.kp1 = parameters(1);
    drone.ki1 = parameters(2);
    drone.kd1 = parameters(3);
    drone.kp2 = parameters(4);
    drone.ki2 = parameters(5);
    drone.kd2 = parameters(6);
    drone.kp3 = parameters(7);
    drone.ki3 = parameters(8);
    drone.kd3 = parameters(9);
    drone.kp4 = parameters(10);
    drone.ki4 = parameters(11);
    drone.kd4 = parameters(12);
    
    G_tracker(counter) = G;
end

splotting.uplup(1:counter, G_tracker(1:counter), 4, "G over iteration", "G");