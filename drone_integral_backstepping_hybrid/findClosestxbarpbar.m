
tjs = trajectorySensitivities;

% Sim Constants
SC.T = 1e-3;
SC.t_final = 5;

[drone, x0, ref] = initializeDrone();
xk = x0;
lambdak = 1/2;

n = drone.numStates;
p = drone.numParams-drone.numStates;

for i = 1:2
    [time, stateProgression, drone, stable] = trajectory(drone, SC, xk, false, ref);
    
    trajSens = droneTrajectorySensitivitiesHybrid.calculateTrajectorySensitivitiesParams(drone, SC, stateProgression, ref);
    secondTrajSensParam = droneTrajectorySensitivitiesHybrid.calculateSecondOrderTrajectorySensitivitiesParams(drone, SC, stateProgression, trajSens, ref);

    G = norms.infimumOfReciprocalOfOneNorm(trajSens);
    tp = trajectorySensitivities.findtp(trajSens);
    DG = trajectorySensitivities.calculate_dG(G, trajSens, secondTrajSensParam, tp, drone);
    
    A = [eye(n)                       transpose(DG*[eye(n); zeros(p+n, n)]); 
         DG*[eye(n); zeros(p+n, n)]     0];

    kp1 = A\[x0; DG*[xk; zeros(p+n, 1)] - G];
    xk = kp1(1:23);
    lambdak = kp1(24);
end

