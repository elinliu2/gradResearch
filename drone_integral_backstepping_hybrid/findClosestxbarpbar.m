
tjs = trajectorySensitivities;

% Sim Constants
SC.T = 1e-3;
SC.t_final = 5;

[drone, x0, ref] = initializeDrone();
xk = x0;
xk_prev = 0;

lambdak = 1/2;

n = drone.numStates;
p = drone.numParams-drone.numStates;
stable = true;
u = 1/2;
while sum(abs(xk - xk_prev), "all") > 1e-2
    [time, stateProgression, drone, stable] = trajectory(drone, SC, xk, false, ref);
    disp("xk:")
    disp(string(transpose(xk)))
    pause(0.01)
    if stable
        trajSens = droneTrajectorySensitivitiesHybrid.calculateTrajectorySensitivitiesParams(drone, SC, stateProgression, ref);
        secondTrajSensParam = droneTrajectorySensitivitiesHybrid.calculateSecondOrderTrajectorySensitivitiesParams(drone, SC, stateProgression, trajSens, ref);
    
        G = norms.infimumOfReciprocalOfOneNorm(trajSens);
        tp = trajectorySensitivities.findtp(trajSens);
        DG = trajectorySensitivities.calculate_dG(G, trajSens, secondTrajSensParam, tp, drone);
        
        A = [eye(n)                       transpose(DG*[eye(n); zeros(p+n, n)]); 
             DG*[eye(n); zeros(p+n, n)]     0];
    
        kp1 = A\[x0; DG*[xk; zeros(p+n, 1)] - G];
        xk_prev = xk;
        xk = kp1(1:23);
        lambdak = kp1(24);
        u = 1/2;
    else 
        xk = xk_prev + u * (xk-xk_prev);
        u = u/2;
    end
end

