function [G, DG] = calcGDG(drone, SC, stateProgression, ref)
    trajSens = droneTrajectorySensitivitiesHybrid.calculateTrajectorySensitivitiesParams(drone, SC, stateProgression, ref);
    
    G = norms.infimumOfReciprocalOfOneNorm(trajSens);
    tp = trajectorySensitivities.findtp(trajSens);

    secondTrajSensParam = droneTrajectorySensitivitiesHybrid.calculateSecondOrderTrajectorySensitivitiesParams(drone, SC, stateProgression(1:1+tp, :), trajSens(:, :, 1:1+tp), ref);

    DG = trajectorySensitivities.calculate_dG(G, trajSens, secondTrajSensParam, tp, drone);   
end

