function [G, DG] = calcGDG(drone, SC, stateProgression, ref)
    trajSens = droneTrajectorySensitivitiesHybrid.calculateTrajectorySensitivitiesParams(drone, SC, stateProgression, ref);
    secondTrajSensParam = droneTrajectorySensitivitiesHybrid.calculateSecondOrderTrajectorySensitivitiesParams(drone, SC, stateProgression, trajSens, ref);

    G = norms.infimumOfReciprocalOfOneNorm(trajSens);
    tp = trajectorySensitivities.findtp(trajSens);
    DG = trajectorySensitivities.calculate_dG(G, trajSens, secondTrajSensParam, tp, drone);   
end

