function [G, DG] = calcGDG(drone, stateProgression, ref, time, index)
    trajSens = droneTrajectorySensitivitiesHybrid.calculateTrajectorySensitivitiesParams(drone, stateProgression, ref, time, index);
    
    G = norms.infimumOfReciprocalOfOneNorm(trajSens);
    tp = trajectorySensitivities.findtp(trajSens);

    secondTrajSensParam = ...
        droneTrajectorySensitivitiesHybrid.calculateSecondOrderTrajectorySensitivitiesParams( ...
        drone, stateProgression(1:1+tp, :), trajSens(:, :, 1:1+tp), ref,  time, tp);

    DG = trajectorySensitivities.calculate_dG(G, trajSens, secondTrajSensParam, tp, drone);   
end

