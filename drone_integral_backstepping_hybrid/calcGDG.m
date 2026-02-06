function [G, DG] = calcGDG(drone, stateProgression, ref, time, index)
    % trajSens = droneTrajectorySensitivitiesHybrid.calculateTrajectorySensitivities(drone, stateProgression, ref, time, index);
    trajSensParam = droneTrajectorySensitivitiesHybrid.calculateTrajectorySensitivitiesParams(drone, stateProgression, ref, time, index);
    trajSens = trajSensParam(:, 1:23, :);
    
    G = norms.infimumOfReciprocalOfOneNorm(trajSens);
    tp = trajectorySensitivities.findtp(trajSens);

    secondTrajSensParam = ...
        droneTrajectorySensitivitiesHybrid.calculateSecondOrderTrajectorySensitivitiesParams( ...
        drone, stateProgression(1:1+tp, :), trajSensParam(:, :, 1:1+tp), ref,  time, tp);

    DG = trajectorySensitivities.calculate_dG(G, trajSens, secondTrajSensParam, tp, drone);   
end

