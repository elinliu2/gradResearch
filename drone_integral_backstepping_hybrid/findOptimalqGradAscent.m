function [z,p, p_iteration, dist_zbar, G_iteration] = findOptimalqGradAscent(drone, SC, z0, ref, dist, zbar)
    zk = zbar;
    zk_prev = z0;
    
    p0 = extractP(drone);
    pk = p0;
    pk_prev = zeros(size(pk));

    n = drone.numStates;
    p = drone.numParams;
    alpha = 2;
    count = 0;
    p_iteration = zeros(p, 7);
    dist_zbar = zeros(1, 7);
    G_iteration = zeros(1, 7);
    u = 0.5;
    P = eye(n);
    Pinv = inv(P);
    epsilon = 1e-8;
    while sum(abs(pk-pk_prev), 'all') > 1e-8
        count = count + 1;
        [time, stateProgression, drone, ~, index] = trajectory(drone, SC, zk, false, ref, dist);
        [G, DG] = calcGDG(drone, stateProgression, ref, time, index);
        vz = transpose(DG*[eye(n); zeros(p, n)]);
        vp = transpose(DG*[zeros(n, p); eye(p)]);
        zk_prev = zk; 
        pk_prev = pk;
        
        pk = pk_prev + alpha*(vz'*(z0-zk_prev) + G - epsilon)/(vz'*Pinv*vz)*vp;
        drone = assignP(drone, pk);
        zk = (Pinv*vz*vz')/(vz'*Pinv*vz)*(zk_prev-z0) + (Pinv*vz*vp')/(vz'*Pinv*vz)*(pk_prev - pk) + z0 - (Pinv*vz)/(vz'*Pinv*vz)*(G-epsilon);
        p_iteration(:, count) = pk;
        
        [~, stateProgression, drone, stable, ~] = trajectory(drone, SC, zk, false, ref, dist);
        countbk = 0;
        while ~stable && sum(abs(pk-pk_prev), 'all') > 1e-8
            pk = pk_prev + u*(pk-pk_prev);
            drone = assignP(drone, pk);
            zk = zk_prev + u*(zk-zk_prev);
            [~, stateProgression, drone, stable] = trajectory(drone, SC, zk, false, ref, dist);
            countbk = countbk + 1;
        end 
        if countbk > 1 
            u = (u)^countbk;
        end

        % zbar = findClosestzbar(drone, SC, z0, ref, dist, zk);
        % dist_zbar(:, count) = (zbar-z0)'*P*(zbar-z0);    
        % p_iteration(:, count) = pk;
        % 
        % [time, ~, drone, ~, index] = trajectory(drone, SC, z0, false, ref, dist);
        % trajSenstest = droneTrajectorySensitivitiesHybrid.calculateTrajectorySensitivitiesParams(drone, stateProgressiontest, ref, time, index);
        % Gtest = norms.infimumOfReciprocalOfOneNorm(trajSenstest);
        % G_iteration(count) = Gtest;

        disp("outer count")
        disp(count)
        disp("backstepping count")
        disp(countbk)
    end
    if ~stable
        z = zk_prev;
        p = pk_prev;
    else 
        z = zk;
        p = pk;
    end

    disp("number of iterations: ")
    disp(count)
end
