function [z,p, p_iteration, dist_zbar] = findOptimalqGradAscent(drone, SC, z0, ref, dist, zbar)
    zk = zbar;
    zk_prev = z0;
    
    p0 = extractP(drone);
    pk = p0;
    pk_prev = zeros(size(pk));

    n = drone.numStates;
    p = drone.numParams;
    alpha = 5e-1;
    count = 0;
    p_iteration = zeros(p, 100);
    dist_zbar = zeros(1, 100);
    u = 0.5;
    P = eye(n);
    P(22, 22) = 10;
    P(23, 23) = 10;
    Pinv = inv(P);
    vp = 1;
    while sum(abs(pk-pk_prev), 'all') > 1e-15
        count = count + 1;
        [~, stateProgression, drone, stable] = trajectory(drone, SC, zk, false, ref, dist);
        
        [G, DG] = calcGDG(drone, SC, stateProgression, ref);
        vz = transpose(DG*[eye(n); zeros(p, n)]);
        vp = transpose(DG*[zeros(n, p); eye(p)]);
        zk_prev = zk; 
        pk_prev = pk;
        
        pk = pk_prev + alpha*(vz'*(z0-zk_prev) + G)/(vz'*Pinv*vz)*vp;
        drone = assignP(drone, pk);
        zk = (Pinv*vz*vz')/(vz'*Pinv*vz)*(zk_prev-z0) + (Pinv*vz*vp')/(vz'*Pinv*vz)*(pk_prev - pk) + z0 - (Pinv*vz)/(vz'*Pinv*vz)*G;
        p_iteration(:, count) = pk;
        
        [~, ~, drone, stable] = trajectory(drone, SC, zk, false, ref, dist);
        countbk = 0;
        while ~stable && sum(abs(pk-pk_prev), 'all') > 1e-15
            pk = pk_prev + u*(pk-pk_prev);
            drone = assignP(drone, pk);
            zk = zk_prev + u*(zk-zk_prev);
            [~, ~, drone, stable] = trajectory(drone, SC, zk, false, ref, dist);
            countbk = countbk + 1;
        end 
        if countbk > 5 
            u = u^countbk;
        end

        zbar = findClosestzbar(drone, SC, z0, ref, dist, zk_prev);
        dist_zbar(:, count) = (zbar-z0)'*P*(zbar-z0);    
        p_iteration(:, count) = pk;
        disp("outer count")
        disp(count)
        disp("backstepping count")
        disp(countbk)
    end
    z = zk_prev;
    p = pk_prev;
    disp("number of iterations: ")
    disp(count)
end
