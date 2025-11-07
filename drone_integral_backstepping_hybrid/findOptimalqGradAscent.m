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
    while (sum(abs(pk - pk_prev), "all") > 1e-5)
        count = count + 1;
        [~, stateProgression, drone, ~] = trajectory(drone, SC, zk, false, ref, dist);
        
        [G, DG] = calcGDG(drone, SC, stateProgression, ref);
        vz = transpose(DG*[eye(n); zeros(p, n)]);
        vp = transpose(DG*[zeros(n, p); eye(p)]);
        zk_prev = zk; 
        pk_prev = pk;
        
        pk = pk_prev + alpha*(vp*vz')/(vz'*vz)*(z0-zk);
        drone = assignP(drone, pk);
        zk = z0 + (1 + alpha*(vp')*vp/(vz'*vz))*(vz*vz')/(vz'*vz)*(zk-z0) - vz/(vz'*vz)*G;
        p_iteration(:, count) = pk;
        u = 1/2;

        [~, ~, drone, stable] = trajectory(drone, SC, zk, false, ref, dist);
        while ~stable
            pk = pk_prev + u*(pk-pk_prev);
            drone = assignP(drone, pk);
            zk = zk_prev + u*(zk-zk_prev);
            u = u/2;
            [~, ~, drone, stable] = trajectory(drone, SC, zk, false, ref, dist);
        end 

        zbar = findClosestzbar(drone, SC, z0, ref, dist, zk_prev);
        dist_zbar(:, count) = 1/2*norm(zbar-z0)^2;    
        p_iteration(:, count) = pk;
        disp(count)
    end
    z = zk_prev;
    p = pk_prev;
    disp("number of iterations: ")
    disp(count)
end
