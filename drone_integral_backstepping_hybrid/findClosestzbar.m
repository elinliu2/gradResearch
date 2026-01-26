function zbar = findClosestzbar(drone, SC, z0, ref, dist, zbar_guess)
    zk = zbar_guess;
    zk_prev = zbar_guess;
    
    n = drone.numStates;
    p = drone.numParams;

    u = 1/2;

    P = eye(n);
    Pinv = inv(P);

    count = 0;
    epsilon = 1e-8;
    while sum(abs(zk - zk_prev), "all") > 1e-8 || count == 0
        count = count + 1;
        [~, stateProgression, drone, stable] = trajectory(drone, SC, zk, false, ref, dist);
        if stable
            [G, DG] = calcGDG(drone, SC, stateProgression, ref);
            vz = transpose(DG*[eye(n); zeros(p, n)]);
            zk_prev = zk;
            zk = z0 + (Pinv*vz*vz')/(vz'*Pinv*vz)*(zk - z0) - Pinv*vz/(vz'*Pinv*vz)*(G-epsilon);
            u = 1/2;
        else 
            zk = zk_prev + u * (zk-zk_prev);
            u = u/2;
        end
    end

    % disp("number of iterations: ")
    % disp(count)
    if ~stable
        zbar = zk_prev;
    else 
        zbar = zk;
    end
end
