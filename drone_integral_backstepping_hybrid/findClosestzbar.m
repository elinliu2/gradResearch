function zbar = findClosestzbar(drone, SC, z0, ref, dist, zbar_guess)
    zk = zbar_guess;
    zk_prev = zbar_guess;
    
    n = drone.numStates;
    p = drone.numParams;

    u = 1/2;

    P = eye(n);
    Pinv = inv(P);

    count = 0;
    epsilon = 1e-7;
    while sum(abs(zk - zk_prev), "all") > 1e-8 || count == 0
        count = count + 1;
        [time, stateProgression, drone, stable, index] = trajectory(drone, SC, zk, false, ref, dist);
        if stable
            [G, DG] = calcGDG(drone, stateProgression, ref, time, index);
            vz = transpose(DG*[eye(n); zeros(p, n)]);
            zk_prev = zk;
            zk = z0 + (Pinv*vz*vz')/(vz'*Pinv*vz)*(zk - z0) - Pinv*vz/(vz'*Pinv*vz)*(G-epsilon);
            u = 1/2;
        else 
            zk = zk_prev + u * (zk-zk_prev);
            u = u/2;
        end
        disp(count)
    end

    disp("number of iterations: ")
    disp(count)
    if ~stable
        zbar = zk_prev;
    else 
        zbar = zk;
    end
end
