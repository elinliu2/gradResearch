function zbar = findClosestzbar(drone, SC, z0, ref, dist)
    zk = z0;
    zk_prev = 0;
    
    n = drone.numStates;
    p = drone.numParams;

    u = 1/2;
    count = 0;
    while sum(abs(zk - zk_prev), "all") > 1e-2
        count = count + 1;
        [~, stateProgression, drone, stable] = trajectory(drone, SC, zk, false, ref, dist);
        if stable
        
            [G, DG] = calcGDG(drone, SC, stateProgression, ref);

            A = [eye(n)                       transpose(DG*[eye(n); zeros(p, n)]); 
                 DG*[eye(n); zeros(p, n)]     0];
        
            kp1 = A\[z0; DG*[zk; zeros(p, 1)] - G];
            zk_prev = zk;
            zk = kp1(1:23);
            u = 1/2;
        else 
            zk = zk_prev + u * (zk-zk_prev);
            u = u/2;
        end
    end

    disp("number of iterations: ")
    disp(count)
    zbar = zk;
end
