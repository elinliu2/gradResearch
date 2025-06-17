function dist = distToRoA(drone, SC, z0, zbar, ref, dist)
    p0 = extractP(drone);
    qk = [zbar; p0];
    qk_prev = 0;
    
    n = drone.numStates;
    p = drone.numParams;

    u = 1/2;
    
    while sum(abs(qk - qk_prev), "all") > 1e-2
        [~, stateProgression, drone, stable] = trajectory(drone, SC, qk(1:n), false, ref, dist);
        if stable
        
            [G, DG] = calcGDG(drone, SC, stateProgression, ref);

            A = [eye(n+p)  transpose(DG); 
                 DG      0];
        
            kp1 = A\[[z0;p0]; DG*qk - G];
            qk_prev = qk;
            qk = kp1(1:n+p);
            assignP(drone, qk(n+1:n+p));
            u = 1/2;
        else 
            qk = qk_prev + u * (qk-qk_prev);
            u = u/2;
            assignP(drone, qk(n+1:n+p));
        end
    end

    dist = norm([z0;p0]-qk);
    
end