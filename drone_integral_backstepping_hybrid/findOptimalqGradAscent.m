function [z,p] = findOptimalqGradAscent(drone, SC, z0, ref, dist, zbar)
    zk = zbar;
    zk_prev = z0;
    
    p0 = extractP(drone);
    pk = p0;

    epsilon = 5e-1;
    stable = true;

    n = drone.numStates;
    p = drone.numParams;

    while sum(abs(zk - zk_prev), "all") > 1e-2 || ~stable
        [~, stateProgression, drone, stable] = trajectory(drone, SC, zk, false, ref, dist);
        if stable
            disp("stable")
            [~, DG] = calcGDG(drone, SC, stateProgression, ref);
            A = [eye(n) (DG*[eye(n); zeros(p, n)])'; DG*[eye(n); zeros(p, n)] 0;];
   
            dzdp = [eye(n) zeros(n, 1)]*(A\[zeros(n, p); -DG*[zeros(n, p); eye(p)]]);
            
            pk_prev = pk;
            d = p0 - pk;
            pdir = ((zk- z0)' * dzdp / norm((zk- z0)' * dzdp))';
            alpha = d'*pdir + sqrt((d'*pdir)^2 - d'*d + epsilon^2);
            pk = pk + alpha * pdir;
            
            drone = assignP(drone, pk);

            zk_prev = zk; 
            zk = [eye(n) zeros(n, 1)]*(A\[z0; DG*[zk_prev;pk_prev]-DG*[zeros(n, p); eye(p)]*pk]);

            u = 1/2;
        else
            disp("unstable")
            zk = zk_prev + u * (zk-zk_prev);
            u = u/2;
        end   
    end
    z = zk;
    p = pk;
end
