function [z,p] = findOptimalqGradAscent(drone, SC, z0, ref, dist, zbar)
    zk = zbar;
    zk_prev = z0;
    
    p0 = extractP(drone);
    pk = p0;

    epsilon = 5e-1;
    stable = true;

    n = drone.numStates;
    p = drone.numParams;
    count = 0;
    while sum(abs(zk - zk_prev), "all") > 1e-2 || ~stable
        count = count + 1;
        [~, stateProgression, drone, stable] = trajectory(drone, SC, zk, false, ref, dist);
        if stable
            DG = droneTrajectorySensitivitiesHybrid.calculate_dG(drone, SC, stateProgression, ref);
            vz = DG*[eye(n);zeros(p, n)];
            vp = DG*[zeros(n, p); eye(p)];
            zk_prev = zk; 
            pk_prev = pk;
            
            zk = z0 + transpose(vz)/(vz*transpose(vz))*(-vz*z0+DG*[zk_prev;pk_prev]-vp*pk);
            d = pk - p0;
            alpha = (-d'*vp' + sqrt( (d'*vp')^2 - (d'*d - epsilon^2)*(vp*vp') ))/(vp*vp')*(vz*vz')/(vz*(z0-zk));
            pk = pk + alpha*(vp'*vz)/(vz*vz')*(z0-zk);
            drone = assignP(drone, pk);
            zk = z0 + transpose(vz)/(vz*transpose(vz))*(-vz*z0+DG*[zk_prev;pk_prev]-vp*pk);

            u = 1/2;
        else
            zk = zk_prev + u * (zk-zk_prev);
            u = u/2;
        end   
    end
    z = zk;
    p = pk;
    disp("number of iterations: ")
    disp(count)
end
