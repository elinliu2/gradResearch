function [z,p] = findOptimalq(drone, SC, z0, ref, zbar)
    zk = zbar;
    zk_prev = 0;
    
    n = drone.numStates;
    p = drone.numParams;
    p0 = extractP(drone);
    pk = p0;
    epsilon = 1e-1;
    
    while sum(abs(zk - zk_prev), "all") > 1e-2
        [~, stateProgression, drone, stable] = trajectory(drone, SC, zk, false, ref);
        if stable
           [~, DG] = calcGDG(drone, SC, stateProgression, ref);
            
            A = [eye(n)                    zeros(n, p)                transpose(DG*[eye(n); zeros(p, n)]); 
                zeros(p, n)                zeros(p, p)                transpose(DG*[zeros(n, p); eye(p)]);
                DG*[eye(n); zeros(p, n)]   DG*[zeros(n, p); eye(p)]   0];
            % b = [z0; zeros(p, 1); DG*[zk; pk]];

            w = null(A);
            
            maxDist = 0;

            for i = 1:size(w, 2)
                vp = pk;
                wp = w(24:41, i);
    
                a2 = wp'*wp;
                a1 = vp'*wp + wp'*vp - p0'*wp - wp'*p0;
                a0 = vp'*vp-p0'*vp-vp'*p0+p0'*p0-epsilon^2;
    
                t = roots([a2 a1 a0]);
    
                z1 = zk + t(1) * w(1:23, i);
                z2 = zk + t(2) * w(1:23, i);
    
                zdist1 = norm(abs(z1(1:n) - z0));
                zdist2 = norm(abs(z2(1:n) - z0));

                if (zdist1 > maxDist)
                    zk = z1(1:n);
                    zk_prev = zk;
                    pk = vp + t(1) * w(24:41, i);
                    drone = assignP(drone, pk);
                    maxDist = zdist1;
                elseif (zdist2 > maxDist)
                    zk = z2(1:n);
                    zk_prev = zk;
                    pk = vp + t(2) * w(24:41, i);
                    drone = assignP(drone, pk);
                    maxDist = zdist2;
                end
            end

            u = 1/2;
        else
            disp("unstable")
            zk = zk_prev + u * (zk-zk_prev);
            u = u/2;
        end
    end

    z = zk;
end

function p = extractP(drone)
    p = zeros(drone.numParams, 1);
    p(1) = drone.c1;
    p(2) = drone.c2;
    p(3) = drone.c3;
    p(4) = drone.c4;
    p(5) = drone.c5;
    p(6) = drone.c6;
    p(7) = drone.c7;
    p(8) = drone.c8;
    p(9) = drone.c9;
    p(10) = drone.c10;
    p(11) = drone.c11;
    p(12) = drone.c12;

    p(13) = drone.lambda1;
    p(14) = drone.lambda2;
    p(15) = drone.lambda3;
    p(16) = drone.lambda4;
    p(17) = drone.lambda5;
    p(18) = drone.lambda6;
end

function drone = assignP(drone, p)
    drone.c1      = p(1);
    drone.c2      = p(2);
    drone.c3      = p(3);
    drone.c4      = p(4);
    drone.c5      = p(5);
    drone.c6      = p(6);
    drone.c7      = p(7);
    drone.c8      = p(8);
    drone.c9      = p(9);
    drone.c10     = p(10);
    drone.c11     = p(11);
    drone.c12     = p(12);
    
    drone.lambda1 = p(13);
    drone.lambda2 = p(14);
    drone.lambda3 = p(15);
    drone.lambda4 = p(16);
    drone.lambda5 = p(17);
    drone.lambda6 = p(18);

end