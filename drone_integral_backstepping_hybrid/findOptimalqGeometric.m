function [z, p] = findOptimalqGeometric(drone, SC, z0, ref, zbar)
    zk = zbar;
    zk_prev = z0;
    
    p0 = extractP(drone);
    pk = p0;
    pk_prev = zeros(1, drone.numParams);
    epsilon = 1e-1;
    stable = true;
    while sum(abs(zk - zk_prev), "all") > 1e-2 || ~stable
        [~, stateProgression, drone, stable] = trajectory(drone, SC, zk, false, ref);
        if stable
            disp("stable")
            [~, DG] = calcGDG(drone, SC, stateProgression, ref);
            q0 = [z0; pk];
            qk = [zk; pk];
            qhat = findProjection(DG, q0, qk);
            qdir = findqdir(qhat, qk);
            alpha = findAlpha(p0, pk, qdir, epsilon);
            zk_prev = zk;
            zk = zk + alpha * qdir(1:23);
            
            pk_prev = pk;
            pk = pk + alpha * qdir(24:41);
            drone = assignP(drone, pk);
            u = 1/2;
        else
            disp("unstable")
            zk = zk_prev + u * (zk - zk_prev);
            pk = pk_prev + u * (pk - pk_prev);
            drone = assignP(drone, pk);
            u = u/2;
        end  
            
    end

    z = zk;
    p = pk;
end

function [alpha] = findAlpha(p0, pk, qdir, epsilon)
    d = p0 - pk;
    pdir = qdir(24:41);
    alpha = d'*pdir + sqrt((d'*pdir)^2 - d'*d + epsilon^2);
end

function [qdir] = findqdir(qhat, qk)
    qdir = -(qhat - qk) / norm(qhat - qk);
end

function [qhat] = findProjection(DG, q0, qk)
    qhat =  qk + dot(DG', q0)/norm(DG)^2 .* DG';
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