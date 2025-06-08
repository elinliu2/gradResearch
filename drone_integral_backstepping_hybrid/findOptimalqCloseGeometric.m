function [z, p] = findOptimalqCloseGeometric(drone, SC, z0, ref, dist, zbar)
    zk = zbar;
    zk_prev = z0;
    
    p0 = extractP(drone);
    pk = p0;
    pk_prev = zeros(1, drone.numParams);
    epsilon = 1e-1;
    stable = true;
    while sum(abs(zk - zk_prev), "all") > 1e-2 || ~stable
        [~, stateProgression, drone, stable] = trajectory(drone, SC, zk, false, ref, dist);
        if stable
            disp("stable")
            [~, DG] = calcGDG(drone, SC, stateProgression, ref);
            q0 = [z0; pk];
            qk = [zk; pk];
            qhat = findProjection(DG, q0, qk);
            qdir = findqdir(qhat, qk);
            alpha = findAlpha(p0, pk, qdir, epsilon);
            zk_prev = zk;
            pk_prev = pk;
            pk = pk + alpha * qdir(24:41);

            drone = assignP(drone, pk);
            zk = findClosestzbar(drone, SC, z0, ref, dist);

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
    qdir = (qhat - qk) / norm(qhat - qk);
end

function [qhat] = findProjection(DG, q0, qk)
    qhat =  qk + dot(DG', q0)/norm(DG)^2 .* DG';
end