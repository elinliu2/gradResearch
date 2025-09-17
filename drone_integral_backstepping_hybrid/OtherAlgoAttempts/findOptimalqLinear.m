function [z,p] = findOptimalqLinear(drone, SC, z0, ref, zbar)
    zk = zbar;
    zk_prev = 0;
    
    n = drone.numStates;
    m = drone.numParams;
    p0 = extractP(drone);
    pk = p0;
    epsilon = 1e-1;
    
    while sum(abs(zk - zk_prev), "all") > 1e-2
        [~, stateProgression, drone, stable] = trajectory(drone, SC, zk, false, ref);
        if stable
           [~, DG] = calcGDG(drone, SC, stateProgression, ref);
           zk_prev = zk;
           [zk, pk, ~, ~] = NewtonRaphson(zk, pk, DG, n, m, z0, p0, epsilon);
 
           assignP(drone, pk)
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

function [zk_1, pk_1, l1_1, l2_1] = NewtonRaphson(zk, pk, DG, n, m, z0, p0, epsilon)
    zk_0 = zk;
    pk_0 = pk + 0.01*rand(m, 1);
    l1_0 = 1/2;
    l2_0 = 1/2;    

    zk_1 = zeros(n, 1);
    pk_1 = zeros(m, 1);
    l1_1 = 0;
    l2_1 = 0;

    while norm([zk_1; pk_1; l1_1; l2_1] - [zk_0; pk_0; l1_0; l2_0]) > 1e-3
        
        Gradient = [
            zk-z0+(l1_0*DG*[eye(n); zeros(m, n)])';
            (l1_0*DG*[zeros(n, m); eye(m)])'+l2_0*2*(pk_0 - p0);
            DG*([zk_0; pk_0] - [zk; pk])];
        
        Hessian = [
            [eye(n) zeros(n, m) (DG*[eye(n); zeros(m, n)])' zeros(n, 1)];
            [zeros(m, n) l2_0*eye(m)  (DG*[zeros(n, m); eye(m)])' pk_0];
            [DG*[eye(n); zeros(m, n)] DG*[zeros(n, m); eye(m)] 0 0]];

        v = Hessian\Gradient;
        w = null(Hessian);

        p = pk_0-p0;
        vp = v(24:41);
        wp = w(24:41, 1);

        t = 1/(2*p'*wp)*(-2*p'*vp+p'*p-epsilon^2);
        x = v + t*w(:, 1);
        next = [zk_0; pk_0; l1_0; l2_0] - x;
        
        zk_0 = zk_1;
        pk_0 = pk_1;
        l1_0 = l1_1;
        l2_0 = l2_1;  

        zk_1 = next(1:23);
        pk_1 = next(24:41);
        l1_1 = next(42);
        l2_1 = next(43);
    end
end