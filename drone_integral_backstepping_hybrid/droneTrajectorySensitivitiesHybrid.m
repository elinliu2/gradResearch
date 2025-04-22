classdef droneTrajectorySensitivitiesHybrid
    methods(Static)
        function trajSens = calculateTrajectorySensitivities(model, SC, stateProgression, ref)
            n = model.numStates;
            trajSens = zeros(n, n, size(stateProgression, 1));
            trajSens(:, :, 1) = model.initialdxdxo;
            for i = 2:size(stateProgression,1)
                dfdx_plus = derivativesDroneIBHybrid.dfdx(model, stateProgression(i, :), ref);
                dfdx = derivativesDroneIBHybrid.dfdx(model, stateProgression(i-1, :), ref);
                dfdy = derivativesDroneIBHybrid.dfdy(model, stateProgression(i-1, :), ref);

                dthetarefdx_plus = derivativesDroneIBHybrid.dthetarefdx(model, stateProgression(i, :), ref);
                dthetarefdx = derivativesDroneIBHybrid.dthetarefdx(model, stateProgression(i-1, :), ref);

                dphirefdx_plus = derivativesDroneIBHybrid.dphirefdx(model, stateProgression(i, :), ref);
                dphirefdx = derivativesDroneIBHybrid.dphirefdx(model, stateProgression(i-1, :), ref);

                dgdx_plus = [dthetarefdx_plus; dphirefdx_plus];
                dgdx = [dthetarefdx; dphirefdx];
        
                dxdzo = trajSens(1:model.numX, :, i-1);
                dydzo = trajSens(model.numX+1:model.numX+model.numY, :, i-1);

                dxdzo_plus = (eye(model.numX) - SC.T/2*dfdx_plus)\(dxdzo + SC.T/2*(dfdx*dxdzo + 2*dfdy*dydzo));
                dy12dzo_plus = 1/SC.T*(dgdx_plus*dxdzo_plus - dgdx*dxdzo);
                dy34dzo_plus = 1/SC.T*(dy12dzo_plus - dydzo(1:2, :));

                trajSens(:, :, i) =  [dxdzo_plus; dy12dzo_plus; dy34dzo_plus];
            end
        end

        function trajSensParams = calculateTrajectorySensitivitiesParameters(model, SC, stateProgression, trajSens, ref)
            n = model.numStates;
            m = model.numParams;
            trajSensParams = zeros(n, n+m, size(stateProgression, 1));
            trajSensParams(:, :, 1) = 0;
            trajSensParams(1:n, 1:n, 1) = eye(n);
            for i = 2:size(stateProgression, 1)
                dfdx = derivativesDroneIBHybrid.dfdx(model, stateProgression(i-1, :), ref);
                dfdx_plus = derivativesDroneIBHybrid.dfdx(model, stateProgression(i, :), ref);
                dfdy = derivativesDroneIBHybrid.dfdy(model, stateProgression(i-1, :), ref);
                
                dfdp = derivativesDroneIBHybrid.dfdp(model, stateProgression(i-1, :), ref);
                dfdp_plus = derivativesDroneIBHybrid.dfdp(model, stateProgression(i, :), ref);

                dthetarefdx_plus = derivativesDroneIBHybrid.dthetarefdx(model, stateProgression(i, :), ref);
                dthetarefdx = derivativesDroneIBHybrid.dthetarefdx(model, stateProgression(i-1, :), ref);

                dphirefdx_plus = derivativesDroneIBHybrid.dphirefdx(model, stateProgression(i, :), ref);
                dphirefdx = derivativesDroneIBHybrid.dphirefdx(model, stateProgression(i-1, :), ref);

                dgdx_plus = [dthetarefdx_plus; dphirefdx_plus];
                dgdx = [dthetarefdx; dphirefdx];

                dthetarefdp_plus = derivativesDroneIBHybrid.dthetarefdp(model, stateProgression(i, :), ref);
                dthetarefdp = derivativesDroneIBHybrid.dthetarefdp(model, stateProgression(i-1, :), ref);

                dphirefdp_plus = derivativesDroneIBHybrid.dphirefdp(model, stateProgression(i, :), ref);
                dphirefdp = derivativesDroneIBHybrid.dphirefdp(model, stateProgression(i-1, :), ref);

                dgdp_plus = [dthetarefdp_plus; dphirefdp_plus];
                dgdp = [dthetarefdp; dphirefdp];

                dxdp = trajSensParams(1:model.numX, n+1:n+m, i-1);
                dydp = trajSensParams(model.numX+1:model.numX+model.numY, n+1:n+m, i-1);

                dxdp_plus = (eye(model.numX) - SC.T/2*dfdx_plus)\(dxdp + SC.T/2*(dfdx*dxdp + 2*dfdy*dydp + dfdp_plus + dfdp));
                dy12dp_plus = 1/SC.T*(dgdx_plus*dxdp_plus + dgdp_plus - dgdx*dxdp - dgdp);
                dy34dp_plus = 1/SC.T*(dy12dp_plus - dydp(1:2, :));

                trajSensParams(:, :, i) =  [ trajSens(:, :, i) [dxdp_plus; dy12dp_plus; dy34dp_plus]];
            end
        end
        
        function d2TrajSens = calculateSecondOrderTrajectorySensitivitiesParams(model, SC, stateProgression, dzdzo_list, dzdq_list, secondTrajSensTestParam, ref)
            numX = model.numX;
            numY = model.numY;
            n = model.numStates;
            m = model.numParams;
            d2TrajSens = zeros(n, n, n+m, size(stateProgression, 1));

            for i = 2:size(stateProgression, 1)
                dxdzo =      dzdzo_list(1:numX, :,       i-1);
                dxdp =       dzdq_list( 1:numX, n+1:n+m, i-1);
                dxdzo_plus = dzdzo_list(1:numX, :,       i);
                dxdp_plus =  dzdq_list( 1:numX, n+1:n+m, i);
                dydzo =      dzdzo_list(numX+1:n, :,       i-1);
                
                d2xdzo2  =  d2TrajSens(1:numX,   1:n, 1:n,     i-1);
                d2ydzo2  =  d2TrajSens(numX+1:n, 1:n, 1:n,     i-1);
                d2y12dzo2 = d2ydzo2(1:2, :, :);
                d2xdzodp =  d2TrajSens(1:numX,   1:n, n+1:n+m, i-1);
                d2ydzodp =  d2TrajSens(numX+1:n, 1:n, n+1:n+m, i-1);

                d2fdx2 = derivativesDroneIBHybrid.d2fdx2(model, stateProgression(i-1, :), ref);
                dfdx = derivativesDroneIBHybrid.dfdx(model, stateProgression(i-1, :), ref);
                
                d2fdy2 = derivativesDroneIBHybrid.d2fdy2(model, stateProgression(i-1, :), ref);
                dfdy = derivativesDroneIBHybrid.dfdy(model, stateProgression(i-1, :), ref);

                d2fdx2_plus = derivativesDroneIBHybrid.d2fdx2(model, stateProgression(i, :), ref);
                dfdx_plus = derivativesDroneIBHybrid.dfdx(model, stateProgression(i, :), ref);
               
                a1 = conts({d2fdx2_plus, dxdzo_plus, dxdzo_plus}, {[2 1], [3 1]}, []);
                b1 = conts({d2fdx2, dxdzo, dxdzo}, {[2 1], [3 1]}, []);
                c1 = conts({dfdx, d2xdzo2}, {[2 1]}, []);
                d1 = conts({d2fdy2, dydzo, dydzo}, {[2 1], [3 1]}, []);
                e1 = conts({dfdy, d2ydzo2}, {[2 1]}, []);
                f1 = d2xdzo2 + SC.T/2*(a1 + b1 + c1 + 2 * (d1 + e1));
                d2xdzo2_plus = conts({inv(eye(numX) - dfdx_plus), f1}, {[2 1]}, []);

                dthetarefdx_plus = derivativesDroneIBHybrid.dthetarefdx(model, stateProgression(i, :), ref);
                dthetarefdx = derivativesDroneIBHybrid.dthetarefdx(model, stateProgression(i-1, :), ref);

                dphirefdx_plus = derivativesDroneIBHybrid.dphirefdx(model, stateProgression(i, :), ref);
                dphirefdx = derivativesDroneIBHybrid.dphirefdx(model, stateProgression(i-1, :), ref);

                dgdx_plus = [dthetarefdx_plus; dphirefdx_plus];
                dgdx = [dthetarefdx; dphirefdx];

                d2thetarefdx2_plus = derivativesDroneIBHybrid.d2thetarefdx2(model, stateProgression(i, :), ref);
                d2thetarefdx2 = derivativesDroneIBHybrid.d2thetarefdx2(model, stateProgression(i-1, :), ref);

                d2phirefdx2_plus = derivativesDroneIBHybrid.d2phirefdx2(model, stateProgression(i, :), ref);
                d2phirefdx2 = derivativesDroneIBHybrid.d2phirefdx2(model, stateProgression(i-1, :), ref);

                d2gdx2_plus = zeros(2, numX, numX);
                d2gdx2_plus(1, :, :) = d2thetarefdx2_plus; 
                d2gdx2_plus(2, :, :) = d2phirefdx2_plus;
                d2gdx2 = zeros(2, numX, numX);
                d2gdx2(1, :, :) = d2thetarefdx2; 
                d2gdx2(2, :, :) = d2phirefdx2;

                % d2xdzo2_plus_test = secondTrajSensTestParam(1:19, 1:23, 1:23, i);
                % d2xdzo2_test = secondTrajSensTestParam(1:19, 1:23, 1:23, i-1);
                
                a2 = conts({d2gdx2_plus, dxdzo_plus, dxdzo_plus}, {[2 1], [3 1]}, []);
                b2 = conts({dgdx_plus, d2xdzo2_plus}, {[2 1]}, []);
                c2 = conts({d2gdx2, dxdzo, dxdzo}, {[2 1], [3 1]}, []);
                d2 = conts({dgdx, d2xdzo2}, {[2 1]}, []);
                d2y12dzo2_plus = 1/SC.T*(a2 + b2 - c2 - d2);
                d2y34dzo2_plus = 1/SC.T*(d2y12dzo2_plus - d2y12dzo2);

                d2ydzo2_plus = [d2y12dzo2_plus; d2y34dzo2_plus];

                d2xdzodp_plus = zeros(numX, n, m);
                d2ydzodp_plus = zeros(numY, n, m);

                d2TrajSens(1:numX,   1:n, 1:n,     i)   = d2xdzo2_plus;
                d2TrajSens(1:numX,   1:n, n+1:n+m, i) = d2xdzodp_plus;
                d2TrajSens(numX+1:n, 1:n, 1:n,     i) = d2ydzo2_plus;
                d2TrajSens(numX+1:n, 1:n, n+1:n+m, i) = d2ydzodp_plus;
            end
        end
    end
end