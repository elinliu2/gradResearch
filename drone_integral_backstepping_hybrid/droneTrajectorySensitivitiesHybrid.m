classdef droneTrajectorySensitivitiesHybrid
    methods(Static)
        function trajSens = calculateTrajectorySensitivitiesParams(model, SC, stateProgression, ref)
            numX = model.numX;
            n = model.numStates;
            m = model.numParams;
            trajSens = zeros(n, n+m, size(stateProgression, 1));
            trajSens(1:n, 1:n, 1) = eye(n);
            for i = 2:size(stateProgression, 1)
                dfdx = derivativesDroneIBHybrid.dfdx(model, stateProgression(i-1, :), ref);
                dfdx_plus = derivativesDroneIBHybrid.dfdx(model, stateProgression(i, :), ref);
                dfdy = derivativesDroneIBHybrid.dfdy(model, stateProgression(i-1, :), ref);
                
                dfdq = derivativesDroneIBHybrid.dfdq(model, stateProgression(i-1, :), ref);
                dfdq_plus = derivativesDroneIBHybrid.dfdq(model, stateProgression(i, :), ref);

                dxdq = trajSens(1:numX, :, i-1);
                dydq = trajSens(numX+1:n, :, i-1);
                dy12dq = dydq(1:2, :);
                
                dxdq_plus = (eye(numX) - SC.T/2*dfdx_plus)\(dxdq + SC.T/2*(dfdq_plus + dfdx*dxdq + dfdq + 2*dfdy*dydq));

                dgdx = derivativesDroneIBHybrid.dgdx(model, stateProgression(i-1, :), ref);
                dgdx_plus = derivativesDroneIBHybrid.dgdx(model, stateProgression(i, :), ref);

                dgdq = derivativesDroneIBHybrid.dgdq(model, stateProgression(i-1, :), ref);
                dgdq_plus = derivativesDroneIBHybrid.dgdq(model, stateProgression(i, :), ref);
                
                dy12dq_plus = 1/SC.T*(dgdx_plus*dxdq_plus + dgdq_plus - dgdx*dxdq - dgdq);
                dy34dq_plus = 1/SC.T*(dy12dq_plus-dy12dq);
                
                trajSens(1:numX, :, i) =  dxdq_plus;
                trajSens(numX+1:numX+2, :, i) =  dy12dq_plus;
                trajSens(numX+3:n, :, i) =  dy34dq_plus;

            end
        end
    
        function secondTrajSens = calculateSecondOrderTrajectorySensitivitiesParams(model, SC, stateProgression, trajSens, ref)
            numX = model.numX;
            n = model.numStates;
            m = model.numParams;
            secondTrajSens =  zeros(n, n, n+m, size(stateProgression, 1));
            for i = 2:size(stateProgression, 1)
                d2fdx2 = derivativesDroneIBHybrid.d2fdx2(model, stateProgression(i-1, :), ref);
                d2fdx2_plus = derivativesDroneIBHybrid.d2fdx2(model, stateProgression(i, :), ref);
                
                d2fdxdy = derivativesDroneIBHybrid.d2fdxdy(model, stateProgression(i-1, :), ref);
                d2fdxdy_plus = derivativesDroneIBHybrid.d2fdxdy(model, stateProgression(i, :), ref);
                
                d2fdxdq = derivativesDroneIBHybrid.d2fdxdq(model, stateProgression(i-1, :), ref);
                d2fdxdq_plus = derivativesDroneIBHybrid.d2fdxdq(model, stateProgression(i, :), ref);

                dfdx = derivativesDroneIBHybrid.dfdx(model, stateProgression(i-1, :), ref);
                dfdx_plus = derivativesDroneIBHybrid.dfdx(model, stateProgression(i, :), ref);

                d2fdydx = derivativesDroneIBHybrid.d2fdydx(model, stateProgression(i-1, :), ref);
                d2fdy2 = derivativesDroneIBHybrid.d2fdy2(model, stateProgression(i-1, :), ref);
                d2fdydq = derivativesDroneIBHybrid.d2fdydq(model, stateProgression(i-1, :), ref);
                dfdy = derivativesDroneIBHybrid.dfdy(model, stateProgression(i-1, :), ref);

                dxdq = trajSens(1:numX, :, i-1);
                dxdq_plus = trajSens(1:numX, :, i);
                dxdz0 = dxdq(:, 1:n, :);
                dxdz0_plus = dxdq_plus(:, 1:n, :);

                dydq = trajSens(numX+1:n, :, i-1);
                dydq_plus = trajSens(numX+1:n, :, i);
                dydz0 = dydq(:, 1:n, :);

                d2xdz0dq = secondTrajSens(1:numX, :, :, i-1);
                d2ydz0dq = secondTrajSens(numX+1:n, :, :, i-1);
                d2y12dz0dq = d2ydz0dq(1:2, :, :, :);

                a1 = conts({d2fdx2_plus, dxdq_plus, dxdz0_plus}, {[3 1],[2 1]}, []);
                b1 = conts({d2fdxdy_plus, dydq_plus, dxdz0_plus}, {[3 1],[2 1]}, []);
                c1 = conts({d2fdxdq_plus, dxdz0_plus}, {[2 1]}, []);

                d1 = conts({d2fdx2, dxdq, dxdz0}, {[3 1],[2 1]}, []);
                e1 = conts({d2fdxdy, dydq, dxdz0}, {[3 1],[2 1]}, []);
                f1 = conts({d2fdxdq, dxdz0}, {[2 1]}, []);
                g1 = conts({dfdx, d2xdz0dq}, {[2 1]}, []);

                h1 = conts({d2fdydx, dxdq, dydz0}, {[3 1],[2 1]}, []);
                i1 = conts({d2fdy2, dydq, dydz0}, {[3 1],[2 1]}, []);
                j1 = conts({d2fdydq, dydz0}, {[2 1]}, []);
                k1 = conts({dfdy, d2ydz0dq}, {[2 1]}, []);

                l1 = d2xdz0dq + SC.T/2*(a1+b1+c1 + d1+e1+f1+g1 + 2*(h1+i1+j1+k1));
                m1 = inv(eye(numX) - SC.T/2*dfdx_plus);

                d2xdz0dq_plus = conts({m1, l1}, {[2 1]}, []);

                d2gdx2 = derivativesDroneIBHybrid.d2gdx2(model, stateProgression(i-1, :), ref);
                d2gdx2_plus = derivativesDroneIBHybrid.d2gdx2(model, stateProgression(i, :), ref);
                d2gdxdq = derivativesDroneIBHybrid.d2gdxdq(model, stateProgression(i-1, :), ref);
                d2gdxdq_plus = derivativesDroneIBHybrid.d2gdxdq(model, stateProgression(i, :), ref);

                dgdx = derivativesDroneIBHybrid.dgdx(model, stateProgression(i-1, :), ref);
                dgdx_plus = derivativesDroneIBHybrid.dgdx(model, stateProgression(i, :), ref);

                
                a2 = conts({d2gdx2_plus, dxdq_plus, dxdz0_plus}, {[3 1],[2 1]}, []);
                b2 = conts({d2gdxdq_plus, dxdz0_plus}, {[2 1]}, []);
                c2 = conts({dgdx_plus, d2xdz0dq_plus}, {[2 1]}, []);

                d2 = conts({d2gdx2, dxdq, dxdz0}, {[3 1],[2 1]}, []);
                e2 = conts({d2gdxdq, dxdz0}, {[2 1]}, []);
                f2 = conts({dgdx, d2xdz0dq}, {[2 1]}, []);
                
                d2y12dz0dq_plus = 1/SC.T*(a2 + b2 + c2 - d2 - e2 - f2);
                d2y34dz0dq_plus = 1/SC.T*(d2y12dz0dq_plus - d2y12dz0dq);

                secondTrajSens(1:numX, :, :, i) = d2xdz0dq_plus;
                secondTrajSens(numX+1: numX+2, :, :, i) = d2y12dz0dq_plus;
                secondTrajSens(numX+3: n, :, :, i) = d2y34dz0dq_plus;
                
            end
        end
    end
end