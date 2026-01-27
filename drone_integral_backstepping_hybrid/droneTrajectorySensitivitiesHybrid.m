classdef droneTrajectorySensitivitiesHybrid
    methods(Static)
        function trajSens = calculateTrajectorySensitivitiesParams(model, stateProgression, ref, time, index)
            numX = model.numX;
            n = model.numStates;
            m = model.numParams;
            trajSens = zeros(n, n+m, index);
            trajSens(1:n, 1:n, 1) = eye(n);
            for i = 2:index
                T = time(i) - time(i-1);
                dfdx = derivativesDroneIBHybrid.dfdx(model, stateProgression(i-1, :), ref);
                dfdx_plus = derivativesDroneIBHybrid.dfdx(model, stateProgression(i, :), ref);
                dfdy = derivativesDroneIBHybrid.dfdy(model, stateProgression(i-1, :), ref);
                
                dfdq = derivativesDroneIBHybrid.dfdq(model, stateProgression(i-1, :), ref);
                dfdq_plus = derivativesDroneIBHybrid.dfdq(model, stateProgression(i, :), ref);

                dxdq = trajSens(1:numX, :, i-1);
                dydq = trajSens(numX+1:n, :, i-1);
                dy12dq = dydq(1:2, :);
                
                dxdq_plus = (eye(numX) - T/2*dfdx_plus)\(dxdq + T/2*(dfdq_plus + dfdx*dxdq + dfdq + 2*dfdy*dydq));

                dgdx = derivativesDroneIBHybrid.dgdx(model, stateProgression(i-1, :), ref);
                dgdx_plus = derivativesDroneIBHybrid.dgdx(model, stateProgression(i, :), ref);

                dgdq = derivativesDroneIBHybrid.dgdq(model, stateProgression(i-1, :), ref);
                dgdq_plus = derivativesDroneIBHybrid.dgdq(model, stateProgression(i, :), ref);
                
                dy12dq_plus = 1/T*(dgdx_plus*dxdq_plus + dgdq_plus - dgdx*dxdq - dgdq);
                dy34dq_plus = 1/T*(dy12dq_plus-dy12dq);
                
                trajSens(1:numX, :, i) =  dxdq_plus;
                trajSens(numX+1:numX+2, :, i) =  dy12dq_plus;
                trajSens(numX+3:n, :, i) =  dy34dq_plus;
            end
        end
    
        function secondTrajSens = calculateSecondOrderTrajectorySensitivitiesParams(model, stateProgression, trajSens, ref, time, index)
            numX = model.numX;
            n = model.numStates;
            m = model.numParams;
            secondTrajSens =  zeros(n, n, n+m, index);
            for i = 2:index
                T = time(i) - time(i-1);
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

                l1 = d2xdz0dq + T/2*(a1+b1+c1 + d1+e1+f1+g1 + 2*(h1+i1+j1+k1));
                m1 = inv(eye(numX) - T/2*dfdx_plus);

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
                
                d2y12dz0dq_plus = 1/T*(a2 + b2 + c2 - d2 - e2 - f2);
                d2y34dz0dq_plus = 1/T*(d2y12dz0dq_plus - d2y12dz0dq);

                secondTrajSens(1:numX, :, :, i) = d2xdz0dq_plus;
                secondTrajSens(numX+1: numX+2, :, :, i) = d2y12dz0dq_plus;
                secondTrajSens(numX+3: n, :, :, i) = d2y34dz0dq_plus; 
            end
        end

        function nextTrajSens = calcTrajSensSingle(model, SC, state, nextState, trajSensSingle, ref)
            numX = model.numX;
            n = model.numStates;
            dfdx = derivativesDroneIBHybrid.dfdx(model, state, ref);
            dfdx_plus = derivativesDroneIBHybrid.dfdx(model, nextState, ref);
            dfdy = derivativesDroneIBHybrid.dfdy(model, state, ref);
            
            dfdq = derivativesDroneIBHybrid.dfdq(model, state, ref);
            dfdq_plus = derivativesDroneIBHybrid.dfdq(model, nextState, ref);

            dxdq = trajSensSingle(1:numX, :);
            dydq = trajSensSingle(numX+1:n, :);
            dy12dq = dydq(1:2, :);
            
            dxdq_plus = (eye(numX) - SC.T/2*dfdx_plus)\(dxdq + SC.T/2*(dfdq_plus + dfdx*dxdq + dfdq + 2*dfdy*dydq));

            dgdx = derivativesDroneIBHybrid.dgdx(model, state, ref);
            dgdx_plus = derivativesDroneIBHybrid.dgdx(model, nextState, ref);

            dgdq = derivativesDroneIBHybrid.dgdq(model, state, ref);
            dgdq_plus = derivativesDroneIBHybrid.dgdq(model, nextState, ref);
            
            dy12dq_plus = 1/SC.T*(dgdx_plus*dxdq_plus + dgdq_plus - dgdx*dxdq - dgdq);
            dy34dq_plus = 1/SC.T*(dy12dq_plus-dy12dq);
            
            nextTrajSens = [dxdq_plus; dy12dq_plus; dy34dq_plus];
        end

        function nextSecondTrajSens = calculateSecondOrderTrajectorySensitivitiesParamsSingle(model, SC, state, nextState, trajSens, nextTrajSens, secondTrajSens, ref)
            numX = model.numX;
            n = model.numStates;
            m = model.numParams;
            nextSecondTrajSens =  zeros(n, n, n+m);

            d2fdx2 = derivativesDroneIBHybrid.d2fdx2(model, state, ref);
            d2fdx2_plus = derivativesDroneIBHybrid.d2fdx2(model, nextState, ref);
            
            d2fdxdy = derivativesDroneIBHybrid.d2fdxdy(model, state, ref);
            d2fdxdy_plus = derivativesDroneIBHybrid.d2fdxdy(model, nextState, ref);
            
            d2fdxdq = derivativesDroneIBHybrid.d2fdxdq(model, state, ref);
            d2fdxdq_plus = derivativesDroneIBHybrid.d2fdxdq(model, nextState, ref);

            dfdx = derivativesDroneIBHybrid.dfdx(model, state, ref);
            dfdx_plus = derivativesDroneIBHybrid.dfdx(model, nextState, ref);

            d2fdydx = derivativesDroneIBHybrid.d2fdydx(model, state, ref);
            d2fdy2 = derivativesDroneIBHybrid.d2fdy2(model, state, ref);
            d2fdydq = derivativesDroneIBHybrid.d2fdydq(model, state, ref);
            dfdy = derivativesDroneIBHybrid.dfdy(model, state, ref);

            dxdq = trajSens(1:numX, :);
            dxdq_plus = nextTrajSens(1:numX, :);
            dxdz0 = dxdq(:, 1:n, :);
            dxdz0_plus = dxdq_plus(:, 1:n, :);

            dydq = trajSens(numX+1:n, :);
            dydq_plus = nextTrajSens(numX+1:n, :);
            dydz0 = dydq(:, 1:n, :);

            d2xdz0dq = secondTrajSens(1:numX, :, :);
            d2ydz0dq = secondTrajSens(numX+1:n, :, :);
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

            d2gdx2 = derivativesDroneIBHybrid.d2gdx2(model, state, ref);
            d2gdx2_plus = derivativesDroneIBHybrid.d2gdx2(model, nextState, ref);
            d2gdxdq = derivativesDroneIBHybrid.d2gdxdq(model, state, ref);
            d2gdxdq_plus = derivativesDroneIBHybrid.d2gdxdq(model, nextState, ref);

            dgdx = derivativesDroneIBHybrid.dgdx(model, state, ref);
            dgdx_plus = derivativesDroneIBHybrid.dgdx(model, nextState, ref);
            
            a2 = conts({d2gdx2_plus, dxdq_plus, dxdz0_plus}, {[3 1],[2 1]}, []);
            b2 = conts({d2gdxdq_plus, dxdz0_plus}, {[2 1]}, []);
            c2 = conts({dgdx_plus, d2xdz0dq_plus}, {[2 1]}, []);

            d2 = conts({d2gdx2, dxdq, dxdz0}, {[3 1],[2 1]}, []);
            e2 = conts({d2gdxdq, dxdz0}, {[2 1]}, []);
            f2 = conts({dgdx, d2xdz0dq}, {[2 1]}, []);
            
            d2y12dz0dq_plus = 1/SC.T*(a2 + b2 + c2 - d2 - e2 - f2);
            d2y34dz0dq_plus = 1/SC.T*(d2y12dz0dq_plus - d2y12dz0dq);

            nextSecondTrajSens(1:numX, :, :) = d2xdz0dq_plus;
            nextSecondTrajSens(numX+1: numX+2, :, :) = d2y12dz0dq_plus;
            nextSecondTrajSens(numX+3: n, :, :) = d2y34dz0dq_plus; 
        end

        function G = calcG(model, SC, stateProgression, ref)
            n = model.numStates;
            m = model.numParams;
            trajSens = zeros(n, n+m, size(stateProgression, 1));
            trajSens(1:n, 1:n, 1) = eye(n);
            
            G = 1/sum(abs(trajSens), 'all');
            for i = 2:size(stateProgression, 1)
                nextTrajSens = calcTrajSensSingle(model, SC, nextState, state, trajSens, ref);
                G_k = 1/sum(abs(nextTrajSens), 'all');
                if G_k < G
                    G = G_k;
                end
                trajSens = nextTrajSens;
            end
        end

        function dG = calculate_dG(model, SC, stateProgression, ref)
            n = model.numStates;
            m = model.numParams;
            trajSens = zeros(n, n+m);
            trajSens(1:n, 1:n, 1) = eye(n);
            secondTrajSens = zeros(n, n, n+m);
            G = 1/sum(abs(trajSens), 'all');

            minTrajSens = trajSens;
            minSecondTrajSensParam = secondTrajSens;
            for i = 2:size(stateProgression, 1)
                state = stateProgression(i-1, :);
                nextState = stateProgression(i, :);
                nextTrajSens = droneTrajectorySensitivitiesHybrid.calcTrajSensSingle(model, SC, state, nextState, trajSens, ref);
                nextSecondTrajSensParam = droneTrajectorySensitivitiesHybrid.calculateSecondOrderTrajectorySensitivitiesParamsSingle(model, SC, state, nextState, trajSens, nextTrajSens, secondTrajSens, ref);
                G_k = 1/sum(abs(nextTrajSens), 'all');
                if G_k < G
                    G = G_k;
                    minTrajSens = nextTrajSens;
                    minSecondTrajSensParam = nextSecondTrajSensParam;
                end
                trajSens = nextTrajSens;
                secondTrajSens = nextSecondTrajSensParam;
            end 

            dGsum = zeros(1,model.numStates + model.numParams);

            for j = 1:model.numStates + model.numParams
                for k = 1:model.numStates
                    dGsum(j) = dGsum(j) + transpose(sign(minTrajSens(:, k))) * minSecondTrajSensParam(:,k,j);
                end
            end

            dG = -1 .* G^2 .* dGsum;
        end
    end
end