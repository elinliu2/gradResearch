classdef droneTrajectorySensitivitiesHybrid
    methods(Static)
        function trajSens = calculateTrajectorySensitivities(model, SC, stateProgression, ref)
            n = model.numStates;
            trajSens = zeros(n, n, size(stateProgression, 1));
            trajSens(:, :, 1) = model.initialdxdxo;
            I = eye(n);
            for i = 2:size(stateProgression,1)
                dfdx_plus = model.dfdx(model, stateProgression(i, :), ref);
                dfdx = model.dfdx(model, stateProgression(i-1, :), ref);

                dthetarefdx_plus = derivativesDroneIBHybrid.dthetarefdx(model, stateProgression(i, :), ref);
                dthetarefdx = derivativesDroneIBHybrid.dthetarefdx(model, stateProgression(i-1, :), ref);

                dphirefdx_plus = derivativesDroneIBHybrid.dphirefdx(model, stateProgression(i, :), ref);
                dphirefdx = derivativesDroneIBHybrid.dphirefdx(model, stateProgression(i-1, :), ref);

                df22dx = [zeros(1, 19) 1 0 0 0];
                df23dx = [zeros(1, 20) 1 0 0];
        
                dFdx = [SC.T/2*(dfdx); -1/SC.T*[dthetarefdx; dphirefdx; df22dx; df23dx] ];
                dFdx_plus = [SC.T/2*(dfdx_plus); 1/SC.T*[dthetarefdx_plus; dphirefdx_plus; df22dx; df23dx] ];

                dxdxo = trajSens(:, :, i-1);
                a = [dxdxo(1:19, :); zeros(4, model.numStates)];

                trajSens(:, :, i) =  (I-dFdx_plus)\(a + dFdx*dxdxo);
            end
        end

        function trajSens = calculateTrajectorySensitivitiesParameters(model, SC, stateProgression, ref)
            n = model.numStates;
            m = model.numParams;
            trajSens = zeros(n, m, size(stateProgression, 1));
            I = eye(n);
            trajSens(:, :, 1) = 0;
            for i = 2:size(stateProgression, 1)
                dfdx = model.dfdx(model, stateProgression(i-1, :), ref);
                dfdx_plus = model.dfdx(model, stateProgression(i, :), ref);
                
                dfdp = model.dfdp(model, stateProgression(i-1, :), ref);
                dfdp_plus = model.dfdp(model, stateProgression(i, :), ref);

                dthetarefdx_plus = derivativesDroneIBHybrid.dthetarefdx(model, stateProgression(i, :), ref);
                dthetarefdx = derivativesDroneIBHybrid.dthetarefdx(model, stateProgression(i-1, :), ref);

                dphirefdx_plus = derivativesDroneIBHybrid.dphirefdx(model, stateProgression(i, :), ref);
                dphirefdx = derivativesDroneIBHybrid.dphirefdx(model, stateProgression(i-1, :), ref);

                dthetarefdp_plus = derivativesDroneIBHybrid.dthetarefdp(model, stateProgression(i, :), ref);
                dthetarefdp = derivativesDroneIBHybrid.dthetarefdp(model, stateProgression(i-1, :), ref);

                dphirefdp_plus = derivativesDroneIBHybrid.dphirefdp(model, stateProgression(i, :), ref);
                dphirefdp = derivativesDroneIBHybrid.dphirefdp(model, stateProgression(i-1, :), ref);

                dFdp = [SC.T/2*dfdx; -1/SC.T*dthetarefdx; -1/SC.T*dphirefdx; zeros(2, model.numStates)];
                dFdp_plus = [SC.T/2*dfdx_plus; 1/SC.T*dthetarefdx_plus; 1/SC.T*dphirefdx_plus; zeros(2, model.numStates)];
                
                dxdp = trajSens(:, :, i-1);
                a = [SC.T/2 * (dfdp + dfdp_plus) + dxdp(1:19, :); 
                    1/SC.T*(dthetarefdp_plus-dthetarefdp); 
                    1/SC.T*(dphirefdp_plus-dphirefdp);
                    zeros(2, model.numParams)];

                trajSens(:, :, i) =  (I-dFdp_plus)\(a + dFdp*dxdp);
            end
        end
        
        function d2TrajSens = calculateSecondOrderTrajectorySensitivitiesParams(model, SC, stateProgression, dxdxo_list, dxdp_list, ref)
            n = model.numStates;
            m = model.numParams;
            d2TrajSens = zeros(n, n, m, size(stateProgression, 1));
            I = eye(n);

            for i = 2:size(stateProgression, 1)
                dxdp = dxdp_list(:, :, i-1);
                dxdx0 = dxdxo_list(:, :, i-1);
                dxdp_plus = dxdp_list(:, :, i);
                dxdx0_plus = dxdxo_list(:, :, i);

                d2fdx2 = model.d2fdx2(model, stateProgression(i-1, :), ref);
                d2fdxdp = model.d2fdxdp(model, stateProgression(i-1, :), ref);
                dfdx = model.dfdx(model, stateProgression(i-1, :), ref);

                d2fdx2_plus = model.d2fdx2(model, stateProgression(i, :), ref);
                d2fdxdp_plus = model.d2fdxdp(model, stateProgression(i, :), ref);
                dfdx_plus = model.dfdx(model, stateProgression(i, :), ref);
                
                b1 = conts({d2fdx2, dxdp, dxdx0}, {[3 1],[2 1]}, []);
                c1 = conts({d2fdxdp, dxdx0}, {[2 1]}, []);
                d1 = conts({d2fdx2_plus, dxdp_plus, dxdx0_plus}, {[3 1],[2 1]}, []);
                e1 = conts({d2fdxdp_plus, dxdx0_plus}, {[2 1]}, []);
                
                d2thetarefdx2 = derivativesDroneIBHybrid.d2thetarefdx2(model, stateProgression(i-1, :), ref);
                d2thetarefdxdp = derivativesDroneIBHybrid.d2thetarefdxdp(model, stateProgression(i-1, :), ref);
                dthetarefdx = derivativesDroneIBHybrid.dthetarefdx(model, stateProgression(i-1, :), ref);
             
                d2thetarefdx2_plus = derivativesDroneIBHybrid.d2thetarefdx2(model, stateProgression(i, :), ref);
                d2thetarefdxdp_plus = derivativesDroneIBHybrid.d2thetarefdxdp(model, stateProgression(i, :), ref);
                dthetarefdx_plus = derivativesDroneIBHybrid.dthetarefdx(model, stateProgression(i, :), ref);

                b2 = conts({d2thetarefdx2, dxdp, dxdx0}, {[2 1],[1 1]}, []);
                c2 = conts({d2thetarefdxdp, dxdx0}, {[1 1]}, []);
                d2 = conts({d2thetarefdx2_plus, dxdp_plus, dxdx0_plus}, {[2 1],[1 1]}, []);
                e2 = conts({d2thetarefdxdp_plus, dxdx0_plus}, {[1 1]}, []);

                d2phirefdx2 = derivativesDroneIBHybrid.d2phirefdx2(model, stateProgression(i-1, :), ref);
                d2phirefdxdp = derivativesDroneIBHybrid.d2phirefdxdp(model, stateProgression(i-1, :), ref);
                dphirefdx = derivativesDroneIBHybrid.dphirefdx(model, stateProgression(i-1, :), ref);

                d2phirefdx2_plus = derivativesDroneIBHybrid.d2phirefdx2(model, stateProgression(i, :), ref);
                d2phirefdxdp_plus = derivativesDroneIBHybrid.d2phirefdxdp(model, stateProgression(i, :), ref);
                dphirefdx_plus = derivativesDroneIBHybrid.dphirefdx(model, stateProgression(i, :), ref);

                b3 = conts({d2phirefdx2, dxdp, dxdx0}, {[2 1],[1 1]}, []);
                c3 = conts({d2phirefdxdp, dxdx0}, {[1 1]}, []);
                d3 = conts({d2phirefdx2_plus, dxdp_plus, dxdx0_plus}, {[2 1],[1 1]}, []);
                e3 = conts({d2phirefdxdp_plus, dxdx0_plus}, {[1 1]}, []);

                d2xdxodp = d2TrajSens(:, :, :, i-1);
                
                a = zeros(model.numStates, model.numStates, model.numParams);
                a(1:19, :, :) = SC.T/2*(b1 + c1 + d1 + e1);
                a(20, :, :) = 1/SC.T*(d2 + e2 - b2 - c2);
                a(21, :, :) = 1/SC.T*(d3 + e3 - b3 - c3);
                a(22:23, :, :) = zeros(2, model.numStates, model.numParams);

                f = [SC.T/2*dfdx; -1/SC.T*dthetarefdx; -1/SC.T*dphirefdx; zeros(2, model.numStates) ];
                g = [SC.T/2*dfdx_plus; 1/SC.T*dthetarefdx_plus; 1/SC.T*dphirefdx_plus; zeros(2, model.numStates) ];

                h = conts({f, d2xdxodp}, {[2 1]}, []);
                d2TrajSens(:, :, :, i) = conts({inv(I - g), (a+h)}, {[2 1]}, []);
            end
        end
    end
end