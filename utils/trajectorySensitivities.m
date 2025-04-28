classdef trajectorySensitivities
    methods(Static)
        function trajSens = calculateTrajectorySensitivities(model, SC, stateProgression, ref)
            n = model.numStates;
            trajSens = zeros(n, n, size(stateProgression, 1));
            trajSens(:, :, 1) = model.initialdxdxo;
            I = eye(n);
            for i = 2:size(stateProgression,1)
                dfdx_plus = model.dfdx(model, stateProgression(i, :), ref);
                dfdx = model.dfdx(model, stateProgression(i-1, :), ref);
        
                dxdxo = trajSens(:, :, i-1);
                trajSens(:, :, i) =  (I-(SC.T/2)*dfdx_plus)\((I + (SC.T/2)*dfdx) * dxdxo);
            end
        end

        function trajSensParams = calculateTrajectorySensitivitiesParameters(model, SC, stateProgression, ref)
            n = model.numStates;
            m = model.numParams;
            trajSensParams = zeros(n, m, size(stateProgression, 1));
            I = eye(n);
            trajSensParams(:, :, 1) = 0;
            for i = 2:size(stateProgression, 1)
                dfdx = model.dfdx(model, stateProgression(i-1, :), ref);
                dfdx_plus = model.dfdx(model, stateProgression(i, :), ref);
                
                dfdp = model.dfdp(model, stateProgression(i-1, :), (i-1)*SC.T, ref);
                dfdp_plus = model.dfdp(model, stateProgression(i, :), i*SC.T, ref);
        
                dxdp = trajSensParams(:, :, i-1);
                trajSensParams(:, :, i) =  (I-(SC.T/2)*dfdx_plus)\((SC.T/2)*(dfdx*dxdp + dfdp + dfdp_plus) + dxdp);
            end
        end
        
        function d2TrajSens = calculateSecondOrderTrajectorySensitivitiesParams(model, SC, stateProgression, dxdxo_list, dxdp_list, ref)
            n = model.numStates;
            m = model.numParams;
            d2TrajSens = zeros(n, n, m, size(stateProgression, 1));
            I = eye(n);

            for i = 2:size(stateProgression, 1)
                d2fdx2 = model.d2fdx2(model, stateProgression(i-1, :), ref);
                dxdp = dxdp_list(:, :, i-1);
                dxdx0 = dxdxo_list(:, :, i-1);
                d2fdxdp = model.d2fdxdp(model, stateProgression(i-1, :), ref);
                dfdx = model.dfdx(model, stateProgression(i-1, :), ref);
                d2xdxodp = d2TrajSens(:, :, :, i-1);

                d2fdx2_plus = model.d2fdx2(model, stateProgression(i, :), ref);
                dxdp_plus = dxdp_list(:, :, i);
                dxdx0_plus = dxdxo_list(:, :, i);
                d2fdxdp_plus = model.d2fdxdp(model, stateProgression(i, :), ref);
                dfdx_plus = model.dfdx(model, stateProgression(i, :), ref);
                
                a = inv(I - (SC.T/2)*dfdx_plus);
                b = conts({d2fdx2, dxdp, dxdx0}, {[3 1],[2 1]}, []);
                c = conts({d2fdxdp, dxdx0}, {[2 1]}, []);
                d = conts({dfdx, d2xdxodp}, {[2 1]}, []);

                e = conts({d2fdx2_plus, dxdp_plus, dxdx0_plus}, {[3 1],[2 1]}, []);
                f = conts({d2fdxdp_plus, dxdx0_plus}, {[2 1]}, []);
                g = (SC.T/2)*(b + c + d + e + f) + d2xdxodp;
                
                d2TrajSens(:, :, :, i) = conts({a, g}, {[2 1]}, []);
            end
        end

        function tp = findtp(trajSens)
            trajSens = abs(trajSens);
            min = 1/sum(trajSens(:, :, 1), 'all');
            tp = 1;
            for i = 2:size(trajSens,3)
                if (1/sum(trajSens(:, :, i), 'all') < min)
                    min = 1/sum(trajSens(:, :, i), 'all');
                    tp = i;
                end
            end
        end

        function dG = calculate_dG(G, trajSens, secondTrajSensParam, tp, model)
            sum = zeros(1,model.numStates + model.numParams);

            for j = 1:model.numStates + model.numParams
                for i = 1:model.numStates
                   sum(j) = sum(j) + transpose(sign(trajSens(:, i,tp))) * secondTrajSensParam(:,i,j,tp);
                end
            end

            dG = -1 .* G^2 .* sum;
        end

        
    end
end