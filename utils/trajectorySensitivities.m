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
            min = 1/norm(trajSens(:,:,1), 1);
            tp = 1;
            for i = 2:size(trajSens,3)
                if (1/norm(trajSens(:,:,i), 1) < min)
                    min = 1./norm(trajSens(:,:,i), 1);
                    tp = i;
                end
            end
        end

        function tp = findMaxTime(trajSens)
            infnorm = max(trajSens, [],"all");
            vecTrajSens = reshape(trajSens,[],1);
            tp = 0;
            itemsPerTimestamp = 1;
            for i = 1:ndims(trajSens)-1
                itemsPerTimestamp = itemsPerTimestamp * size(trajSens, i);
            end
            for i=1:length(vecTrajSens)
                if vecTrajSens(i) == infnorm
                    tp = idivide(length(vecTrajSens),int16(itemsPerTimestamp));
                    return
                end
            end
        end

        function dG = calculate_dG(G, trajSens, d2xdxodp, tp, model)
            sum = zeros(1,model.numParams);
            for k = 1:model.numStates
                for l = 1:model.numStates
                    for i = 1:model.numParams
                        sum(i) = sum(i) + sign(trajSens(k,l,tp)) * d2xdxodp(k,l,i,tp);
                    end
                end
            end
            dG = -1 * G^2 * sum;
        end

        function dGTest = calculate_dGTest(model, SC, initialState)
            delta = 1e-5;
            oldmodel = model;
            tjs = trajectorySensitivities;
            dGTest = zeros(1, model.numParams);
            for index = 1:model.numParams
                model = model.nudgeParams(model, delta, index);
                [~, state_index_plus] = trajectory(model, SC, initialState, false);
                trajSens_plus = tjs.calculateTrajectorySensitivities(model, SC, state_index_plus);
                Gplus = norms.infimumOfReciprocalOfElementOneNorm(trajSens_plus);
                model = oldmodel;

                model = model.nudgeParams(model, -1*delta, index);
                [~, state_index_minus] = trajectory(model, SC, initialState, false);
                trajSens_minus = tjs.calculateTrajectorySensitivities(model, SC, state_index_minus);
                Gminus = norms.infimumOfReciprocalOfElementOneNorm(trajSens_minus);
                model = oldmodel;

                dGTest(index) = (Gplus-Gminus)/(2*delta);
            end
        end
    end
end