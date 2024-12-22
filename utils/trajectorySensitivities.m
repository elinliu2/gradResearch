classdef trajectorySensitivities
    methods(Static)
        function trajSens = calculateTrajectorySensitivities(model, SC, stateProgression)
            n = model.numStates;
            trajSens = zeros(n, n, length(stateProgression));
            trajSens(:, :, 1) = model.initialdxdxo;
            I = eye(n);
            for i = 2:length(stateProgression)
                dfdx_plus = model.dfdx(model, stateProgression(i, :));
                dfdx = model.dfdx(model, stateProgression(i-1, :));
        
                dxdxo = trajSens(:, :, i-1);
                trajSens(:, :, i) =  (I-(SC.T/2)*dfdx_plus)\((I + (SC.T/2)*dfdx) * dxdxo);
            end
        end

        function trajSens = calculateTrajectorySensitivitiesParameters(model, SC, stateProgression)
            n = model.numStates;
            m = model.numParams;
            trajSens = zeros(n, m, length(stateProgression));
            I = eye(n);
            trajSens(:, :, 1) = 0;
            for i = 2:length(stateProgression)
                dfdx = model.dfdx(model, stateProgression(i-1, :));
                dfdx_plus = model.dfdx(model, stateProgression(i, :));
                
                dfdp = model.dfdp(model, stateProgression(i-1, :));
                dfdp_plus = model.dfdp(model, stateProgression(i, :));
        
                dxdp = trajSens(:, :, i-1);
                trajSens(:, :, i) =  (I-(SC.T/2)*dfdx_plus)\((SC.T/2)*(dfdx*dxdp + dfdp + dfdp_plus) + dxdp);
            end
        end
        
        function d2TrajSens = calculateSecondOrderTrajectorySensitivitiesParams(model, SC, stateProgression, dxdxo_list, dxdp_list)
            n = model.numStates;
            m = model.numParams;
            d2TrajSens = zeros(n, n, m, length(stateProgression));
            I = eye(n);

            for i = 2:length(stateProgression)
                d2fdx2 = model.d2fdx2(model, stateProgression(i-1, :));
                dxdp = dxdp_list(:, :, i-1);
                dxdx0 = dxdxo_list(:, :, i-1);
                d2fdxdp = model.d2fdxdp(model, stateProgression(i-1, :));
                dfdx = model.dfdx(model, stateProgression(i-1, :));
                d2xdxodp = d2TrajSens(:, :, :, i-1);


                d2fdx2_plus = model.d2fdx2(model, stateProgression(i, :));
                dxdp_plus = dxdp_list(:, :, i);
                dxdx0_plus = dxdxo_list(:, :, i);
                d2fdxdp_plus = model.d2fdxdp(model, stateProgression(i, :));
                dfdx_plus = model.dfdx(model, stateProgression(i, :));
                
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
        
        function trajSensTest = calculateTrajectorySensitivitiesTest(model, SC, initialState, numSteps)
            delta = 1e-5;
            trajSensTest = zeros(model.numStates, model.numStates, numSteps);

            for index = 1:model.numStates
                state = initialState;
                state(index) = initialState(index) + delta;
                [~, state_index_plus] = trajectory(model, SC, state, false);
                state(index) = initialState(index) - delta;
                [~, state_index_minus]  = trajectory(model, SC, state, false);
    
                for i = 1:numSteps
                    for j = 1:model.numStates
                        dxdstate = (state_index_plus(i,j)-state_index_minus(i,j))/(2*delta);
                        trajSensTest(j, index, i) = dxdstate;
                    end
                end
            end            
        end

        function trajSensTestParams = calculateTrajectorySensitivitiesTestParams(model, SC, initialState, numSteps)
            delta = 1e-5;
            oldmodel = model;
            trajSensTestParams = zeros(model.numStates, model.numParams, numSteps);
            for index = 1:model.numParams
                model = model.nudgeParams(model, delta, index);
                [~, state_index_plus] = trajectory(model, SC, initialState, false);
                model = oldmodel;
                model = model.nudgeParams(model, -1*delta, index);
                [~, state_index_minus] = trajectory(model, SC, initialState, false);
                model = oldmodel;

                for i = 1:numSteps
                    for j = 1:model.numStates
                        dxdstate = (state_index_plus(i, j)-state_index_minus(i, j))/(2*delta);
                        trajSensTestParams(j, index, i) = dxdstate;
                    end
                end
            end
        end
        
        function secondTrajSensTest = calculateSecondOrderTrajectorySensitivitiesTest(model, SC, initialState, numSteps)
            delta = 1e-5;
            tjs = trajectorySensitivities;

            secondTrajSensTest = zeros(model.numStates, model.numStates, model.numStates, numSteps);

            fns = fieldnames(initialState);

            for index = 1:model.numParams
                state = initialState;
                state.(fns{index}) = initialState.(fns{index}) + delta;
                [~, state_index_plus] = trajectory(model, SC, state, false);
                trajSens_plus = tjs.calculateTrajectorySensitivities(state_index_plus, model);
            
                state.(fns{index}) = initialState.(fns{index}) - delta;
                [~, state_index_minus]  = trajectory(model, SC, state, false);
                trajSens_minus = tjs.calculateTrajectorySensitivities(state_index_minus, model);
                
                trajSens1 = (trajSens_plus - trajSens_minus)/(2*delta);
                secondTrajSensTest(:, :, index, :) = trajSens1;
            end 
        end

        function secondTrajSensTest = calculateSecondOrderTrajectorySensitivitiesTestParams(model, SC, initialState, numSteps)
            delta = 1e-5;
            secondTrajSensTest = zeros(model.numStates, model.numStates, model.numParams, numSteps);
            oldmodel = model;
            tjs = trajectorySensitivities;
            for index = 1:model.numParams
                model = model.nudgeParams(model, delta, index);
                [~, state_index_plus] = trajectory(model, SC, initialState, false);
                trajSens_plus = tjs.calculateTrajectorySensitivities(model, SC, state_index_plus);
                model = oldmodel;

                model = model.nudgeParams(model, -1*delta, index);
                [~, state_index_minus] = trajectory(model, SC, initialState, false);
                trajSens_minus = tjs.calculateTrajectorySensitivities(model, SC, state_index_minus);
                model = oldmodel;

                dtrajSens = (trajSens_plus - trajSens_minus)/(2*delta);

                secondTrajSensTest(:, :, index, :) = dtrajSens;
            end
        end
        
        function tp = findtp(trajSens)
            min = 1/norm(trajSens(:,:,1), 1);
            tp = 1;
            for i = 2:length(trajSens)
                if (1/norm(trajSens(:,:,i), 1) < min)
                    min = 1./norm(trajSens(:,:,i), 1);
                    tp = i;
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