classdef trajectorySensitivities
    methods(Static)
        function trajSens = calculateTrajectorySensitivities(state, SC)
            n = SC.numStates;
            trajSens = zeros(n, n, length(state));
            trajSens(:, :, 1) = SC.initialdxdxo;
            I = eye(n);
            for i = 2:length(state)
                dfdx_plus = SC.df(SC, state(i, :));
                dfdx = SC.df(SC, state(i-1, :));
        
                dxdxo = trajSens(:, :, i-1);
                trajSens(:, :, i) =  (I-(SC.T/2)*dfdx_plus)\((I + (SC.T/2)*dfdx) * dxdxo);
            end
        end

        function trajSens = calculateTrajectorySensitivitiesParameters(state, SC)
            n = SC.numStates;
            m = SC.numParams;
            trajSens = zeros(n, m, length(state));
            I = eye(n);
            trajSens(:, :, 1) = 0;
            for i = 2:length(state)
                dfdx = SC.df(SC, state(i-1, :));
                dfdx_plus = SC.df(SC, state(i, :));
                
                dfdp = SC.dfdp(SC, state(i-1, :));
                dfdp_plus = SC.dfdp(SC, state(i, :));
        
                dxdp = trajSens(:, :, i-1);
                trajSens(:, :, i) =  (I-(SC.T/2)*dfdx_plus)\((SC.T/2)*(dfdx*dxdp + dfdp + dfdp_plus) + dxdp);
            end
        end
        
        function d2TrajSens = calculateSecondOrderTrajectorySensitivities(state, trajSens, SC)
            n = SC.numStates;
            d2TrajSens = zeros(n, n, n, length(state));
            I = eye(n);
        
            for i = 2:length(state)
                dfdx_plus = SC.df(SC, state(i, :));
        
                d2xdx0 = d2TrajSens(:, :, :, i-1);
        
                d2fdx2_plus =  SC.d2fdx2(SC, state(i, :));
                dxdx0_plus = trajSens(:,:,i);
        
                d2fdx2 =  SC.d2fdx2(SC, state(i-1, :));
                dxdx0 = trajSens(:,:,i-1);
        
                dfdx =  SC.df(SC, state(i-1, :));
                
                a = inv(I - (SC.T/2)*dfdx_plus);
                b = conts({d2fdx2_plus, dxdx0_plus, dxdx0_plus}, {[3 1],[2 1]}, []);
                c = conts({d2fdx2, dxdx0, dxdx0}, {[3 1],[2 1]}, []);
                d = conts({d2xdx0, dfdx}, {[1 2]}, []); % how to know which dimension
                
                e = d2xdx0 + (SC.T/2)*(b+c+d);

                d2TrajSens(:, :, :, i) = conts({a, e}, {[2 1]}, []);
            end
        end

        function d2TrajSens = calculateSecondOrderTrajectorySensitivitiesParams(state, dxdxo_list, dxdp_list, SC)
            n = SC.numStates;
            m = SC.numParams;
            d2TrajSens = zeros(n, n, m, length(state));
            I = eye(n);

            for i = 2:length(state)
                d2fdx2 = SC.d2fdx2(SC, state(i-1, :));
                dxdp = dxdp_list(:, :, i-1);
                dxdx0 = dxdxo_list(:, :, i-1);
                d2fdxdp = SC.d2fdxdp(SC, state(i-1, :));
                dfdx = SC.df(SC, state(i-1, :));
                d2xdxodp = d2TrajSens(:, :, :, i-1);


                d2fdx2_plus = SC.d2fdx2(SC, state(i, :));
                dxdp_plus = dxdp_list(:, :, i);
                dxdx0_plus = dxdxo_list(:, :, i);
                d2fdxdp_plus = SC.d2fdxdp(SC, state(i, :));
                dfdx_plus = SC.df(SC, state(i, :));
                
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
        
        function trajSensTest = calculateTrajectorySensitivitiesTest(initialState, SC, numSteps)
            delta = 1e-5;
            sCtrl = stonkCtrllers;
            trajSensTest = zeros(SC.numStates, SC.numStates, numSteps);

            for index = 1:SC.numStates
                state = initialState;
                state(index) = initialState(index) + delta;
                [~, state_index_plus] = pendulum_trajectory_v2(SC, state, SC.controller, @sCtrl.dummyD, false);
                state(index) = initialState(index) - delta;
                [~, state_index_minus]  = pendulum_trajectory_v2(SC, state, SC.controller, @sCtrl.dummyD, false);
    
                for i = 1:numSteps
                    for j = 1:SC.numStates
                        dxdstate = (state_index_plus(i,j)-state_index_minus(i,j))/(2*delta);
                        trajSensTest(j, index, i) = dxdstate;
                    end
                end
            end            
        end

        function trajSensTestParams = calculateTrajectorySensitivitiesTestParams(state, SC, numSteps)
            delta = 1e-5;
            oldSC = SC;
            trajSensTestParams = zeros(SC.numStates, SC.numParams, numSteps);
            for index = 1:SC.numParams
                SC = SC.nudgeParams(SC, delta, index);
                [~, state_index_plus] = pendulum_trajectory_v2(SC, state, SC.controller, @stonkCtrllers.dummyD, false);
                SC = oldSC;
                SC = SC.nudgeParams(SC, -1*delta, index);
                [~, state_index_minus] = pendulum_trajectory_v2(SC, state, SC.controller, @stonkCtrllers.dummyD, false);
                SC = oldSC;

                for i = 1:numSteps
                    for j = 1:SC.numStates
                        dxdstate = (state_index_plus(i, j)-state_index_minus(i, j))/(2*delta);
                        trajSensTestParams(j, index, i) = dxdstate;
                    end
                end
            end
        end
        
        function secondTrajSensTest = calculateSecondOrderTrajectorySensitivitiesTest(initialState, SC, numSteps)
            delta = 1e-5;
            tjs = trajectorySensitivities;
            sCtrl = stonkCtrllers;

            secondTrajSensTest = zeros(SC.numStates, SC.numStates, SC.numStates, numSteps);

            fns = fieldnames(initialState);

            for index = 1:SC.numParams
                state = initialState;
                state.(fns{index}) = initialState.(fns{index}) + delta;
                [~, state_index_plus] = pendulum_trajectory_v2(SC, state, SC.controller, @sCtrl.dummyD, false);
                trajSens_plus = tjs.calculateTrajectorySensitivities(state_index_plus, SC);
            
                state.(fns{index}) = initialState.(fns{index}) - delta;
                [~, state_index_minus]  = pendulum_trajectory_v2(SC, state, SC.controller, @sCtrl.dummyD, false);
                trajSens_minus = tjs.calculateTrajectorySensitivities(state_index_minus, SC);
                
                trajSens1 = (trajSens_plus - trajSens_minus)/(2*delta);
                secondTrajSensTest(:, :, index, :) = trajSens1;
            end 
        end

        function secondTrajSensTest = calculateSecondOrderTrajectorySensitivitiesTestParams(state, SC, numSteps)
            delta = 1e-5;
            secondTrajSensTest = zeros(SC.numStates, SC.numStates, SC.numParams, numSteps);
            oldSC = SC;
            tjs = trajectorySensitivities;
            for index = 1:SC.numParams
                SC = SC.nudgeParams(SC, delta, index);
                [~, state_index_plus] = pendulum_trajectory_v2(SC, state, SC.controller, @stonkCtrllers.dummyD, false);
                trajSens_plus = tjs.calculateTrajectorySensitivities(state_index_plus, SC);
                SC = oldSC;

                SC = SC.nudgeParams(SC, -1*delta, index);
                [~, state_index_minus] = pendulum_trajectory_v2(SC, state, SC.controller, @stonkCtrllers.dummyD, false);
                trajSens_minus = tjs.calculateTrajectorySensitivities(state_index_minus, SC);
                SC = oldSC;

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

        function dG = calculate_dG(G, trajSens, d2xdxodp, tp, SC)
            sum = zeros(1,SC.numParams);
            for k = 1:SC.numStates
                for l = 1:SC.numStates
                    for i = 1:SC.numParams
                        sum(i) = sum(i) + sign(trajSens(k,l,tp)) * d2xdxodp(k,l,i,tp);
                    end
                end
            end
            dG = -1 * G^2 * sum;
        end

        function dGTest = calculate_dGTest(state, SC)
            delta = 1e-5;
            oldSC = SC;
            tjs = trajectorySensitivities;
            dGTest = zeros(1, SC.numParams);
            for index = 1:SC.numParams
                SC = SC.nudgeParams(SC, delta, index);
                [~, state_index_plus] = pendulum_trajectory_v2(SC, state, SC.controller, @stonkCtrllers.dummyD, false);
                trajSens_plus = tjs.calculateTrajectorySensitivities(state_index_plus, SC);
                Gplus = normalStuff.infimumOfReciprocalOfElementOneNorm(trajSens_plus);
                SC = oldSC;

                SC = SC.nudgeParams(SC, -1*delta, index);
                [~, state_index_minus] = pendulum_trajectory_v2(SC, state, SC.controller, @stonkCtrllers.dummyD, false);
                trajSens_minus = tjs.calculateTrajectorySensitivities(state_index_minus, SC);
                Gminus = normalStuff.infimumOfReciprocalOfElementOneNorm(trajSens_minus);
                SC = oldSC;

                dGTest(index) = (Gplus-Gminus)/(2*delta);
            end
        end
    end
end