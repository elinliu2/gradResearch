classdef trajectorySensitivitiesTest
    methods(Static)
        function trajSensTest = calculateTrajectorySensitivitiesParams(model, SC, initialState, numSteps, ref)
            delta = 1e-5;
            oldmodel = model;
            trajSensTest = zeros(model.numStates, model.numParams, numSteps);
            for index = 1:model.numStates+model.numParams
                [model, state] = model.nudgeParams(oldmodel, initialState, delta, index);
                [~, state_index_plus] = model.modelTrajectory(model, SC, state, false, ref);
                [model, state] = model.nudgeParams(oldmodel, initialState, -1*delta, index);
                [~, state_index_minus] = model.modelTrajectory(model, SC, state, false, ref);
                trajSensTest(:, index, :) = transpose(state_index_plus-state_index_minus)/(2*delta);
            end
        end

        function secondTrajSensTest = calculateSecondOrderTrajectorySensitivitiesParams(model, SC, initialState, numSteps, ref)
            delta = 1e-6;
            secondTrajSensTest = zeros(model.numStates, model.numStates, model.numStates + model.numParams, numSteps);
            oldmodel = model;
            n = model.numStates;
            for index = 1:model.numStates + model.numParams
                [model, state] = model.nudgeParams(oldmodel, initialState, delta, index);
                [~, state_index_plus] = model.modelTrajectory(model, SC, state, false, ref);
                trajSens_plus = model.calculateTrajectorySensitivitiesParams(model, SC, state_index_plus, ref);
               
                [model, state] = model.nudgeParams(oldmodel, initialState, -1*delta, index);
                [~, state_index_minus] = model.modelTrajectory(model, SC, state, false, ref);
                trajSens_minus = model.calculateTrajectorySensitivitiesParams(model, SC, state_index_minus, ref);

                dtrajSens = (trajSens_plus(:, 1:n, :) - trajSens_minus(:, 1:n, :))/(2*delta);

                secondTrajSensTest(:, :, index, :) = dtrajSens;
            end
        end

        function dGTest = calculate_dG(model, SC, initialState, ref)
            delta = 1e-5;
            oldmodel = model;
            dGTest = zeros(1, model.numParams);
            for index = 1:model.numStates + model.numParams
                [model, state] = model.nudgeParams(oldmodel, initialState, delta, index);
                [~, state_index_plus] = model.modelTrajectory(model, SC, state, false, ref);
                trajSens_plus = model.calculateTrajectorySensitivitiesParams(model, SC, state_index_plus, ref);
                Gplus = norms.infimumOfReciprocalOfOneNorm(trajSens_plus);

                [model, state] = model.nudgeParams(oldmodel, initialState, -1*delta, index);
                [~, state_index_minus] = model.modelTrajectory(model, SC, state, false, ref);
                trajSens_minus = model.calculateTrajectorySensitivitiesParams(model, SC, state_index_minus, ref);
                Gminus = norms.infimumOfReciprocalOfOneNorm(trajSens_minus);


                dGTest(index) = (Gplus-Gminus)/(2*delta);
            end
        end
    end
end