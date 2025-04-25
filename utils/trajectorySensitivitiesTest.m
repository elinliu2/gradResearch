classdef trajectorySensitivitiesTest
    methods(Static)
        function trajSensTest = calculateTrajectorySensitivitiesParams(model, SC, initialState, numSteps, ref)
            delta = 1e-5;
            oldmodel = model;
            trajSensTest = zeros(model.numStates, model.numParams, numSteps);
            for index = 1:model.numStates+model.numParams
                [model, state] = model.nudgeParams(oldmodel, initialState, delta, index);
                [~, state_index_plus] = trajectory(model, SC, state, false, ref);
                [model, state] = model.nudgeParams(oldmodel, initialState, -1*delta, index);
                [~, state_index_minus] = trajectory(model, SC, state, false, ref);
                trajSensTest(:, index, :) = transpose(state_index_plus-state_index_minus)/(2*delta);
            end
        end

        function secondTrajSensTest = calculateSecondOrderTrajectorySensitivitiesParams(model, SC, initialState, numSteps, ref)
            delta = 1e-5;
            secondTrajSensTest = zeros(model.numStates, model.numStates, model.numStates + model.numParams, numSteps);
            oldmodel = model;
            for index = 1:model.numStates + model.numParams
                [model, state] = model.nudgeParams(oldmodel, initialState, delta, index);
                [~, state_index_plus] = trajectory(model, SC, state, false, ref);
                trajSens_plus = droneTrajectorySensitivitiesHybrid.calculateTrajectorySensitivities(model, SC, state_index_plus, ref);
               
                [model, state] = model.nudgeParams(oldmodel, initialState, -1*delta, index);
                [~, state_index_minus] = trajectory(model, SC, state, false, ref);
                trajSens_minus = droneTrajectorySensitivitiesHybrid.calculateTrajectorySensitivities(model, SC, state_index_minus, ref);

                dtrajSens = (trajSens_plus - trajSens_minus)/(2*delta);

                secondTrajSensTest(:, :, index, :) = dtrajSens;
            end
        end

    end
end