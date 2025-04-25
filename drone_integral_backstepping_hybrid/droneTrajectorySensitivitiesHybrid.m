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
    end
end