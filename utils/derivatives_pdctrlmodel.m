classdef derivatives_pdctrlmodel
    methods (Static)
        function dfdx = df(SC, state)
            x1 = state.x1;
            dfdx = [0 1; -(SC.g/SC.l)*cos(x1)-SC.pid_p/(SC.m*SC.l) -(SC.D+SC.pid_d)/(SC.m*SC.l)];
        end


        function SC = assignDerivatives(SC)
            SC.df = @derivatives_pdctrlmodel.df;
            SC.numStates = 2;
            SC.initialdxdxo = eye(SC.numStates);
            SC.numParams = 2;
            SC.controller = @proposedCtrllers.ParameterizedPIDCtrl;
        end
    end
end