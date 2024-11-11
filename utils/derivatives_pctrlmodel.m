classdef derivatives_pctrlmodel
    methods (Static)
        function dfdx = df(SC, state)
            x1 = state.x1;
            dfdx = [0 1; -(SC.g/SC.l)*cos(x1)-SC.pid_p/(SC.m*SC.l) -SC.D/(SC.m*SC.l)];
        end


        function SC = assignDerivatives(SC)
            SC.df = @derivatives_pctrlmodel.df;
            SC.numStates = 2;
            SC.initialdxdxo = eye(SC.numStates);
            SC.numParams = 2;
            SC.controller = @proposedCtrllers.ParameterizedPIDCtrl;
        end
    end
end