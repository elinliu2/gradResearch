classdef derivatives_pdctrlmodel
    methods (Static)
        function dfdx = df(SC, state)
            x1 = state(1);
            dfdx = [0 1; -(SC.g/SC.l)*cos(x1)-SC.pid_p/(SC.m*SC.l) -(SC.D+SC.pid_d)/(SC.m*SC.l)];
        end
        function dfdp = dfdp(SC, state)
            x1 = state(1);
            x2 = state(2);
            dfdp = [0 0; -x1/(SC.m*SC.l) -x2/(SC.m*SC.l)];
        end
        function d2fdx2 = d2fdx2(SC, state)
            x1 = state(1);
            d2fdx2 = zeros(2,2,2);
            d2fdx2(2, 1, 1) = SC.g/SC.l*sin(x1);
        end
        function d2fdxdp = d2fdxdp(SC, ~)
            d2fdxdp = zeros(2,2,2);
            d2fdxdp(2,1,1) = -1/(SC.m*SC.l);
            d2fdxdp(2,2,2) = -1/(SC.m*SC.l);
        end
        function SC = nudgeParams(SC, delta, index)
            if index == 1
                SC.pid_p = SC.pid_p + delta;
            elseif index == 2
                SC.pid_d = SC.pid_d + delta;
            end
        end
        function SC = assignDerivatives(SC)
            SC.df = @derivatives_pdctrlmodel.df;
            SC.dfdp = @derivatives_pdctrlmodel.dfdp;
            SC.d2fdx2 = @derivatives_pdctrlmodel.d2fdx2;
            SC.d2fdxdp = @derivatives_pdctrlmodel.d2fdxdp;
            SC.controller = @proposedCtrllers.ParameterizedPIDCtrl;
            SC.nudgeParams = @derivatives_pdctrlmodel.nudgeParams;

            SC.numStates = 2;
            SC.numParams = 2;
            SC.initialdxdxo = eye(2);
        end
    end
end