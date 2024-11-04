classdef derivatives_pidmodel
    methods (Static)
        function dfdx = df(SC, state)
            dfdx = [0 1 0 0 0; 
                -(SC.g/SC.l)*cos(state.x1) -SC.D/(SC.m*SC.l) SC.pid_i/(SC.m*SC.l) SC.pid_p/(SC.m*SC.l) SC.pid_d/(SC.m*SC.l);
                0 0 0 1 0;
                0 0 0 0 1;
                (SC.g/SC.l)*cos(state.x1) SC.D/(SC.m*SC.l) -SC.pid_i/(SC.m*SC.l) -SC.pid_p/(SC.m*SC.l) -SC.pid_d/(SC.m*SC.l)];
        end

        function dfdp = dfdp(SC, state)
            dfdp = zeros(5, 3);
            dfdp(2, :) = [state.ei/(SC.m*SC.l) state.e/(SC.m*SC.l) state.ed/(SC.m*SC.l)];
        end
        
        function d2fdx2 = d2fdx2(SC, state)
            d2fdx2 = zeros(5,5,5);
            d2fdx2(2, 1, 1) = SC.g/SC.l * sin(state.x1);
        end

        function d2fdxdp = d2fdxdp(SC, ~)
            d2fdxdp = zeros(5,5,3);
            d2fdxdp(2, 3, 1) = 1/(SC.m*SC.l);
            d2fdxdp(2, 3, 2) = 1/(SC.m*SC.l);
            d2fdxdp(2, 3, 3) = 1/(SC.m*SC.l);
        end

        function SC = nudgeParams(SC, delta, index)
            if index == 1
                SC.pid_i = SC.pid_i + delta;
            elseif index == 2
                SC.pid_p = SC.pid_p + delta;
            elseif index == 3
                SC.pid_d = SC.pid_d + delta;
            end
        end

        function SC = assignDerivatives(SC)
            SC.df = @derivatives_pidmodel.df;
            SC.dfdp = @derivatives_pidmodel.dfdp;
            SC.d2fdx2 = @derivatives_pidmodel.d2fdx2;
            SC.d2fdxdp = @derivatives_pidmodel.d2fdxdp;
            SC.numStates = 5;
           
            SC.numParams = 3;
            SC.controller = @proposedCtrllers.ParameterizedPIDCtrl;
            SC.nudgeParams = @derivatives_pidmodel.nudgeParams;

            SC.initialdxdxo = [ 1 0 0 0 0 ; 
                                0 1 0 0 0 ;
                                0 0 0 0 0 ;
                                -1 0 0 0 0;
                                0 0 0 0 0];
        end
    end
end