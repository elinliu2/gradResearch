classdef derivatives_x1x2ulmodel
    methods (Static)
        function dfdx = df(SC, state)
            x1 = state.x1;
            dfdx = [0 1; -(SC.g/SC.l)*cos(x1) -SC.D/(SC.m*SC.l)];
        end

        function dfdp = dfdp(SC, state)
            x1 = state.x1;
            x2 = state.x2;
            dfdp = [0 0; 1/(SC.m*SC.l) ((SC.g*sin(x1) + SC.D*x2/SC.m - SC.u/SC.m)/(SC.l^2)) ];
        end
        
        function d2fdx2 = d2fdx2(SC, state)
            d2fdx2 = zeros(2,2,2);
            d2fdx2(2, 1, 1) = SC.g/SC.l*sin(state.x1);
        end

        function d2fdxdp = d2fdxdp(SC, state)
            x1 = state.x1;
            d2fdxdp = zeros(2,2,2);
            d2fdxdp(2, 1, 2) = SC.g/(SC.l^2)*cos(x1);
            d2fdxdp(2, 2, 2) = SC.D/(SC.m*SC.l^2);

        end

        function SC = nudgeParams(SC, delta, index)
            if index == 1
                SC.u = SC.u + delta;
            elseif index == 2
                SC.l = SC.l + delta;
            end
        end

        function SC = assignDerivatives(SC)
            SC.df = @derivatives_x1x2ulmodel.df;
            SC.dfdp = @derivatives_x1x2ulmodel.dfdp;
            SC.d2fdx2 = @derivatives_x1x2ulmodel.d2fdx2;
            SC.d2fdxdp = @derivatives_x1x2ulmodel.d2fdxdp;
            SC.controller = @stonkCtrllers.dummyU;
            SC.nudgeParams = @derivatives_x1x2ulmodel.nudgeParams;

            SC.numStates = 2;
            SC.numParams = 2;
            SC.initialdxdxo = eye(SC.numStates);
        end
    end
end