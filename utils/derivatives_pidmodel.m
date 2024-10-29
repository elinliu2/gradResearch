classdef derivatives_pidmodel
    methods (Static)
        function dfdx = df(SC, state)
            dfdx = [0 1 0 0 0; 
                -(SC.g/SC.l)*cos(state.x1) -SC.D/(SC.m*SC.l) SC.pid_i/(SC.m*SC.l) SC.pid_p/(SC.m*SC.l) SC.pid_d/(SC.m*SC.l);
                0 0 0 1 0;
                0 0 0 0 1;
                0 -1 0 0 0];
        end

        function dfdp = dfdp(SC, state)
            dfdp = zeros(3, 3);
            dfdp(2, :) = [state.ei/(SC.m*SC.l) state.e/(SC.m*SC.l) state.ed/(SC.m*SC.l)];
        end
        
        function d2fdx2 = d2fdx2(SC, x1)
            d2fdx2 = zeros(5,5,2);
            d2fdx2(2, 1, 1) = SC.g/SC.l * sin(x1);
        end

        function d2fdxdp = d2fdxdp(SC, x1)
            d2fdxdp = zeros(5,5,2);
        end

        function SC = assignDerivatives(SC)
            SC.df = @derivatives_pidmodel.df;
            SC.dfdp = @derivatives_pidmodel.dfdp;
            SC.d2fdx2 = @derivatives_pidmodel.d2fdx2;
            SC.d2fdxdp = @derivatives_pidmodel.d2fdxdp;
        end
    end
end