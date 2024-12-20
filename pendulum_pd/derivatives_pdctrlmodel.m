classdef derivatives_pdctrlmodel
    methods (Static)
        function dfdx = df(model, state)
            x1 = state(1);
            dfdx = [0 1; -(model.g/model.l)*cos(x1)-model.kp/(model.m*model.l) -(model.D+model.kd)/(model.m*model.l)];
        end
        function dfdp = dfdp(model, state)
            x1 = state(1);
            x2 = state(2);
            dfdp = [0 0; -x1/(model.m*model.l) -x2/(model.m*model.l)];
        end
        function d2fdx2 = d2fdx2(model, state)
            x1 = state(1);
            d2fdx2 = zeros(2,2,2);
            d2fdx2(2, 1, 1) = model.g/model.l*sin(x1);
        end
        function d2fdxdp = d2fdxdp(model, ~)
            d2fdxdp = zeros(2,2,2);
            d2fdxdp(2,1,1) = -1/(model.m*model.l);
            d2fdxdp(2,2,2) = -1/(model.m*model.l);
        end
        function model = nudgeParams(model, delta, index)
            if index == 1
                model.kp = model.kp + delta;
            elseif index == 2
                model.kd = model.kd + delta;
            end
        end
    end
end