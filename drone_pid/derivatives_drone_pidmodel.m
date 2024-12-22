classdef derivatives_drone_pidmodel
    methods (Static)
        function dfdx = dfdx(model, state)
            
        end
        function dfdp = dfdp(model, state)
            
        end
        function d2fdx2 = d2fdx2(model, state)
            
        end
        function d2fdxdp = d2fdxdp(model, state)
            
        end
        function model = nudgeParams(model, delta, index)
            if index == 1
                model.kp = model.kp + delta;
            elseif index == 2
                model.ki = model.ki + delta;
            elseif index == 3
                model.kd = model.kd + delta;
            end
        end
    end
end

