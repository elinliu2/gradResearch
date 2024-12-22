classdef DroneControllers
   methods (Static)
        function [u, SC] = dummyU(~, ~, ~, ~, ~, SC)
            u = SC.u;
        end
        function d = dummyD(disturbance, time)
            if (time < 5)
                d = disturbance;
            else 
                d = 0;
            end
        end
        function [u, model] = ParameterizedPIDCtrl(x1, x2, model)
            u = model.kp * -1*x1 + model.kd * -1*x2;
        end
   end
end

