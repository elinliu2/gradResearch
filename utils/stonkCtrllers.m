classdef stonkCtrllers
    methods (Static)
        function [u, SC] = dummyU(~, ~, ~, ~, ~, SC)
            u = SC.u;
        end
        
        function d = dummyD(SC, time)
            if (time < 5)
                d = SC.disturbance;
            else 
                d = 0;
            end
        end
    end
end