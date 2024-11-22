classdef stonkCtrllers
    methods (Static)
        function [u, SC] = dummyU(~, ~, ~, ~, ~, SC)
            u = SC.u;
        end
        
        function d = dummyD(SC, time)
            if (time > 0.05)
                d = 0;
            else 
                d = SC.disturbance;
            end 
        end
    end
end