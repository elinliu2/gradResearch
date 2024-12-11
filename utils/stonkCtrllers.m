classdef stonkCtrllers
    methods (Static)
        function [u, SC] = dummyU(~, ~, ~, ~, ~, SC)
            u = SC.u;
        end
        
        function d = dummyD(SC, time)
            d = SC.disturbance;
        end
    end
end