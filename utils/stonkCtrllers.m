classdef stonkCtrllers
    methods (Static)
        function u = dummyU(~, ~, ~, ~, ~, SC)
            u = SC.u;
        end
        
        function d = dummyD(SC)
            d = SC.disturbance;
        end
    end
end