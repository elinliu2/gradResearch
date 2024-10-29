classdef splatting 
    methods (Static)
        function trajSens = slapPack(dx1dx1, dx1dx2, dx2dx1, dx2dx2)
            trajSens = zeros(2, 2, length(dx1dx1));
            for i = 1:length(dx1dx1)
                trajSens(:, :, i) = [dx1dx1(i) dx1dx2(i); dx2dx1(i) dx2dx2(i) ];
            end
        end 
        
        function secondTrajSens = slapPile(trajSens1, trajSens2)
            secondTrajSens = zeros(2, 2, 2, length(trajSens1));
            for i = 1:length(trajSens1)
                temp = trajSens1(:, :, i);
                temp(:, :, 2) = trajSens2(:, :, i);
                secondTrajSens(:, :, :, i) = temp;
            end
        end

        function [x1, x2] = extractX1andX2FromState(state)
            x1 = zeros([1 length(state)]);
            x2 = zeros([1 length(state)]);
            for i = 1:length(state)
                x1(i) = state(i).x1;
                x2(i) = state(i).x2;
            end
        end
    end
end