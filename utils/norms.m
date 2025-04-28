classdef norms
    methods (Static)

        function min = infimumOfReciprocalOfOneNorm(trajSens)
            trajSens = abs(trajSens);
            min = 1/sum(trajSens(:, :, 1), 'all');

            for i = 2:size(trajSens, 3)
                if (1/sum(trajSens(:, :, i), 'all') < min)
                    min = 1/sum(trajSens(:, :, i), 'all');
                end
            end
        end
        
        function max = supremum(trajSens)
            max = sum(trajSens(:, :, 1), 'all');
            for i = 2:size(trajSens, 3)
                if (sum(trajSens(:, :, i), 'all') > max)
                    max = sum(trajSens(:, :, i), 'all');
                end
            end
        end
        
        function max = secondOrderMax(secondTrajSens)
            max = sum(secondTrajSens(:, :, :, 1), 'all');
            for i = 2:size(secondTrajSens, 4)
                if (sum(secondTrajSens(:, :, :, i), 'all') > max)
                    max = sum(secondTrajSens(:, :, :, i), 'all');
                end
            end
        end
    end
end