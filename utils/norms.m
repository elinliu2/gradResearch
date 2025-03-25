classdef norms
    methods (Static)

        function min = infimumOfReciprocalOfOneNorm(trajSens)
            min = 1/norm(trajSens(:,:, 1), 1);

            for i = 2:length(trajSens)
                if (1/norm(trajSens(:, :,i), 1) < min)
                    min = 1/norm(trajSens(:,:,i), 1);
                end
            end
        end

        function min = infimumOfReciprocalOfElementOneNorm(trajSens)
            min = 1/norms.elementOneNorm(trajSens(:,:, 1));

            for i = 2:size(trajSens, 3)
                if (1/norms.elementOneNorm(trajSens(:, :,i)) < min)
                    min = 1/norms.elementOneNorm(trajSens(:,:,i));
                end
            end
        end

        function min = infimum(trajSens)
            min = norm(trajSens(:,:, 1), 1);

            for i = 2:length(trajSens)
                if (norm(trajSens(:, :,i), 1) < min)
                    min = norm(trajSens(:,:,i), 1);
                end
            end
        end
        
        function max = supremum(trajSens)
            max = norm(trajSens(:,:, 1), 1);
            for i = 2:size(trajSens, 3)
                if (norm(trajSens(:,:, i), 1) > max)
                    max = norm(trajSens(:,:,i), 1);
                end
            end
        end
        
        function max = secondOrderMax(secondTrajSens)
            max = norm(secondTrajSens(:, :, :, 1), "fro");
            for i = 1:size(secondTrajSens, 4)
                if (norm(secondTrajSens(:, :, :, i), "fro") > max)
                    max = norm(secondTrajSens(:, :, :, i), "fro");
                end
            end
        end

        function trajNorm = trajOneNorm(trajSens)
            trajNorm = 0;
            for i = 1:length(trajSens)
                trajNorm = trajNorm + norm(trajSens(:,:,i), 1);
            end
        end

        function sum = elementOneNorm(a)
            sum = 0;
            for i = 1:size(a, 1)
                for j = 1:size(a, 2)
                   sum = sum + abs(a(i, j));
                end
            end
        end

    end
end