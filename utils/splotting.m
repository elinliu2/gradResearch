classdef splotting
    methods (Static)
        function splot(figNum, x1, x2, func, title)
            figure(figNum)
            scatter3(x1, x2, func, 10, func/max(func), "filled");
            hold on
            xlabel("x1");
            ylabel("x2");
            zlabel(title);
            hold off
        end
        
        function squiver(figNum, x1_val, x2_val, U, V, title)
            figure(figNum);
            hold off
            quiver(x1_val, x2_val, U, V);
            hold on
        
            xlabel("x1");
            ylabel("x2");
            zlabel(title);
        end

        function StatePlot(time, x1, x2, figNum, plotTitle)
            figure(figNum);
            hold off
            scatter(x1, x2, 1, time);
            hold on
            title(plotTitle);
            xlabel("x1");
            ylabel("x2");
            % xlim([-8 8]);
            % ylim([-8 8]);
        end

        function uplup(time, u, figNum, plotTitle, ylab)
            figure(figNum);
            hold off
            plot(time, u);
            hold on
            title(plotTitle);
            xlabel("time");
            ylabel(ylab);
        end

        function trajectoryPlot3d(time, x, y, z, figNum, plotTitle)
            figure(figNum);
            hold off
            scatter3(x, y, z, 10, time);
            hold on
            title(plotTitle);
            xlabel("x");
            ylabel("y");
            zlabel("z");
        end

        function plotDroneXYZ(time, stateProgression, figNum)
            figure(figNum);
            subplot(3, 1, 1);
            plot(time,stateProgression(:,1));
            subplot(3, 1, 2);
            plot(time,stateProgression(:,2));
            subplot(3, 1, 3);
            plot(time,stateProgression(:,3));
        end
        
        function plotDroneAngles(time, stateProgression, figNum)
            figure(figNum);
            subplot(3, 1, 1);
            plot(time,stateProgression(:,7));
            subplot(3, 1, 2);
            plot(time,stateProgression(:,8));
            subplot(3, 1, 3);
            plot(time,stateProgression(:,9));
        end

    end
end