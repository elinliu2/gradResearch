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

        function StatePlot(x1, x2, time, figNum, plotTitle)
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
    end
end