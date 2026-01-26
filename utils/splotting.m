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
            scatter3(x, y, z, 10, time,'filled');
            % plot3(x, y, z, 'k--', 'LineWidth', 0.5);
            % hold on
            % 
            % % Drone arm length
            % L = 0.1;
            % 
            % xf = x(length(time), 1);
            % yf = y(length(time), 1);
            % zf = z(length(time), 1);
            % 
            % % Define arm endpoints (before rotation)
            % arms = [ xf+L  yf    zf;
            %          xf-L  yf    zf;
            %          xf    yf+L  zf;
            %          xf    yf-L  zf];
            % 
            % % Rotation matrices
            % Rx = @(a)[1 0 0; 0 cosd(a) -sind(a); 0 sind(a) cosd(a)];
            % Ry = @(a)[cosd(a) 0 sind(a); 0 1 0; -sind(a) 0 cosd(a)];
            % Rz = @(a)[cosd(a) -sind(a) 0; sind(a) cosd(a) 0; 0 0 1];
            % 
            % % Combined rotation
            % R = Rz(0) * Ry(0) * Rx(0);
            % 
            % % Apply rotation
            % arms_rot = (R * arms')';
            % 
            % % Plot drone
            % plot3([arms_rot(1,1) arms_rot(2,1)], [arms_rot(1,2) arms_rot(2,2)], [arms_rot(1,3) arms_rot(2,3)], 'k-', 'LineWidth', 3);
            % plot3([arms_rot(3,1) arms_rot(4,1)], [arms_rot(3,2) arms_rot(4,2)], [arms_rot(3,3) arms_rot(4,3)], 'k-', 'LineWidth', 3);
            % 
            % % Center and rotors
            % plot3(xf, yf, zf, 'ko','MarkerFaceColor','k');
            % for i=1:4
            %     plot3(arms_rot(i,1), arms_rot(i,2), arms_rot(i,3), 'ko','MarkerFaceColor','k','MarkerSize',7.5);
            % end
        
            title(plotTitle);
            xlabel("x [m]");
            ylabel("y [m]");
            zlabel("z [m]");
        end

        function trajectoryPlots()
            stateNoDist = load("stateInitialParams.mat").stateProgression;
            stateDistInitial = load("stateDistInitialParams.mat").stateProgression;
            stateDistNew = load("stateDistNewParams.mat").stateProgression;

            states = zeros(3, 100001, 23);
            states(1, :, :) = stateNoDist; 
            states(2, :, :) = stateDistInitial; 
            states(3, :, :) = stateDistNew; 
            hold off;
            plot3(states(1, 1, 1), states(1, 1, 2), states(1, 1, 3), 'o', ...
              'MarkerSize', 10, ...
              'MarkerFaceColor', 'g', ...
              'MarkerEdgeColor', 'k', ...
              'LineWidth', 0.5);
            hold on;
            for plot_index = 1:3
                x = states(plot_index, :, 1);
                y = states(plot_index, :, 2);
                z = states(plot_index, :, 3);

                if plot_index == 1 
                    plot3(x, y, z, 'r:', 'LineWidth', 1.5);
                    hold on;
                elseif plot_index == 2
                    plot3(x, y, z, 'k-', 'LineWidth', 1);
                else 
                    plot3(x, y, z, 'b--', 'LineWidth', 1.5);
                end
            
            end

            plot3(1, 1, 100, 'p', ...
              'MarkerSize', 15, ...
              'MarkerFaceColor', 'y', ...
              'MarkerEdgeColor', 'k', ...
              'LineWidth', 0.5);
        
            title("Drone Trajectory");
            xlabel("x [m]");
            xlim([-5 2]);
            ylabel("y [m]");
            ylim([-5 2]);
            zlabel("z [m]");
            zlim([99 101]);
            grid on;

            for plot_index = 1:3
                x = states(plot_index, :, 1);
                y = states(plot_index, :, 2);
                z = states(plot_index, :, 3);
                dx = diff(x);
                dy = diff(y);
                dz = diff(z);

                dirs = [dx(:), dy(:), dz(:)];
                dirs = dirs ./ vecnorm(dirs,2,2);

                L = 0.15;
                R = 0.25 * L;                                 % radius
                
                [CX, CY, CZ] = cylinder([R 0], 20);
                CZ = CZ * L;

                if plot_index == 1
                    % idx = 2000;
                    % for k = idx
                    %     p = [x(k+1), y(k+1), z(k+1)];
                    %     v = dirs(k,:);
                    % 
                    %     Rm = vrrotvec2mat(vrrotvec([0 0 1], v));
                    %     P  = Rm * [CX(:)'; CY(:)'; CZ(:)'];
                    % 
                    %     surf( reshape(P(1,:),size(CX)) + p(1), ...
                    %           reshape(P(2,:),size(CY)) + p(2), ...
                    %           reshape(P(3,:),size(CZ)) + p(3), ...
                    %           'FaceColor','r', ...
                    %           'EdgeColor','none');
                    % end
                elseif plot_index == 2
                    idx = [156, 167];
                    for k = idx
                        p = [x(k+1), y(k+1), z(k+1)];
                        v = dirs(k,:);
                    
                        Rm = vrrotvec2mat(vrrotvec([0 0 1], v));
                        P  = Rm * [CX(:)'; CY(:)'; CZ(:)'];
                    
                        surf( reshape(P(1,:),size(CX)) + p(1), ...
                              reshape(P(2,:),size(CY)) + p(2), ...
                              reshape(P(3,:),size(CZ)) + p(3), ...
                              'FaceColor','k', ...
                              'EdgeColor','none');
                    end
                else 
                    % idx = 500;
                    % for k = idx
                    %     p = [x(k+1), y(k+1), z(k+1)];
                    %     v = dirs(k,:);
                    % 
                    %     Rm = vrrotvec2mat(vrrotvec([0 0 1], v));
                    %     P  = Rm * [CX(:)'; CY(:)'; CZ(:)'];
                    % 
                    %     surf( reshape(P(1,:),size(CX)) + p(1), ...
                    %           reshape(P(2,:),size(CY)) + p(2), ...
                    %           reshape(P(3,:),size(CZ)) + p(3), ...
                    %           'FaceColor','b', ...
                    %           'EdgeColor','none');
                    % end
                end
            end

            for plot_index = 2:3
                 % Drone arm length
                L = 0.1;
                if plot_index == 3
                    xf = 0.63;
                    yf = 0.63;
                    zf = 99.95;
                else 
                    xf = -1.21;
                    yf = -2.79;
                    zf = 100.138;
                end
                % Define arm endpoints (before rotation)
                arms_body = [ L   0   0;
                             -L   0   0;
                              0   L   0;
                              0  -L   0];

                Rx = @(a)[1 0 0;
                          0 cosd(a) -sind(a);
                          0 sind(a)  cosd(a)];

                Ry = @(a)[ cosd(a) 0 sind(a);
                           0       1 0;
                          -sind(a) 0 cosd(a)];

                Rz = @(a)[cosd(a) -sind(a) 0;
                          sind(a)  cosd(a) 0;
                          0        0       1];

                R = Rz(90) * Ry(30) * Rx(30);

                arms_rot   = (R * arms_body')';
                arms_world = arms_rot + [xf yf zf];

                % arms
                plot3(arms_world([1 2],1),arms_world([1 2],2),arms_world([1 2],3),'k-','LineWidth',3)
                plot3(arms_world([3 4],1),arms_world([3 4],2),arms_world([3 4],3),'k-','LineWidth',3)

                % center (hollow)
                plot3(xf,yf,zf,'ko','MarkerSize',6,'LineWidth',1.5)

                for i=1:4
                    plot3(arms_world(i,1), arms_world(i,2), arms_world(i,3), 'ko', 'MarkerSize',10);
                end
            end

            

            legend("Initial Condition", "Initial Parameters No Disturbance", "Initial Parameters With Disturbance", "Locally Optimal Parameters With Disturbance", "Goal Point")
        end

        function trajectoryPlotsstateDistInitial()
            stateDistInitial = load("stateDistInitialParams.mat").stateProgression;
            plot3(stateDistInitial(1, 1), stateDistInitial(1, 2), stateDistInitial(1, 3), 'o', ...
              'MarkerSize', 10, ...
              'MarkerFaceColor', 'g', ...
              'MarkerEdgeColor', 'k', ...
              'LineWidth', 0.5);
            hold on;
            x = stateDistInitial(1:300, 1);
            y = stateDistInitial(1:300, 2);
            z = stateDistInitial(1:300, 3);
            plot3(x, y, z, 'k-', 'LineWidth', 1.5);

            plot3(1, 1, 100, 'p', ...
              'MarkerSize', 10, ...
              'MarkerFaceColor', 'y', ...
              'MarkerEdgeColor', 'k', ...
              'LineWidth', 0.5);
        
            title("Drone Trajectory with Initial Parameters and Finite Time Disturbance");
            xlabel("x [m]");
            xlim([-20 5]);
            ylabel("y [m]");
            ylim([-30 5]);
            zlabel("z [m]");
            zlim([95 110]);
            grid on;

           
            % arrows
            x = stateDistInitial(:, 1);
            y = stateDistInitial(:, 2);
            z = stateDistInitial(:, 3);
            dx = diff(x);
            dy = diff(y);
            dz = diff(z);

            dirs = [dx(:), dy(:), dz(:)];
            dirs = dirs ./ vecnorm(dirs,2,2);

            L = 1.5;
            R = 0.25 * L;                                 % radius
            
            [CX, CY, CZ] = cylinder([R 0], 20);
            CZ = CZ * L;

            
            idx = [174, 201];
            for k = idx
                p = [x(k+1), y(k+1), z(k+1)];
                v = dirs(k,:);
            
                Rm = vrrotvec2mat(vrrotvec([0 0 1], v));
                P  = Rm * [CX(:)'; CY(:)'; CZ(:)'];
            
                surf( reshape(P(1,:),size(CX)) + p(1), ...
                      reshape(P(2,:),size(CY)) + p(2), ...
                      reshape(P(3,:),size(CZ)) + p(3), ...
                      'FaceColor','k', ...
                      'EdgeColor','none');
            end
               
            % drone
            xf = -8.28808;
            yf = -16.9257;
            zf = 102.987;
               
            % Define arm endpoints (before rotation)
            arms_body = [ L   0   0;
                         -L   0   0;
                          0   L   0;
                          0  -L   0];

            Rx = @(a)[1 0 0;
                      0 cosd(a) -sind(a);
                      0 sind(a)  cosd(a)];

            Ry = @(a)[ cosd(a) 0 sind(a);
                       0       1 0;
                      -sind(a) 0 cosd(a)];

            Rz = @(a)[cosd(a) -sind(a) 0;
                      sind(a)  cosd(a) 0;
                      0        0       1];

            R = Rz(0) * Ry(10) * Rx(10);

            arms_rot   = (R * arms_body')';
            arms_world = arms_rot + [xf yf zf];

            % arms
            plot3(arms_world([1 2],1),arms_world([1 2],2),arms_world([1 2],3),'k-','LineWidth',3)
            plot3(arms_world([3 4],1),arms_world([3 4],2),arms_world([3 4],3),'k-','LineWidth',3)

            % center (hollow)
            plot3(xf,yf,zf,'ko','MarkerSize',6,'LineWidth',1.5)

            for i=1:4
                plot3(arms_world(i,1), arms_world(i,2), arms_world(i,3), 'ko', 'MarkerSize',10);
            end
            legend("Start Point", "Initial Parameters With Disturbance", "Goal Point")
        end

        function trajectoryPlotsstateDistNew()
            stateDistNew = load("stateDistNewParams.mat").stateProgression;

            plot3(stateDistNew(1, 1), stateDistNew(1, 2), stateDistNew(1, 3), 'o', ...
              'MarkerSize', 10, ...
              'MarkerFaceColor', 'g', ...
              'MarkerEdgeColor', 'k', ...
              'LineWidth', 0.5);
            hold on;
            x = stateDistNew(:, 1);
            y = stateDistNew(:, 2);
            z = stateDistNew(:, 3);
            plot3(x, y, z, 'k-', 'LineWidth', 1.5);

            plot3(1, 1, 100, 'p', ...
              'MarkerSize', 10, ...
              'MarkerFaceColor', 'y', ...
              'MarkerEdgeColor', 'k', ...
              'LineWidth', 0.5);
        
            title("Drone Trajectory with Locally Optimal Parameters and Finite Time Disturbance");
            xlabel("x [m]");
            xlim([-0.5 2]);
            ylabel("y [m]");
            ylim([-0.5 2]);
            zlabel("z [m]");
            zlim([99 101]);
            grid on;

           
            % arrows
            x = stateDistNew(:, 1);
            y = stateDistNew(:, 2);
            z = stateDistNew(:, 3);
            dx = diff(x);
            dy = diff(y);
            dz = diff(z);

            dirs = [dx(:), dy(:), dz(:)];
            dirs = dirs ./ vecnorm(dirs,2,2);

            L = 0.2;
            R = 0.25 * L;                                 % radius
            
            [CX, CY, CZ] = cylinder([R 0], 20);
            CZ = CZ * L;

            
            idx = [1159, 3759];
            for k = idx
                p = [x(k+1), y(k+1), z(k+1)];
                v = dirs(k,:);
            
                Rm = vrrotvec2mat(vrrotvec([0 0 1], v));
                P  = Rm * [CX(:)'; CY(:)'; CZ(:)'];
            
                surf( reshape(P(1,:),size(CX)) + p(1), ...
                      reshape(P(2,:),size(CY)) + p(2), ...
                      reshape(P(3,:),size(CZ)) + p(3), ...
                      'FaceColor','k', ...
                      'EdgeColor','none');
            end
               
            % drone
            xf = 0.8378;
            yf = 0.7748;
            zf = 99.9446;
            L = 0.1;
            % Define arm endpoints (before rotation)
            arms_body = [ L   0   0;
                         -L   0   0;
                          0   L   0;
                          0  -L   0];

            Rx = @(a)[1 0 0;
                      0 cosd(a) -sind(a);
                      0 sind(a)  cosd(a)];

            Ry = @(a)[ cosd(a) 0 sind(a);
                       0       1 0;
                      -sind(a) 0 cosd(a)];

            Rz = @(a)[cosd(a) -sind(a) 0;
                      sind(a)  cosd(a) 0;
                      0        0       1];

            R = Rz(0) * Ry(10) * Rx(10);

            arms_rot   = (R * arms_body')';
            arms_world = arms_rot + [xf yf zf];

            % arms
            plot3(arms_world([1 2],1),arms_world([1 2],2),arms_world([1 2],3),'k-','LineWidth',3)
            plot3(arms_world([3 4],1),arms_world([3 4],2),arms_world([3 4],3),'k-','LineWidth',3)

            % center (hollow)
            plot3(xf,yf,zf,'ko','MarkerSize',6,'LineWidth',1.5)

            for i=1:4
                plot3(arms_world(i,1), arms_world(i,2), arms_world(i,3), 'ko', 'MarkerSize',10);
            end
            legend("Start Point", "Locally Optimal Parameters With Disturbance", "Goal Point")
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

        function trajectoryVideo(time, stateProgression)
            figure('Units','normalized','OuterPosition',[0 0 1 1]);
            % Load STL file
            fv = stlread('./drone.stl'); 
            scale = 0.02;

            % Apply scaling
            scaledPoints = fv.Points * scale;
        
            % === Create transform handle ===
            t = hgtransform;
        
            % === Attach STL patch to the transform ===
            h_drone = patch('Faces', fv.ConnectivityList, ...
                'Vertices', scaledPoints, ...
                'FaceColor',[0.3 0.3 0.8], ...
                'EdgeColor','none', ...
                'Parent', t);
        
            % === Lighting & figure settings ===
            view(3)               % force 3D perspective view
            grid on;
            camlight; lighting gouraud;
            title("Drone Sim");
            xlabel("x");
            ylabel("y");
            zlabel("z");
        
            % === Fixed reference marker ===
            hold on;
        
            % === Axis limits so both target + drone are visible ===
            xlim([min(stateProgression(:,1))-0.5, max(stateProgression(:,1))+0.5]);
            ylim([min(stateProgression(:,2))-0.5, max(stateProgression(:,2))+0.5]);
            zlim([min([stateProgression(:,3); 100])-1, max([stateProgression(:,3); 100])+1]);
 
            % === Annotation for time ===
            h_annot = annotation('textbox',[.9 .5 .1 .2], ...
                                 'String',0, ...
                                 'EdgeColor','none');

            writerObj = VideoWriter('drone.mp4','MPEG-4');
            writerObj.FrameRate=60;
            open(writerObj);

            for i = 1 :16:length(time)
                t_now = time(i);
                plot3(stateProgression(1:i,1), stateProgression(1:i,2), stateProgression(1:i,3), 'k--', 'LineWidth', 0.5);

                % --- Translation (move drone along trajectory)
                tx = stateProgression(i,1);
                ty = stateProgression(i,2);
                tz = stateProgression(i,3);
        
                yaw   = stateProgression(i,7);  
                pitch = stateProgression(i,8);
                roll  = stateProgression(i,9);
        
                % --- Build transformation matrix ---
                M = makehgtform('translate',[tx ty tz], ...
                                'zrotate',yaw, ...
                                'yrotate',pitch, ...
                                'xrotate',roll);
        
                % --- Apply to hgtransform ---
                t.Matrix = M;
        
                % --- Update time annotation ---
                set(h_annot,'String',sprintf('t=%.2f s',t_now));
                
                % pad = 0.5;  
                % xlim([tx - pad, tx + pad]);
                % ylim([ty - pad, ty + pad]);
                % zlim([tz - pad, tz + pad]);
                % --- Render and capture frame ---
                drawnow;
                frame = getframe(gcf);
                writeVideo(writerObj, frame);
            end
            
            close(writerObj);

        end

    end
end