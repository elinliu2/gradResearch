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
            % scatter3(x, y, z, 10, time,'filled');
            plot3(x, y, z, 'k--', 'LineWidth', 0.5);
            hold on

            % Drone arm length
            L = 0.1;

            xf = x(length(time), 1);
            yf = y(length(time), 1);
            zf = z(length(time), 1);
        
            % Define arm endpoints (before rotation)
            arms = [ xf+L  yf    zf;
                     xf-L  yf    zf;
                     xf    yf+L  zf;
                     xf    yf-L  zf];
        
            % Rotation matrices
            Rx = @(a)[1 0 0; 0 cosd(a) -sind(a); 0 sind(a) cosd(a)];
            Ry = @(a)[cosd(a) 0 sind(a); 0 1 0; -sind(a) 0 cosd(a)];
            Rz = @(a)[cosd(a) -sind(a) 0; sind(a) cosd(a) 0; 0 0 1];
        
            % Combined rotation
            R = Rz(0) * Ry(0) * Rx(0);
        
            % Apply rotation
            arms_rot = (R * arms')';
        
            % Plot drone
            plot3([arms_rot(1,1) arms_rot(2,1)], [arms_rot(1,2) arms_rot(2,2)], [arms_rot(1,3) arms_rot(2,3)], 'k-', 'LineWidth', 3);
            plot3([arms_rot(3,1) arms_rot(4,1)], [arms_rot(3,2) arms_rot(4,2)], [arms_rot(3,3) arms_rot(4,3)], 'k-', 'LineWidth', 3);
        
            % Center and rotors
            plot3(xf, yf, zf, 'ko','MarkerFaceColor','k');
            for i=1:4
                plot3(arms_rot(i,1), arms_rot(i,2), arms_rot(i,3), 'ko','MarkerFaceColor','k','MarkerSize',7.5);
            end
        
            title(plotTitle);
            xlabel("x [m]");
            ylabel("y [m]");
            zlabel("z [m]");
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