classdef derivatives 
    methods (Static)
        function dfdx = derivativef(SC, x1)
            dfdx = [0 1; -(SC.g/SC.l)*cos(x1) -SC.D/(SC.m*SC.l)];
        end

        function dfdp = derivativefdp(SC, x1, x2)
            dfdp = [0 0; 1/(SC.m*SC.l) ((SC.g*sin(x1) + SC.D*x2/SC.m - SC.u/SC.m)/(SC.l^2)) ];
        end
        
        function d2fdx2 = derivative2fdx2(SC, x1)
            d2fdx2 = zeros(2,2,2);
            d2fdx2(2, 1, 1) = SC.g/SC.l * sin(x1);
        end

        function d2fdxdp = derivative2fdxdp(SC, x1)
            d2fdxdp = zeros(2,2,2);
            d2fdxdp(2, 1, 2) = SC.g/(SC.l^2)*cos(x1);
            d2fdxdp(2, 2, 2) = SC.D/(SC.m*SC.l^2);

        end
        
        function trajSens = calculateTrajectorySensitivities(x1, x2, SC)
            trajSens = zeros(2, 2, length(x2));
            trajSens(:, :, 1) = eye(2);
            I = eye(2);
            for i = 2:length(x2)
                dfdx_plus = derivatives.derivativef(SC, x1(i));
                dfdx = derivatives.derivativef(SC, x1(i-1));
        
                dxdxo = trajSens(:, :, i-1);
                trajSens(:, :, i) =  (I-(SC.T/2)*dfdx_plus)\((I + (SC.T/2)*dfdx) * dxdxo);
            end
        end

        function trajSens = calculateTrajectorySensitivitiesParameters(x1, x2, SC)
            trajSens = zeros(2, 2, length(x2));
            I = eye(2);
            trajSens(:, :, 1) = 0;
            for i = 2:length(x2)
                dfdx = derivatives.derivativef(SC, x1(i-1));
                dfdx_plus = derivatives.derivativef(SC, x1(i));
                
                dfdp = derivatives.derivativefdp(SC, x1(i-1), x2(i-1));
                dfdp_plus = derivatives.derivativefdp(SC, x1(i), x2(i));
        
                dxdp = trajSens(:, :, i-1);
                trajSens(:, :, i) =  (I-(SC.T/2)*dfdx_plus)\((SC.T/2)*(dfdx*dxdp + dfdp + dfdp_plus) + dxdp);
            end
        end
        
        function d2TrajSens = calculateSecondOrderTrajectorySensitivities(x1, x2, trajSens, SC)
            d2TrajSens = zeros(2, 2, 2, length(x2));
            I = eye(2);
        
            for i = 2:length(x2)
                dfdx_plus = derivatives.derivativef(SC, x1(i));
        
                d2xdx0 = d2TrajSens(:, :, :, i-1);
        
                d2fdx2_plus =  derivatives.derivative2fdx2(SC, x1(i));
                dxdx0_plus = trajSens(:,:,i);
        
                d2fdx2 =  derivatives.derivative2fdx2(SC, x1(i-1));
                dxdx0 = trajSens(:,:,i-1);
        
                dfdx =  derivatives.derivativef(SC, x1(i-1));
                
                a = inv(I - (SC.T/2)*dfdx_plus);
                b = conts({d2fdx2_plus, dxdx0_plus, dxdx0_plus}, {[3 1],[2 1]}, []);
                c = conts({d2fdx2, dxdx0, dxdx0}, {[3 1],[2 1]}, []);
                d = conts({d2xdx0, dfdx}, {[1 2]}, []); % how to know which dimension
                
                e = d2xdx0 + (SC.T/2)*(b+c+d);

                d2TrajSens(:, :, :, i) = conts({a, e}, {[2 1]}, []);
            end
        end

        function d2TrajSens = calculateSecondOrderTrajectorySensitivitiesParams(x1, x2, dxdxo_list, dxdp_list, SC)
            d2TrajSens = zeros(2, 2, 2, length(x2));
            I = eye(2);

            for i = 2:length(x2)
                
                d2fdx2 = derivatives.derivative2fdx2(SC, x1(i-1));
                dxdp = dxdp_list(:, :, i-1);
                dxdx0 = dxdxo_list(:, :, i-1);
                d2fdxdp = derivatives.derivative2fdxdp(SC, x1(i-1));
                dfdx = derivatives.derivativef(SC, x1(i-1));
                d2xdxodp = d2TrajSens(:, :, :, i-1);


                d2fdx2_plus = derivatives.derivative2fdx2(SC, x1(i));
                dxdp_plus = dxdp_list(:, :, i);
                dxdx0_plus = dxdxo_list(:, :, i);
                d2fdxdp_plus = derivatives.derivative2fdxdp(SC, x1(i));
                dfdx_plus = derivatives.derivativef(SC, x1(i));
                
                a = inv(I - (SC.T/2)*dfdx_plus);
                b = conts({d2fdx2, dxdp, dxdx0}, {[3 1],[2 1]}, []);
                c = conts({d2fdxdp, dxdx0}, {[2 1]}, []);
                d = conts({dfdx, d2xdxodp}, {[2 1]}, []);

                e = conts({d2fdx2_plus, dxdp_plus, dxdx0_plus}, {[3 1],[2 1]}, []);
                f = conts({d2fdxdp_plus, dxdx0_plus}, {[2 1]}, []);
                g = (SC.T/2)*(b + c + d + e + f) + d2xdxodp;
                
                d2TrajSens(:, :, :, i) = conts({a, g}, {[2 1]}, []);
            end

            
        end
        
        function trajSensTest = calculateTrajectorySensitivitiesTest(i, j, SC)
            delta = 1e-5;
            sCtrl = stonkCtrllers;
            [x1_dx1_plus, x2_dx1_plus] = pendulum_trajectory(SC, i + delta, j, @sCtrl.dummyU, @sCtrl.dummyD, false);
            [x1_dx1_minus, x2_dx1_minus]  = pendulum_trajectory(SC, i - delta, j, @sCtrl.dummyU, @sCtrl.dummyD, false);
        
            [x1_dx2_plus, x2_dx2_plus] = pendulum_trajectory(SC, i, j + delta, @sCtrl.dummyU, @sCtrl.dummyD, false);
            [x1_dx2_minus, x2_dx2_minus]  = pendulum_trajectory(SC, i, j - delta, @sCtrl.dummyU, @sCtrl.dummyD, false);
        
            dx1dx1 = (x1_dx1_plus-x1_dx1_minus)/(2*delta);
            dx2dx1 = (x2_dx1_plus-x2_dx1_minus)/(2*delta);
        
            dx1dx2 = (x1_dx2_plus-x1_dx2_minus)/(2*delta);
            dx2dx2 = (x2_dx2_plus-x2_dx2_minus)/(2*delta);
            
            splats = splatting;
            trajSensTest = splats.slapPack(dx1dx1, dx1dx2, dx2dx1, dx2dx2);
        end

        function trajSensTest = calculateTrajectorySensitivitiesTestParams(i, j, xe_1, xe_2, SC)
            delta = 1e-5;
            sCtrl = stonkCtrllers;
            
            SC.u = i + delta;
            SC.l = j;
            [x1_du_plus, x2_du_plus] = pendulum_trajectory(SC, xe_1, xe_2, @sCtrl.dummyU, @sCtrl.dummyD, false);
            SC.u = i - delta;
            [x1_du_minus, x2_du_minus]  = pendulum_trajectory(SC, xe_1, xe_2, @sCtrl.dummyU, @sCtrl.dummyD, false);
        
            SC.u = i;
            SC.l = j + delta;
            [x1_dl_plus, x2_dl_plus] = pendulum_trajectory(SC, xe_1, xe_2, @sCtrl.dummyU, @sCtrl.dummyD, false);
            SC.l = j - delta;
            [x1_dl_minus, x2_dl_minus]  = pendulum_trajectory(SC, xe_1, xe_2, @sCtrl.dummyU, @sCtrl.dummyD, false);
        
            dx1du = (x1_du_plus-x1_du_minus)/(2*delta);
            dx2du = (x2_du_plus-x2_du_minus)/(2*delta);
        
            dx1dl = (x1_dl_plus-x1_dl_minus)/(2*delta);
            dx2dl = (x2_dl_plus-x2_dl_minus)/(2*delta);
            
            splats = splatting;
            trajSensTest = splats.slapPack(dx1du, dx1dl, dx2du, dx2dl);
        end
        
        function secondTrajSensTest = calculateSecondOrderTrajectorySensitivitiesTest(i, j, SC)
            delta = 1e-5;
            d = derivatives;
            sCtrl = stonkCtrllers;
            [x1_pdeltax1, x2_pdeltax1] = pendulum_trajectory(SC, i + delta, j, @sCtrl.dummyU, @sCtrl.dummyD, false);
            [x1_mdeltax1, x2_mdeltax1] = pendulum_trajectory(SC, i - delta, j, @sCtrl.dummyU, @sCtrl.dummyD, false);
            trajSens1_plus = d.calculateTrajectorySensitivities(x1_pdeltax1, x2_pdeltax1, SC);
            trajSens1_minus = d.calculateTrajectorySensitivities(x1_mdeltax1, x2_mdeltax1, SC);
            trajSens1 = (trajSens1_plus - trajSens1_minus)/(2*delta);
        
            [x1_pdeltax2, x2_pdeltax2] = pendulum_trajectory(SC, i, j + delta, @sCtrl.dummyU, @sCtrl.dummyD, false);
            [x1_mdeltax2, x2_mdeltax2] = pendulum_trajectory(SC, i, j - delta, @sCtrl.dummyU, @sCtrl.dummyD, false);
            trajSens2_plus = d.calculateTrajectorySensitivities(x1_pdeltax2, x2_pdeltax2, SC);
            trajSens2_minus = d.calculateTrajectorySensitivities(x1_mdeltax2, x2_mdeltax2, SC);
            trajSens2 = (trajSens2_plus - trajSens2_minus)/(2*delta);
        
            splats = splatting;
            secondTrajSensTest = splats.slapPile(trajSens1, trajSens2);
        end

        function secondTrajSensTest = calculateSecondOrderTrajectorySensitivitiesTestParams(i, j, SC)
            delta = 1e-5;
            d = derivatives;
            sCtrl = stonkCtrllers;
            ICx = 0;
            ICy = 0;

            SC.u = i+delta;
            SC.l = j;
            [x1_pdeltax1, x2_pdeltax1] = pendulum_trajectory(SC, ICx, ICy, @sCtrl.dummyU, @sCtrl.dummyD, false);
            SC.u = i-delta;
            [x1_mdeltax1, x2_mdeltax1] = pendulum_trajectory(SC, ICx, ICy, @sCtrl.dummyU, @sCtrl.dummyD, false);
            trajSens1_plus = d.calculateTrajectorySensitivities(x1_pdeltax1, x2_pdeltax1, SC);
            trajSens1_minus = d.calculateTrajectorySensitivities(x1_mdeltax1, x2_mdeltax1, SC);
            trajSens1 = (trajSens1_plus - trajSens1_minus)/(2*delta);
        
            SC.u = i;
            SC.l = j+delta;
            [x1_pdeltax2, x2_pdeltax2] = pendulum_trajectory(SC, ICx, ICy, @sCtrl.dummyU, @sCtrl.dummyD, false);
            SC.l = j-delta;
            [x1_mdeltax2, x2_mdeltax2] = pendulum_trajectory(SC, ICx, ICy, @sCtrl.dummyU, @sCtrl.dummyD, false);
            trajSens2_plus = d.calculateTrajectorySensitivities(x1_pdeltax2, x2_pdeltax2, SC);
            trajSens2_minus = d.calculateTrajectorySensitivities(x1_mdeltax2, x2_mdeltax2, SC);
            trajSens2 = (trajSens2_plus - trajSens2_minus)/(2*delta);
        
            splats = splatting;
            secondTrajSensTest = splats.slapPile(trajSens1, trajSens2);

        end

        function [xe_1, xe_2] = calculateEqmPoint(SC)
            xe_1 = asin(SC.u/(SC.g*SC.m));
            xe_2 = 0;
        end

        function tp = findtp(trajSens)
            min = 1/norm(trajSens(:,:,1), 1);
            tp = 1;
            for i = 2:length(trajSens)
                if (1/norm(trajSens(:,:,i), 1) < min)
                    min = 1./norm(trajSens(:,:,i), 1);
                    tp = i;
                end
            end
        end

        function dG = calculate_dG(G, trajSens, d2xdxodp, tp)
            sum = zeros(1,2);
            for k = 1:2 
                for l = 1:2
                    sum = sum + sign(trajSens(k,l,tp))* d2xdxodp(:,k,l,tp);
                end
            end
            dG = -1 * G^2 * sum;
        end
        
        function dGtest = calculate_dGtest(i, j, IC, SC)
            delta = 1e-5;
            d = derivatives;
            n = normalStuff;
            sCtrl = stonkCtrllers;
            [x1_pdeltax1, x2_pdeltax1] = pendulum_trajectory(SC, i + delta, j, @sCtrl.dummyU, @sCtrl.dummyD, false);
            [x1_mdeltax1, x2_mdeltax1] = pendulum_trajectory(SC, i - delta, j, @sCtrl.dummyU, @sCtrl.dummyD, false);
            trajSens1_plus = d.calculateTrajectorySensitivities(x1_pdeltax1, x2_pdeltax1, IC, SC);
            trajSens1_minus = d.calculateTrajectorySensitivities(x1_mdeltax1, x2_mdeltax1, IC, SC);
            trajSens1 = (n.infimum(1./trajSens1_plus) - n.infimum(1./trajSens1_minus))/(2*delta);
        
            [x1_pdeltax2, x2_pdeltax2] = pendulum_trajectory(SC, i, j + delta, @sCtrl.dummyU, @sCtrl.dummyD, false);
            [x1_mdeltax2, x2_mdeltax2] = pendulum_trajectory(SC, i, j - delta, @sCtrl.dummyU, @sCtrl.dummyD, false);
            trajSens2_plus = d.calculateTrajectorySensitivities(x1_pdeltax2, x2_pdeltax2, IC, SC);
            trajSens2_minus = d.calculateTrajectorySensitivities(x1_mdeltax2, x2_mdeltax2, IC, SC);
            trajSens2 = (n.infimum(1./trajSens2_plus) - n.infimum(1./trajSens2_minus))/(2*delta);
        
            dGtest(1) = trajSens1;
            dGtest(2) = trajSens2;
        end
    end
end