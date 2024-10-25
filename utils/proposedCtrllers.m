classdef proposedCtrllers
    methods (Static)
        function u = controller(x1_o, x2_o, e, ed, ei, SC)

            deriv = derivatives;
            norms = normalStuff;
            sCtrl = stonkCtrllers;
        
            [x1, x2] = pendulum_trajectory(SC, x1_o, x2_o, @sCtrl.dummyU, @sCtrl.dummyD, false);
            trajSens = deriv.calculateTrajectorySensitivities(x1, x2, SC);
            Gx = norms.infimumOfReciprocalOfOneNorm(trajSens);
            
            trajSensTestParam = deriv.calculateTrajectorySensitivitiesTestParams(SC.u, SC.l, x1_o, x2_o, SC);
        
            secondTrajSensParam = deriv.calculateSecondOrderTrajectorySensitivitiesTestParams(SC.u, SC.l, x1_o, x2_o, SC);
           
            tp = deriv.findtp(trajSensTestParam);
            dGp = deriv.calculate_dG(trajSensTestParam, secondTrajSensParam, tp);
        
            
            % if G = 0 we want to heavily weigh the stabilizing controller 
            sf = 1;
            u = (sf/(Gx + sf)) * stabilityController(Gx, dGp)  ...
                + (1 - sf/(Gx + sf)) * PIDCtrl(x1_o, x2_o, e, ed, ei, SC);
        
        end
        
        function u = idkWhatImDoing(x1_o, x2_o, ~, ~, ~, SC)
            deriv = derivatives;
            norms = normalStuff;
            sCtrl = stonkCtrllers;
            oldSC = SC;
            SC.u = 0;
            SC.disturbance = 0;
        
            %Gx
            [x1, x2] = pendulum_trajectory(SC, x1_o, x2_o, @sCtrl.dummyU, @sCtrl.dummyD, false);
            trajSens = deriv.calculateTrajectorySensitivities(x1, x2, SC);
            Gx = norms.infimumOfReciprocalOfOneNorm(trajSens);
        
            if Gx < SC.GSEP
                % calcuate x1_plus, x2_plus, for one time stamp when I increase u + 1
                SC.u = 1;
                SC.t_final = SC.T*2;
                [x1plus, x2plus] = pendulum_trajectory(SC, x1_o, x2_o, @sCtrl.dummyU, @sCtrl.dummyD, false);
            
                SC.u = -1;
                [x1neg, x2neg] = pendulum_trajectory(SC, x1_o, x2_o, @sCtrl.dummyU, @sCtrl.dummyD, false);
            
                SC.u = 0; % calculate trajectory sensitivities with u = 0?
                SC.t_final = oldSC.t_final;
                [x1, x2] = pendulum_trajectory(SC, x1plus(2), x2plus(2), @sCtrl.dummyU, @sCtrl.dummyD, false);
                trajSens = deriv.calculateTrajectorySensitivities(x1, x2, SC); 
                Gxplus = norms.infimumOfReciprocalOfOneNorm(trajSens);
            
                [x1, x2] = pendulum_trajectory(SC, x1neg(2), x2neg(2), @sCtrl.dummyU, @sCtrl.dummyD, false);
                trajSens = deriv.calculateTrajectorySensitivities(x1, x2, SC); 
                Gxneg = norms.infimumOfReciprocalOfOneNorm(trajSens);
            
                if Gxplus > Gx
                    u = 1/Gx;
                elseif Gxneg > Gx
                    u = -1/Gx;
                else 
                    u = 0;
                end 
            else 
                u = 0;
            end
        end
        
        function u = idkWhatImDoingv2(x1_o, x2_o, ~, ~, ~, SC)
            deriv = derivatives;
            norms = normalStuff;
            sCtrl = stonkCtrllers;
            oldSC = SC;
        
            %Gx
            SC.u = 0; % calculate trajectory sensitivities with u = 0?
            [x1, x2] = pendulum_trajectory(SC, x1_o, x2_o, @sCtrl.dummyU, @sCtrl.dummyD, false);
            trajSens = deriv.calculateTrajectorySensitivities(x1, x2, SC);
            Gx1 = norms.infimumOfReciprocalOfOneNorm(trajSens);

            if ~isfield(SC, "GSEP")
                % G(SEP)
                [x1SEP, x2SEP] = pendulum_trajectory(SC, 0, 0, @sCtrl.dummyU, @sCtrl.dummyD, false);
                trajSens = deriv.calculateTrajectorySensitivities(x1SEP, x2SEP, SC);
                SC.GSEP = norms.infimumOfReciprocalOfOneNorm(trajSens);
            end
        
            if Gx1 < SC.GSEP
                % calcuate x1_plus, x2_plus, for one time stamp when I increase u + 1
                SC.u = 1;
                SC.t_final = SC.T*2;
            
                [x1plus, x2plus] = pendulum_trajectory(SC, x1_o, x2_o, @sCtrl.dummyU, @sCtrl.dummyD, false);
                SC.u = 0; % calculate trajectory sensitivities with u = 0?
                SC.t_final = oldSC.t_final;
                [x1, x2] = pendulum_trajectory(SC, x1plus(2), x2plus(2), @sCtrl.dummyU, @sCtrl.dummyD, false);
                trajSens = deriv.calculateTrajectorySensitivities(x1, x2, SC); 
                Gx2 = norms.infimumOfReciprocalOfOneNorm(trajSens);
                
                u = proposedCtrllers.findMaxU(x1_o, x2_o, SC, 0, 1, Gx1, Gx2);
                
            else
                u=0;
            end
            
        
        end
        
        
        function u = stabilityController(Gx, dGdp)
           if dGdp(1) < 0 
                u = -1/Gx;
            else 
                u = 1/Gx;
           end
        end
        
        function u = PIDCtrl(~, ~, e, ed, ei, ~)
            p = 2;
            i = 0.05;
            d = 2;
            u = p * e + i * ei + d * ed;
        end
        
        function u = findMaxU(x1_o, x2_o, SC, u1, u2, Gx1, Gx2)
            deriv = derivatives;
            norms = normalStuff;
            sCtrl = stonkCtrllers;
            oldSC = SC;
        
            if Gx2 > Gx1 % Go in the direction of u2
                uold = u2;
                Gxold = Gx2;
                u3 = u2 + 2*(u2-u1);
            elseif Gx1 > Gx2 % Go in the direction of u1
                uold = u1;
                Gxold = Gx1;
                u3 = u1 + 2*(u1-u2);
            else 
                u = u1;
                return;
            end
        
            SC.u = u3;
            SC.t_final = SC.T*2;
            [x1_3, x2_3] = pendulum_trajectory(SC, x1_o, x2_o, @sCtrl.dummyU, @sCtrl.dummyD, false);
            SC.u = 0; % calculate trajectory sensitivities with u = 0?
            SC.t_final = oldSC.t_final;
            [x1, x2] = pendulum_trajectory(SC, x1_3(2), x2_3(2), @sCtrl.dummyU, @sCtrl.dummyD, false);
            trajSens = deriv.calculateTrajectorySensitivities(x1, x2, SC); 
            Gx3 = norms.infimumOfReciprocalOfOneNorm(trajSens);
            
            if Gx3 > Gx2 && Gx3 > Gx1
                % Keep going in the direction we were going
                u = proposedCtrllers.findMaxU(x1_o, x2_o, SC, uold, u3, Gxold, Gx3);
            else 
                % binary search in the range we were just in
                u = proposedCtrllers.binarySearchMax(x1_o, x2_o, SC, uold, u3, Gxold, Gx3);
            end
        
        end
        
        % yes i know i could do something more sophisticated... in the future lol
        function u = binarySearchMax(x1_o, x2_o, SC, u1, u2, Gx1, Gx2) 
            tol = 0.05;
            deriv = derivatives;
            norms = normalStuff;
            sCtrl = stonkCtrllers;
        
            if Gx1/Gx2 < (1+tol) && Gx1/Gx2 > (1-tol)
                if Gx1 > Gx2
                    u = u1;
                else
                    u = u2;
                end
                return
            else
                t_final = SC.t_final;
                SC.u = (u1+u2)/2;
                SC.t_final = SC.T*2;
                [x1_3, x2_3] = pendulum_trajectory(SC, x1_o, x2_o, @sCtrl.dummyU, @sCtrl.dummyD, false);
                SC.u = 0; % calculate trajectory sensitivities with u = 0?
                SC.t_final = t_final;
                [x1, x2] = pendulum_trajectory(SC, x1_3(2), x2_3(2), @sCtrl.dummyU, @sCtrl.dummyD, false);
                trajSens = deriv.calculateTrajectorySensitivities(x1, x2, SC); 
                Gx3 = norms.infimumOfReciprocalOfOneNorm(trajSens);
        
                if (Gx1>Gx2)
                    u = proposedCtrllers.binarySearchMax(x1_o, x2_o, SC, u1, (u1+u2)/2, Gx1, Gx3); 
                else
                    u = proposedCtrllers.binarySearchMax(x1_o, x2_o, SC, u2, (u1+u2)/2, Gx2, Gx3); 
                end
            end
        end
    end
end

