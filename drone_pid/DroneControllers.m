classdef DroneControllers
   methods (Static)
        function [F1, F2, F3, F4, model] = dummyU(~, model)
            F1 = 3.45;
            F2 = 3.44;
            F3 = 3.45;
            F4 = 3.44;
        end
        function d = dummyD(disturbance, time)
            if (time < 5)
                d = disturbance;
            else 
                d = 0;
            end
        end
        function [F1, F2, F3, F4, model] = ParameterizedPIDCtrl(state, model)
            x_ref = 1;
            y_ref = 1; 
            z_ref = 101; 

            z_delta = z_ref - state(3);

            theta_ref = atan2(state(1) - x_ref, z_delta);
            phi_ref = atan2(state(2) - y_ref, z_delta);
            psi_ref = 0; 

            xz = (state(1) - x_ref) / (z_delta);
            yz = (state(2) - y_ref) / (z_delta);
            
            theta_ref_dot = (1/(1+xz^2)) * ( z_delta*state(4) - (state(1) - x_ref)*-state(6) ) / z_delta^2 ;
            phi_ref_dot =   (1/(1+yz^2)) * ( z_delta*state(5) - (state(2) - y_ref)*-state(6) ) / z_delta^2 ;

            u_theta = model.kp1 * (theta_ref - state(7)) + model.ki1 * state(13) + model.kd1 * ( theta_ref_dot - state(10));
            u_phi = model.kp2 * (phi_ref - state(8)) + model.ki2 * -1 * state(14) + model.kd2 * ( phi_ref_dot - state(11));
            u_psi = model.kp3 * (psi_ref - state(9)) + model.ki3 * -1 * state(15) + model.kd3 * (-1 * state(12));
            u_z = model.kp4 * (z_ref - state(3)) + model.ki4 * -1 * state(15) + model.kd4 * (-state(6));

            F1 = u_theta + u_psi + u_z;
            F2 = - u_theta + u_psi + u_z;
            F3 = u_phi - u_psi + u_z;
            F4 = - u_phi - u_psi + u_z;
        end
   end
end

