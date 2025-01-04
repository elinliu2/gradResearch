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
        function [F1, F2, F3, F4, model] = ParameterizedPIDCtrl(state, ref, model)
            delta_x = ref(1) - state(1);
            delta_y = ref(2) - state(2);
            delta_z = ref(3) - state(3);

            if(delta_z > 0.6)
                theta_ref = atan2(delta_x, delta_z);
                phi_ref = atan2(-delta_y, delta_z);
            else 
                theta_ref = 0;
                phi_ref = 0;
            end

            psi_ref = 0; 
            theta_ref_dot = (1/(1+( delta_x/delta_z)^2)) * ( delta_z*-state(4) - ( delta_x)*-state(6) ) / delta_z^2 ;
            phi_ref_dot =   (1/(1+(-delta_y/delta_z)^2)) * ( delta_z*state(5) - ( -delta_y)*-state(6) ) / delta_z^2 ;

            % thrust_err = norm(ref-state(1:3));
            % thrust_err_dot = -1/sqrt(thrust_err) * (delta_x * state(4) + delta_y * state(5) + delta_z * state(6));

            thrust_err = delta_z;
            thrust_err_dot = -state(6);

            u_theta = model.kp1 * (theta_ref - state(7)) + model.ki1 * state(13) + model.kd1 * ( theta_ref_dot - state(10));
            u_phi = model.kp2 * (phi_ref - state(8)) + model.ki2 * state(14) + model.kd2 * ( phi_ref_dot - state(11));
            u_psi = model.kp3 * (psi_ref - state(9)) + model.ki3 * state(15) + model.kd3 * (-1 * state(12));
            u_thrust = model.kp4 * thrust_err + model.ki4 * state(16) + model.kd4 * thrust_err_dot;

            F1 = u_theta + u_psi + u_thrust;
            F2 = - u_theta + u_psi + u_thrust;
            F3 = u_phi - u_psi + u_thrust;
            F4 = - u_phi - u_psi + u_thrust;
        end
   end
end

