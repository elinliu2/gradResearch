classdef DroneControllers
   methods (Static)
        function [F1, F2, F3, F4, model] = dummyU(~, ~, model)
            F1 = 3.44;
            F2 = 3.43;
            F3 = 3.43;
            F4 = 3.43;
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
          
            psi_ref = 0; 

            theta_ref = atan2(delta_x, delta_z);
            phi_ref = atan2(-delta_y, delta_z);

            theta_ref_dot = (1/(1+( delta_x/delta_z)^2)) * ( delta_z*-state(4) - ( delta_x)*-state(6) ) / delta_z^2 ;
            phi_ref_dot =   (1/(1+( -delta_y/delta_z)^2)) * ( delta_z*state(5) - (-delta_y)*-state(6) ) / delta_z^2 ;
            
            thrust_err = delta_z;
            thrust_err_dot = -state(6);

            U_theta = model.kp1 * (theta_ref - state(7)) + model.ki1 * state(13) + model.kd1 * ( theta_ref_dot - state(10));
            U_phi = model.kp2 * (phi_ref - state(8)) + model.ki2 * state(14) + model.kd2 * ( phi_ref_dot - state(11));
            U_psi = model.kp3 * (psi_ref - state(9)) + model.ki3 * state(15) + model.kd3 * (-1 * state(12));
            U_thrust = model.kp4 * thrust_err + model.ki4 * state(16) + model.kd4 * thrust_err_dot;


            F1 = U_theta + U_psi + U_thrust;
            F2 = - U_theta + U_psi + U_thrust;
            F3 = U_phi - U_psi + U_thrust;
            F4 = - U_phi - U_psi + U_thrust;
        end

        function [F1, F2, F3, F4, model] = CircleControl(state, t, model)
            % quat = normalize(quaternion([0 pi/6 0], 'euler', 'ZYX', 'frame'));
            % quat_delta = normalize(quaternion([0 0 t], 'euler', 'ZYX', 'frame'));
            % new_quat = quat_delta*quat;
            % eul = quat2eul(new_quat);
            % 
            % quat_delta_old = normalize(quaternion([0 0 t-1e-3], 'euler', 'ZYX', 'frame'));
            % new_quat_old = quat_delta_old*quat;
            % eul_old = quat2eul(new_quat_old);
            % 
            % theta_ref = eul(1);
            % phi_ref = eul(2);
            % psi_ref = abs(eul(3));
            % 
            % eul_deriv = (eul-eul_old)/1e-3;
            % 
            % eul_deriv(3) = abs(eul(3)) - abs(eul_old(3));
            % 
            % theta_ref_dot = eul_deriv(1);
            % phi_ref_dot = eul_deriv(2);            
            % psi_ref_dot = eul_deriv(3);

            theta_ref = pi/6;
            phi_ref = 0;
            psi_ref = t;

            theta_ref_dot = 0;
            phi_ref_dot = 0;
            psi_ref_dot = 1;

            u_theta = model.kp1 * (theta_ref - state(7)) + model.ki1 * state(13) + model.kd1 * ( theta_ref_dot - state(10));
            u_phi = model.kp2 * (phi_ref - state(8)) + model.ki2 * state(14) + model.kd2 * ( phi_ref_dot - state(11));
            u_psi = model.kp3 * (psi_ref - state(9)) + model.ki3 * state(15) + model.kd3 * (psi_ref_dot - state(12));
      
            u_thrust = 3.961;

            F1 = u_theta + u_psi + u_thrust;
            F2 = - u_theta + u_psi + u_thrust;
            F3 = u_phi - u_psi + u_thrust;
            F4 = - u_phi - u_psi + u_thrust;
        end

        function [F1, F2, F3, F4, model] = pleaseJustGoInAStraightLine(state, ~, model)
            delta_z = 101 - state(3);

            theta_ref = pi/4;
            phi_ref = pi/4;
            psi_ref = 0;

            theta_ref_dot = 0;
            phi_ref_dot = 0;
            psi_ref_dot = 0;

            thrust_err = delta_z;
            thrust_err_dot = -state(6);

            u_theta = model.kp1 * (theta_ref - state(7)) + model.ki1 * state(13) + model.kd1 * ( theta_ref_dot - state(10));
            u_phi = model.kp2 * (phi_ref - state(8)) + model.ki2 * state(14) + model.kd2 * ( phi_ref_dot - state(11));
            u_psi = model.kp3 * (psi_ref - state(9)) + model.ki3 * state(15) + model.kd3 * (psi_ref_dot - state(12));
            u_thrust = model.kp4 * thrust_err + model.ki4 * state(16) + model.kd4 * thrust_err_dot;

            F1 = u_theta + u_psi + u_thrust;
            F2 = - u_theta + u_psi + u_thrust;
            F3 = u_phi - u_psi + u_thrust;
            F4 = - u_phi - u_psi + u_thrust;
        end

        function [F1, F2, F3, F4, model] = droneGoBrrr(state, t, model)
            delta_x = (0.1*sin(0.1*t) - state(1));
            delta_y = -(0.1*cos(0.1*t) - state(2));
            delta_z = 101 - state(3);

            delta_x_dot = (0.01*cos(t) - state(4));
            delta_y_dot = -(-0.01*sin(t) - state(5));

            psi_ref = 0;
            psi_ref_dot = 0;

            thrust_err = delta_z;
            thrust_err_dot = -state(6);

            u_theta = model.kp1 * delta_x + model.ki1 * state(13) + model.kd1 * delta_x_dot;
            u_phi = model.kp2 * delta_y + model.ki2 * state(14) + model.kd2 * delta_y_dot;
            u_psi = model.kp3 * (psi_ref - state(9)) + model.ki3 * state(15) + model.kd3 * (psi_ref_dot - state(12));
            u_thrust = 5;
            
            F1 = u_theta + u_psi + u_thrust;
            F2 = - u_theta + u_psi + u_thrust;
            F3 = u_phi - u_psi + u_thrust;
            F4 = - u_phi - u_psi + u_thrust;
        end

   end
end

