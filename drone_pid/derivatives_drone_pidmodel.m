classdef derivatives_drone_pidmodel
    methods (Static)
        function dfdx = dfdx(model, state)
            dfdx = zeros(16, 16);

            dfdx(1, 4) = 1;
            dfdx(2, 5) = 1;
            dfdx(3, 6) = 1;

            x7 = state(7);
            x8 = state(8);
            x9 = state(9);

            u_thrust = 3.961;

            dfdx(4, 7) = 4*u_thrust*cos(x7)*cos(x8)*cos(x9)/model.m;
            dfdx(4, 8) = 4*u_thrust*(cos(x8)*sin(x9) - sin(x7)*sin(x8)*cos(x9))/model.m;
            dfdx(4, 9) = 4*u_thrust*(sin(x8)*cos(x9) - sin(x7)*cos(x8)*sin(x9))/model.m;

            dfdx(5, 7) = 4*u_thrust*cos(x7)*cos(x8)*sin(x9)/model.m;
            dfdx(5, 8) = 4*u_thrust*(-cos(x8)*cos(x9) - sin(x7)*sin(x8)*sin(x9))/model.m;
            dfdx(5, 9) = 4*u_thrust*( sin(x8)*sin(x9) + sin(x7)*cos(x8)*cos(x9))/model.m;

            dfdx(6, 7) = 4*u_thrust*-sin(x7)*cos(x8)/model.m;
            dfdx(6, 8) = 4*u_thrust*cos(x7)*-sin(x8)/model.m;

            dfdx(7, 10) = 1;
            dfdx(8, 11) = 1;
            dfdx(9, 12) = 1;

            dfdx(10, 7) = model.l * 2 * model.kp1 * -1/model.Ja;
            dfdx(10, 10) = model.l * 2 * model.kd1 * -1/model.Ja;
            dfdx(10, 13) = model.l * 2 * model.ki1/model.Ja;
            
            dfdx(11, 8) = model.l * 2 * model.kp2 * -1/model.Jb;
            dfdx(11, 11) = model.l * 2 * model.kd2 * -1/model.Jb;
            dfdx(11, 14) = model.l * 2 * model.ki2/model.Jb;

            dfdx(12, 9) = model.K_psi * 4 * model.kp3 * -1/model.Jc;
            dfdx(12, 12) = model.K_psi * 4 * model.kd3 * -1/model.Jc;
            dfdx(12, 15) = model.K_psi * 4 * model.ki3/model.Jc;

            dfdx(13, 7) = -1;
            dfdx(14, 8) = -1;
            dfdx(15, 9) = -1;

            dfdx(16, 3) = -1;
            
        end
        function dfdp = dfdp(model, state)
            
        end
        function d2fdx2 = d2fdx2(model, state)
            
        end
        function d2fdxdp = d2fdxdp(model, state)
            
        end
        function model = nudgeParams(model, delta, index)
            if index == 1
                model.kp = model.kp + delta;
            elseif index == 2
                model.ki = model.ki + delta;
            elseif index == 3
                model.kd = model.kd + delta;
            end
        end
    end
end

