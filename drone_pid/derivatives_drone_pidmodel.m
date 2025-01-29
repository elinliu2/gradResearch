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
        function dfdp = dfdp(model, state, t)

            dfdp = zeros(16, 12);
            dfdp(10, 1) = model.l/model.Ja * 2 * (pi/6*sin(t) - state(7));
            dfdp(10, 2) = model.l/model.Ja * 2 * state(13);
            dfdp(10, 3) = model.l/model.Ja * 2 * (pi/6*cos(t) - state(10));
    
            dfdp(11, 4) = model.l/model.Jb * 2 * (pi/6*cos(t) - state(8));
            dfdp(11, 5) = model.l/model.Jb * 2 * state(14);
            dfdp(11, 6) = model.l/model.Jb * 2 * (-pi/6*sin(t) - state(11));
    
            dfdp(12, 7) = model.K_psi/model.Jc * 4 * (t - state(9));
            dfdp(12, 8) = model.K_psi/model.Jc * 4 * state(15);
            dfdp(12, 9) = model.K_psi/model.Jc * 4 * (1 - state(12));
            
        end
        function d2fdx2 = d2fdx2(model, state)
            d2fdx2 = zeros(16, 16, 16);
            
            x7 = state(7);
            x8 = state(8);
            x9 = state(9);

            u_thrust = 3.961;

            d2fdx2(4, 7, 7) = 4*u_thrust*-sin(x7)*cos(x8)*cos(x9)/model.m;
            d2fdx2(4, 7, 8) = 4*u_thrust*cos(x7)*-sin(x8)*cos(x9)/model.m;
            d2fdx2(4, 7, 9) = 4*u_thrust*cos(x7)*cos(x8)*-sin(x9)/model.m;


            d2fdx2(4, 8, 7) = 4*u_thrust*(-cos(x7)*sin(x8)*cos(x9))/model.m;
            d2fdx2(4, 8, 8) = 4*u_thrust*(-sin(x8)*sin(x9) - sin(x7)*cos(x8)*cos(x9))/model.m;
            d2fdx2(4, 8, 9) = 4*u_thrust*(cos(x8)*cos(x9) - sin(x7)*sin(x8)*-sin(x9))/model.m;

            d2fdx2(4, 9, 7) = 4*u_thrust*(-cos(x7)*cos(x8)*sin(x9))/model.m;
            d2fdx2(4, 9, 8) = 4*u_thrust*(cos(x8)*cos(x9) - sin(x7)*-sin(x8)*sin(x9))/model.m;
            d2fdx2(4, 9, 9) = 4*u_thrust*(sin(x8)*-sin(x9) - sin(x7)*cos(x8)*cos(x9))/model.m;

            d2fdx2(5, 7, 7) = 4*u_thrust*-sin(x7)*cos(x8)*sin(x9)/model.m;
            d2fdx2(5, 7, 8) = 4*u_thrust*cos(x7)*-sin(x8)*sin(x9)/model.m;
            d2fdx2(5, 7, 9) = 4*u_thrust*cos(x7)*cos(x8)*cos(x9)/model.m;

            d2fdx2(5, 8, 7) = 4*u_thrust*(-cos(x7)*sin(x8)*sin(x9))/model.m;
            d2fdx2(5, 8, 8) = 4*u_thrust*(sin(x8)*cos(x9) - sin(x7)*cos(x8)*sin(x9))/model.m;
            d2fdx2(5, 8, 9) = 4*u_thrust*(-cos(x8)*-sin(x9) - sin(x7)*sin(x8)*cos(x9))/model.m;

            d2fdx2(5, 9, 7) = 4*u_thrust*( cos(x7)*cos(x8)*cos(x9))/model.m;
            d2fdx2(5, 9, 8) = 4*u_thrust*( cos(x8)*sin(x9) + sin(x7)*-sin(x8)*cos(x9))/model.m;
            d2fdx2(5, 9, 9) = 4*u_thrust*( sin(x8)*cos(x9) + sin(x7)*cos(x8)*-sin(x9))/model.m;

            d2fdx2(6, 7, 7) = 4*u_thrust*-cos(x7)*cos(x8)/model.m;
            d2fdx2(6, 7, 8) = 4*u_thrust*-sin(x7)*-sin(x8)/model.m;

            d2fdx2(6, 8, 7) = 4*u_thrust*-sin(x7)*-sin(x8)/model.m;
            d2fdx2(6, 8, 8) = 4*u_thrust*cos(x7)*-cos(x8)/model.m;
        end
        function d2fdxdp = d2fdxdp(model, ~)
            d2fdxdp = zeros(16, 16, 12);

            d2fdxdp(10, 7, 1) = model.l * 2 * -1/model.Ja;
            d2fdxdp(10, 13, 2) = model.l * 2 /model.Ja;
            d2fdxdp(10, 10, 3) = model.l * 2  * -1/model.Ja;
            
            
            d2fdxdp(11, 8, 4) = model.l * 2  * -1/model.Jb;
            d2fdxdp(11, 14, 5) = model.l * 2 /model.Jb;
            d2fdxdp(11, 11, 6) = model.l * 2 * -1/model.Jb;
            

            d2fdxdp(12, 9, 7) = model.K_psi * 4  * -1/model.Jc;
            d2fdxdp(12, 15, 8) = model.K_psi * 4 /model.Jc;
            d2fdxdp(12, 12, 9) = model.K_psi * 4  * -1/model.Jc;
        end
        function model = nudgeParams(model, delta, index)
            if index == 1
                model.kp1 = model.kp1 + delta;
            elseif index == 2
                model.ki1 = model.ki1 + delta;
            elseif index == 3
                model.kd1 = model.kd1 + delta;
            elseif index == 4
                model.kp2 = model.kp2 + delta;
            elseif index == 5
                model.ki2 = model.ki2 + delta;
            elseif index == 6
                model.kd2 = model.kd2 + delta;
            elseif index == 7
                model.kp3 = model.kp3 + delta;
            elseif index == 8
                model.ki3 = model.ki3 + delta;
            elseif index == 9
                model.kd3 = model.kd3 + delta;
            elseif index == 10
                model.kp4 = model.kp4 + delta;
            elseif index == 11
                model.ki4 = model.ki4 + delta;
            elseif index == 12
                model.kd4 = model.kd4 + delta;
            end
            
        end
    end
end

