classdef derivatives_drone_IB_syms_dfdx
    methods (Static)
        function dfdx = dfdx(model, state, ref)
            dfdx = zeros(model.numStates, model.numStates);
            dfdx(1, 4) = 1;
            dfdx(2, 5) = 1;
            dfdx(3, 6) = 1;

            dfdx(4, :) = derivatives_drone_IB_syms_dfdx.df4dx(model, state, ref);
            dfdx(5, :) = derivatives_drone_IB_syms_dfdx.df5dx(model, state, ref);
            dfdx(6, :) = derivatives_drone_IB_syms_dfdx.df6dx(model, state, ref);

            dfdx(7, 10) = 1;
            dfdx(8, 11) = 1;
            dfdx(9, 12) = 1;

            dfdx(10, :) = derivatives_drone_IB_syms_dfdx_10.df10dx(model, state, ref);
            dfdx(11, :) = derivatives_drone_IB_syms_dfdx_11.df11dx(model, state, ref);
            dfdx(12, :) = derivatives_drone_IB_syms_dfdx.df12dx(model, state, ref);

            dfdx(13, :) = derivatives_drone_IB_syms_dfdx.df13dx(model, state, ref);
            dfdx(14, :) = derivatives_drone_IB_syms_dfdx.df14dx(model, state, ref);
            dfdx(15, 9) = -1;
            dfdx(16, 3) = -1;
            dfdx(17, 1) = -1;
            dfdx(18, 2) = -1;
            dfdx(19, 20) = 1;
            dfdx(20, :) = derivatives_drone_IB_syms_dfdx_20.df20dx(model, state, ref);
            dfdx(21, 23) = 1;
            dfdx(22, 24) = 1;
            dfdx(23, :) = derivatives_drone_IB_syms_dfdx_23.df23dx(model, state, ref);
            dfdx(24, :) = derivatives_drone_IB_syms_dfdx_24.df24dx(model, state, ref);
        end

        function df4dx = df4dx(model, state, ~)
            df4dx = zeros(1, model.numStates);
            
            m = model.m;

            x7 = state(7);
            x8 = state(8);
            x9 = state(9);
            x19 = state(19);

            df4dx(7)=(4*x19*cos(x7)*cos(x8)*cos(x9))/m;
            df4dx(8)=(4*x19*(cos(x8)*sin(x9) - cos(x9)*sin(x7)*sin(x8)))/m;
            df4dx(9)=(4*x19*(cos(x9)*sin(x8) - cos(x8)*sin(x7)*sin(x9)))/m;
            df4dx(19)=(4*(sin(x8)*sin(x9) + cos(x8)*cos(x9)*sin(x7)))/m;
        end

        function df5dx = df5dx(model, state, ~)
            df5dx = zeros(1, model.numStates);
            
            m = model.m;

            x7 = state(7);
            x8 = state(8);
            x9 = state(9);
            x19 = state(19);

            df5dx(7)=(4*x19*cos(x7)*cos(x8)*sin(x9))/m;
            df5dx(8)=-(4*x19*(cos(x8)*cos(x9) + sin(x7)*sin(x8)*sin(x9)))/m;
            df5dx(9)=(4*x19*(sin(x8)*sin(x9) + cos(x8)*cos(x9)*sin(x7)))/m;
            df5dx(19)=-(4*(cos(x9)*sin(x8) - cos(x8)*sin(x7)*sin(x9)))/m;
        end

        function df6dx = df6dx(model, state, ~)
            df6dx = zeros(1, model.numStates);

            m = model.m;

            x7 = state(7);
            x8 = state(8);
            x19 = state(19);

            df6dx(7)=-(4*x19*cos(x8)*sin(x7))/m;
            df6dx(8)=-(4*x19*cos(x7)*sin(x8))/m;
            df6dx(19)=(4*cos(x7)*cos(x8))/m;
        end

        function df12dx = df12dx(model, ~, ~)
            df12dx  = zeros(1, model.numStates);

            c5 = model.c5;
            c6 = model.c6;
            lambda3 = model.lambda3;

            df12dx(9)=c5^2 - lambda3 - c5*(c5 + c6) - 1;
            df12dx(12)=- c5 - c6;
            df12dx(15)=lambda3*(c5 + c6) - c5*lambda3;
        end

        function df13dx = df13dx(model, ~, ~)
            df13dx = zeros(1, model.numStates);
            df13dx(7)=-1;
            df13dx(21)=1;
        end

        function df14dx = df14dx(model, ~, ~)
            df14dx = zeros(1, model.numStates);
            df14dx(8)=-1;
            df14dx(22)=1;
        end
    end
end