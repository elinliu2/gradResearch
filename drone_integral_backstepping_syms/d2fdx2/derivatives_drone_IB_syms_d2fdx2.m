classdef derivatives_drone_IB_syms_d2fdx2
    methods (Static)
        function d2fdx2 = d2fdx2(model, state, ref)
            d2fdx2 = zeros(model.numStates, model.numStates, model.numStates);

            d2fdx2(4, :, :) = derivatives_drone_IB_syms_d2fdx2.d2f4dx2(model, state, ref);
            d2fdx2(5, :, :) = derivatives_drone_IB_syms_d2fdx2.d2f5dx2(model, state, ref);
            d2fdx2(6, :, :) = derivatives_drone_IB_syms_d2fdx2.d2f6dx2(model, state, ref);

            % d2fdx2(10, :, :) = derivatives_drone_IB_syms_d2fdx2_10.d2f10dx2(model, state, ref);
            d2fdx2(11, :, :) = derivatives_drone_IB_syms_d2fdx2_11.d2f11dx2(model, state, ref);
            % d2fdx2(20, :, :) = derivatives_drone_IB_syms.d2f20dx2(model, state, ref);
            % d2fdx2(23, :, :) = derivatives_drone_IB_syms.d2f23dx2(model, state, ref);
            % d2fdx2(24, :, :) = derivatives_drone_IB_syms.d2f24dx2(model, state, ref);
        end

        function d2f4dx2 = d2f4dx2(model, state, ~)
            d2f4dx2 = zeros(model.numStates, model.numStates);

            
            x7 = state(7);
            x8 = state(8);
            x9 = state(9);
            x19 = state(19);
            
            
            m = model.m;

            d2f4dx2(7,7)=-(4*x19*cos(x8)*cos(x9)*sin(x7))/m;
           d2f4dx2(7,8)=-(4*x19*cos(x7)*cos(x9)*sin(x8))/m;
           d2f4dx2(7,9)=-(4*x19*cos(x7)*cos(x8)*sin(x9))/m;
           d2f4dx2(7,19)=(4*cos(x7)*cos(x8)*cos(x9))/m;
           d2f4dx2(8,7)=-(4*x19*cos(x7)*cos(x9)*sin(x8))/m;
           d2f4dx2(8,8)=-(4*x19*(sin(x8)*sin(x9) + cos(x8)*cos(x9)*sin(x7)))/m;
           d2f4dx2(8,9)=(4*x19*(cos(x8)*cos(x9) + sin(x7)*sin(x8)*sin(x9)))/m;
           d2f4dx2(8,19)=(4*(cos(x8)*sin(x9) - cos(x9)*sin(x7)*sin(x8)))/m;
           d2f4dx2(9,7)=-(4*x19*cos(x7)*cos(x8)*sin(x9))/m;
           d2f4dx2(9,8)=(4*x19*(cos(x8)*cos(x9) + sin(x7)*sin(x8)*sin(x9)))/m;
           d2f4dx2(9,9)=-(4*x19*(sin(x8)*sin(x9) + cos(x8)*cos(x9)*sin(x7)))/m;
           d2f4dx2(9,19)=(4*(cos(x9)*sin(x8) - cos(x8)*sin(x7)*sin(x9)))/m;
           d2f4dx2(19,7)=(4*cos(x7)*cos(x8)*cos(x9))/m;
           d2f4dx2(19,8)=(4*cos(x8)*sin(x9) - 4*cos(x9)*sin(x7)*sin(x8))/m;
           d2f4dx2(19,9)=(4*cos(x9)*sin(x8) - 4*cos(x8)*sin(x7)*sin(x9))/m;
        end
                
        function d2f5dx2 = d2f5dx2(model, state, ~)
            d2f5dx2 = zeros(model.numStates, model.numStates);

           
            x7 = state(7);
            x8 = state(8);
            x9 = state(9);
            x19 = state(19);
            
            c7 = model.c7; 
            c8 = model.c8;
            lambda4 = model.lambda4;
            m = model.m;

            d2f5dx2(3,7)=(sin(x7)*(cos(x9)*sin(x8) - cos(x8)*sin(x7)*sin(x9))*(lambda4 - c7^2 + c7*(c7 + c8) + 1))/(cos(x7)^2*cos(x8)) - sin(x9)*(lambda4 - c7^2 + c7*(c7 + c8) + 1);
            d2f5dx2(7,7)=-(4*x19*cos(x8)*sin(x7)*sin(x9))/m;
           d2f5dx2(7,8)=-(4*x19*cos(x7)*sin(x8)*sin(x9))/m;
           d2f5dx2(7,9)=(4*x19*cos(x7)*cos(x8)*cos(x9))/m;
           d2f5dx2(7,19)=(4*cos(x7)*cos(x8)*sin(x9))/m;
           d2f5dx2(8,7)=-(4*x19*cos(x7)*sin(x8)*sin(x9))/m;
           d2f5dx2(8,8)=(4*x19*(cos(x9)*sin(x8) - cos(x8)*sin(x7)*sin(x9)))/m;
           d2f5dx2(8,9)=(4*x19*(cos(x8)*sin(x9) - cos(x9)*sin(x7)*sin(x8)))/m;
           d2f5dx2(8,19)=-(4*(cos(x8)*cos(x9) + sin(x7)*sin(x8)*sin(x9)))/m;
           d2f5dx2(9,7)=(4*x19*cos(x7)*cos(x8)*cos(x9))/m;
           d2f5dx2(9,8)=(4*x19*(cos(x8)*sin(x9) - cos(x9)*sin(x7)*sin(x8)))/m;
           d2f5dx2(9,9)=(4*x19*(cos(x9)*sin(x8) - cos(x8)*sin(x7)*sin(x9)))/m;
           d2f5dx2(9,19)=(4*(sin(x8)*sin(x9) + cos(x8)*cos(x9)*sin(x7)))/m;
           d2f5dx2(19,7)=(4*cos(x7)*cos(x8)*sin(x9))/m;
           d2f5dx2(19,8)=-(4*cos(x8)*cos(x9) + 4*sin(x7)*sin(x8)*sin(x9))/m;
           d2f5dx2(19,9)=(4*sin(x8)*sin(x9) + 4*cos(x8)*cos(x9)*sin(x7))/m;
        end

        function d2f6dx2 = d2f6dx2(model, state, ~)
            d2f6dx2 = zeros(model.numStates, model.numStates);

           
            x7 = state(7);
            x8 = state(8);
            x19 = state(19);

            m = model.m;

            d2f6dx2(7,7)=-(4*x19*cos(x7)*cos(x8))/m;
           d2f6dx2(7,8)=(4*x19*sin(x7)*sin(x8))/m;
           d2f6dx2(7,19)=-(4*cos(x8)*sin(x7))/m;
           d2f6dx2(8,7)=(4*x19*sin(x7)*sin(x8))/m;
           d2f6dx2(8,8)=-(4*x19*cos(x7)*cos(x8))/m;
           d2f6dx2(8,19)=-(4*cos(x7)*sin(x8))/m;
           d2f6dx2(19,7)=-(4*cos(x8)*sin(x7))/m;
           d2f6dx2(19,8)=-(4*cos(x7)*sin(x8))/m;
        end
    end
end