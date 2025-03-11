classdef derivatives_drone_IB_syms_d2fdx2_20
    methods (Static)
        function d2f20dx2 = d2f20dx2(model, state, ref)
            d2f20dx2 = zeros(model.numStates, model.numStates);
            
            d2f20dx2(1,:)=derivatives_drone_IB_syms_d2fdx2_20_1.d2f20dx2_1(model, state, ref);
            d2f20dx2(2,:)=derivatives_drone_IB_syms_d2fdx2_20_2.d2f20dx2_2(model, state, ref);
            d2f20dx2(3,:)=derivatives_drone_IB_syms_d2fdx2_20_3.d2f20dx2_3(model, state, ref);
            d2f20dx2(4,:)=derivatives_drone_IB_syms_d2fdx2_20_4.d2f20dx2_4(model, state, ref);
            d2f20dx2(5,:)=derivatives_drone_IB_syms_d2fdx2_20_5.d2f20dx2_5(model, state, ref);
            d2f20dx2(6,:)=derivatives_drone_IB_syms_d2fdx2_20_6.d2f20dx2_6(model, state, ref);
            d2f20dx2(7,:)=derivatives_drone_IB_syms_d2fdx2_20_7.d2f20dx2_7(model, state, ref);
            d2f20dx2(8,:)=derivatives_drone_IB_syms_d2fdx2_20_8.d2f20dx2_8(model, state, ref);
            d2f20dx2(9,:)=derivatives_drone_IB_syms_d2fdx2_20_9.d2f20dx2_9(model, state, ref);
            d2f20dx2(10,:)=derivatives_drone_IB_syms_d2fdx2_20_10.d2f20dx2_10(model, state, ref);
            d2f20dx2(11,:)=derivatives_drone_IB_syms_d2fdx2_20_11.d2f20dx2_11(model, state, ref);
            d2f20dx2(12,:)=derivatives_drone_IB_syms_d2fdx2_20_12.d2f20dx2_12(model, state, ref);
            d2f20dx2(13,:)=derivatives_drone_IB_syms_d2fdx2_20_13.d2f20dx2_13(model, state, ref);
            d2f20dx2(14,:)=derivatives_drone_IB_syms_d2fdx2_20_14.d2f20dx2_14(model, state, ref);
            d2f20dx2(16,:)=derivatives_drone_IB_syms_d2fdx2_20_16.d2f20dx2_16(model, state, ref);
            d2f20dx2(17,:)=derivatives_drone_IB_syms_d2fdx2_20_17.d2f20dx2_17(model, state, ref);
            d2f20dx2(18,:)=derivatives_drone_IB_syms_d2fdx2_20_18.d2f20dx2_18(model, state, ref);
            d2f20dx2(19,:)=derivatives_drone_IB_syms_d2fdx2_20_19.d2f20dx2_19(model, state, ref);
            d2f20dx2(20,:)=derivatives_drone_IB_syms_d2fdx2_20_20.d2f20dx2_20(model, state, ref);
            d2f20dx2(21,:)=derivatives_drone_IB_syms_d2fdx2_20_21.d2f20dx2_21(model, state, ref);
            d2f20dx2(22,:)=derivatives_drone_IB_syms_d2fdx2_20_22.d2f20dx2_22(model, state, ref);
            d2f20dx2(23,:)=derivatives_drone_IB_syms_d2fdx2_20_23.d2f20dx2_23(model, state, ref);
            d2f20dx2(24,:)=derivatives_drone_IB_syms_d2fdx2_20_24.d2f20dx2_24(model, state, ref);           
        end
    end
end