function plot_all_the_plots(SC, model)
    total = int32((SC.stop-SC.start)/SC.step)+1;
    
    % initialize arrays
    x1_val =        zeros([1 total*total]);
    x2_val =        zeros([1 total*total]);

    z =             zeros([1 total*total]);
    ztest =         zeros([1 total*total]);
    z_diff =        zeros([1 total*total]);
    
    zparam =        zeros([1 total*total]);
    zparamtest =    zeros([1 total*total]);
    zparam_diff =   zeros([1 total*total]);

    z2param =       zeros([1 total*total]);
    z2paramtest =   zeros([1 total*total]);
    z2param_diff =  zeros([1 total*total]);

    zdG =           zeros([1 total*total]);
    zdGTest =       zeros([1 total*total]);
    zdG_diff =      zeros([1 total*total]);
    U =             zeros([1 total*total]);
    L =             zeros([1 total*total]);

    Utest =             zeros([1 total*total]);
    Ltest =             zeros([1 total*total]);

    RoA =           zeros([1 total*total]);

    % classes so I can call functions
    tjs = trajectorySensitivities;
    splots = splotting;    

    % calculate trajectory sensitivities for different initial conditions
    index = 1; 
    for i = SC.start:SC.step:SC.stop
        for j = SC.start:SC.step:SC.stop

            initialState = model.initializeStates(model, i, j);

            % calculate first order derivatives wrt x0
            [~, stateProgression, ~] = trajectory(model, SC, initialState, false);
            trajSens = tjs.calculateTrajectorySensitivities(model, SC, stateProgression);
            trajSensTest = tjs.calculateTrajectorySensitivitiesTest(model, SC, initialState, length(stateProgression));
            G = norms.infimumOfReciprocalOfElementOneNorm(trajSens);
        
            % % calculate first order trajectory sensitivities wrt p
            % trajSensParam = tjs.calculateTrajectorySensitivitiesParameters(model, SC, stateProgression);
            % trajSensTestParam = tjs.calculateTrajectorySensitivitiesTestParams(model, SC, initialState, length(stateProgression));
            % 
            % % calculate second derivative wrt p
            % secondTrajSensParam = tjs.calculateSecondOrderTrajectorySensitivitiesParams(model, SC, stateProgression, trajSens, trajSensParam);
            % secondTrajSensTestParam = tjs.calculateSecondOrderTrajectorySensitivitiesTestParams(model, SC, initialState, length(stateProgression));

            % tp = tjs.findtp(trajSens);
            % 
            % dGp = tjs.calculate_dG(G, trajSens, secondTrajSensParam, tp, model);
            % dGpTest = tjs.calculate_dGTest(model, SC, initialState);

            z(index) = G;
            ztest(index) = norms.infimumOfReciprocalOfElementOneNorm(trajSensTest);
            z_diff(index) = norms.supremum(trajSens-trajSensTest);

            % zparam(index) = norms.supremum(trajSensParam);
            % zparamtest(index) = norms.supremum(trajSensTestParam);
            % zparam_diff(index) = norms.supremum(trajSensParam-trajSensTestParam);
            % 
            % z2param(index) = norms.secondOrderMax(secondTrajSensParam);
            % z2paramtest(index) = norms.secondOrderMax(secondTrajSensTestParam);
            % z2param_diff(index) = norms.secondOrderMax(secondTrajSensParam-secondTrajSensTestParam);

            % zdG(index) = norm(dGp);
            % zdGTest(index) = norm(dGpTest);
            % zdG_diff(index) = norm(dGp-dGpTest);
            % U(index) = dGp(1)/zdG(index);
            % L(index) = dGp(2)/zdG(index);
            % 
            % Utest(index) = dGpTest(1)/zdGTest(index);
            % Ltest(index) = dGpTest(2)/zdGTest(index);
            
            x1_val(index) = i;
            x2_val(index) = j;

            % if (abs(stateProgression(length(stateProgression), 1)) < 0.1) 
            %     RoA(index) = 1;
            % else 
            %     RoA(index) = -1;
            % end 
            index = index+1;
        end 
    end 
    
    index = index-1;

    splots.splot(1, x1_val(:, 1:index), x2_val(:, 1:index), z(:, 1:index), "infimum inverse trajectory sensitivity");
    splots.splot(2, x1_val(:, 1:index), x2_val(:, 1:index), ztest(:, 1:index), "infimum inverse trajectory sensitivity test");
    splots.splot(3, x1_val(:, 1:index), x2_val(:, 1:index), z_diff(:, 1:index), "max diff");
    disp("max z diff: " + max(z_diff));
    % 
    % splots.splot(7, x1_val(:, 1:index), x2_val(:, 1:index), zparam(:, 1:index), "supremum trajectory sensitivity param");
    % splots.splot(8, x1_val(:, 1:index), x2_val(:, 1:index), zparamtest(:, 1:index), "supremum trajectory sensitivity test param");
    % splots.splot(9, x1_val(:, 1:index), x2_val(:, 1:index), zparam_diff(:, 1:index), "supremum diff trajectory sensitivity param");
    % disp("max zparam diff: " + max(zparam_diff));
    % 
    % splots.splot(10, x1_val(:, 1:index), x2_val(:, 1:index), z2param(:, 1:index), "second order trajectory sensitivity param");
    % splots.splot(11, x1_val(:, 1:index), x2_val(:, 1:index), z2paramtest(:, 1:index), "second order trajectory sensitivity test param");
    % splots.splot(12, x1_val(:, 1:index), x2_val(:, 1:index), z2param_diff(:, 1:index), "supremum diff trajectory sensitivity param");
    % disp("max z2param diff: " + max(z2param_diff));

    % splots.splot(13, x1_val(:, 1:index), x2_val(:, 1:index), zdG(:, 1:index), "dG");
    % splots.splot(14, x1_val(:, 1:index), x2_val(:, 1:index), zdGTest(:, 1:index), "dGtest");
    % splots.splot(15, x1_val(:, 1:index), x2_val(:, 1:index), zdG_diff(1:index), "dG diff");
    % disp("max zdG_diff diff: " + max(zdG_diff));
    % 
    % splots.squiver(16, x1_val(:, 1:index), x2_val(:, 1:index), U(:, 1:index), L(:, 1:index), "dG direction");
    % splots.squiver(17, x1_val(:, 1:index), x2_val(:, 1:index), Utest(:, 1:index), Ltest(:, 1:index), "dG direction test");
    % 
    % splots.splot(18, x1_val(:, 1:index), x2_val(:, 1:index), RoA(:, 1:index), "RoA");
end