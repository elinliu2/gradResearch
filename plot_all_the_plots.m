function plot_all_the_plots(SC)
    total = int32((SC.stop-SC.start)/SC.step)+1;
    
    % initialize arrays
    x1_val =        zeros([1 total*total]);
    x2_val =        zeros([1 total*total]);

    z =             zeros([1 total*total]);
    ztest =         zeros([1 total*total]);
    z_diff =        zeros([1 total*total]);
    % 
    % zparam =        zeros([1 total*total]);
    % zparamtest =    zeros([1 total*total]);
    % zparam_diff =   zeros([1 total*total]);
    % 
    z2 =            zeros([1 total*total]);
    z2test =        zeros([1 total*total]);
    z2_diff =       zeros([1 total*total]);
    % 
    % z2param =       zeros([1 total*total]);
    % z2paramtest =   zeros([1 total*total]);
    % z2param_diff =  zeros([1 total*total]);
    zdG =           zeros([1 total*total]);
    U =             zeros([1 total*total]);
    L =             zeros([1 total*total]);

    % classes so I can call functions
    tjs = trajectorySensitivities;
    norms = normalStuff;
    splots = splotting;    
    sCtrl = stonkCtrllers;

    % calculate trajectory sensitivities for different initial conditions
    index = 1; 
    for i = SC.start:SC.step:SC.stop
        for j = SC.start:SC.step:SC.stop
            x1_0 = i;
            x2_0 = j;

            % % calculate first order derivatives wrt x0
            [stateProgression, ~, ~] = pendulum_trajectory(SC, x1_0, x2_0, @sCtrl.dummyU, @sCtrl.dummyD, false);
            trajSens = tjs.calculateTrajectorySensitivities(stateProgression, SC);
            trajSensTest = tjs.calculateTrajectorySensitivitiesTest(i, j, SC);
    
            % % calculate second derivative wrt x0
            secondTrajSens = tjs.calculateSecondOrderTrajectorySensitivities(stateProgression, trajSens, SC);
            secondTrajSensTest = tjs.calculateSecondOrderTrajectorySensitivitiesTest(i, j, SC);
            
            G = norms.infimumOfReciprocalOfOneNorm(trajSens);
            
            % % calculate first order trajectory sensitivities wrt p
            % trajSensParam = tjs.calculateTrajectorySensitivitiesParameters(stateProgression, SC);
            % % trajSensTestParam = deriv.calculateTrajectorySensitivitiesTestParams(i, j, 0, 0, SC);
            % 
            % % % calculate second derivative wrt p
            % secondTrajSensParam = tjs.calculateSecondOrderTrajectorySensitivitiesParams(stateProgression, trajSens, trajSensParam, SC);
            % % secondTrajSensTestParam = deriv.calculateSecondOrderTrajectorySensitivitiesTestParams(i, j, SC);
            % 
            % tp = tjs.findtp(trajSens);
            % 
            % dGp = tjs.calculate_dG(G, trajSens, secondTrajSensParam, tp);

            z(index) = G;
            ztest(index) = norms.infimumOfReciprocalOfOneNorm(trajSensTest);
            z_diff(index) = norms.infimum(1./trajSens);

            z2(index) = norms.secondOrderMax(secondTrajSens);
            z2test(index) = norms.secondOrderMax(secondTrajSensTest);
            z2_diff(index) = norms.secondOrderMax(secondTrajSens-secondTrajSensTest);
            
            % zparam(index) = norms.infimumOfReciprocalOfOneNorm(trajSensParam);
            % zparamtest(index) = norms.infimumOfReciprocalOfOneNorm(trajSensTestParam);
            % zparam_diff(index) = norms.supremum(trajSensParam-trajSensTestParam);
            % 
            % z2param(index) = norms.secondOrderMax(secondTrajSensParam);
            % z2paramtest(index) = norms.secondOrderMax(secondTrajSensTestParam);
            % z2param_diff(index) = norms.secondOrderMax(secondTrajSensParam-secondTrajSensTestParam);

            % zdG(index) = norm(dGp);
            % U(index) = dGp(1);
            % L(index) = dGp(2);
            
            x1_val(index) = i;
            x2_val(index) = j;
            index = index+1;
        end 
    end 
    
    index = index-1;
    splots.splot(1, x1_val(:, 1:index), x2_val(:, 1:index), z(:, 1:index), "infimum inverse trajectory sensitivity");
    splots.splot(2, x1_val(:, 1:index), x2_val(:, 1:index), ztest(:, 1:index), "infimum inverse trajectory sensitivity test");
    splots.splot(3, x1_val(:, 1:index), x2_val(:, 1:index), z_diff(:, 1:index), "wrong infimum");
    disp("max z diff: " + max(z_diff));

    splots.splot(4, x1_val, x2_val, z2, "second order trajectory sensitivity");
    splots.splot(5, x1_val, x2_val, z2test, "second order trajectory sensitivity test");
    splots.splot(6, x1_val, x2_val, z2_diff, "max diff second order trajectory sensitivity");
    disp("max z2 diff: " + max(z2_diff));
    % 
    % splots.splot(7, x1_val(:, 1:index), x2_val(:, 1:index), zparam(:, 1:index), "infimum inverse trajectory sensitivity param");
    % splots.splot(8, x1_val(:, 1:index), x2_val(:, 1:index), zparamtest(:, 1:index), "infimum inverse trajectory sensitivity test param");
    % splots.splot(9, x1_val(:, 1:index), x2_val(:, 1:index), zparam_diff(:, 1:index), "supremum diff trajectory sensitivity param");
    % 
    % splots.splot(10, x1_val(:, 1:index), x2_val(:, 1:index), z2param(:, 1:index), "second order trajectory sensitivity param");
    % zlim([0, 600]);
    % splots.splot(11, x1_val(:, 1:index), x2_val(:, 1:index), z2paramtest(:, 1:index), "second order trajectory sensitivity test param");
    % zlim([0, 600]);
    % splots.splot(12, x1_val(:, 1:index), x2_val(:, 1:index), z2param_diff(:, 1:index), "supremum diff trajectory sensitivity param");
    % zlim([0, 30]);

    % splots.splot(13, x1_val(:, 1:index), x2_val(:, 1:index), zdG(:, 1:index), "dG");
    % splots.squiver(14, x1_val(:, 1:index), x2_val(:, 1:index), U(:, 1:index), L(:, 1:index), "u l dG");


end