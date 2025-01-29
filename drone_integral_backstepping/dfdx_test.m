model.g = 9.8;
model.m = 1.4; % kg

model.l = 0.2; % m
model.K_psi = 4; % Nm

model.Ja = 0.03; % kgm^2
model.Jb = 0.03; % kgm^2
model.Jc = 0.04; % kgm^2

model.kp1 = 5;
model.ki1 = 0.5;
model.kd1 = 10;

model.kp2 = 5;
model.ki2 = 1;
model.kd2 = 5;

model.kp3 = 5;
model.ki3 = 0.5;
model.kd3 = 10;

model.kp4 = 20;
model.ki4 = 3;
model.kd4 = 20;

initialState = zeros(16, 1);
initialState(3) = 100; % Drone starts 100 m in the air
initialState(8) = pi/6;

dfdx = zeros(16, 16);

delta = 1e-5;
for i = 1:16
    deltaStatePlus = initialState;
    deltaStatePlus(i) = deltaStatePlus(i) + delta;
    deltaStateMinus = initialState;
    deltaStateMinus(i) = deltaStateMinus(i) - delta;
    dfdx(:,i) = transpose(dynamics(deltaStatePlus, model)-dynamics(deltaStateMinus, model))/(2*delta);
end


function f = dynamics(state, model)
    t=0;
    pdot_x = state(4);
    pdot_y = state(5);
    pdot_z = state(6);
    
    theta = state(7); % theta (pitch)
    phi = state(8); % phi (roll)
    psi = state(9); % psi (yaw)
    
    dot_a = state(10);
    dot_b = state(11);
    dot_c = state(12);
    
    dx1 = pdot_x; 
    dx2 = pdot_y;
    dx3 = pdot_z;
    
    theta_ref = pi/6*sin(t);
    phi_ref = pi/6*cos(t);
    psi_ref = t;
    
    theta_ref_dot = pi/6*cos(t);
    phi_ref_dot = -pi/6*sin(t);
    psi_ref_dot = 1;
    
    u_theta = model.kp1 * (theta_ref - state(7)) + model.ki1 * state(13) + model.kd1 * ( theta_ref_dot - state(10));
    u_phi = model.kp2 * (phi_ref - state(8)) + model.ki2 * state(14) + model.kd2 * ( phi_ref_dot - state(11));
    u_psi = model.kp3 * (psi_ref - state(9)) + model.ki3 * state(15) + model.kd3 * (psi_ref_dot - state(12));
    
    u_thrust = 3.961;
    
    F1 = u_theta + u_psi + u_thrust;
    F2 = - u_theta + u_psi + u_thrust;
    F3 = u_phi - u_psi + u_thrust;
    F4 = - u_phi - u_psi + u_thrust;
    
    dx4 = 1/(model.m) * (F1+F2+F3+F4) * (sin(theta)*cos(phi)*cos(psi) + sin(phi)*sin(psi));
    dx5 = 1/(model.m) * (F1+F2+F3+F4) * (sin(theta)*cos(phi)*sin(psi) - sin(phi)*cos(psi));
    dx6 = 1/(model.m) * (F1+F2+F3+F4) * (cos(theta)*cos(phi)) - model.g;
    
    dx7 = dot_a;
    dx8 = dot_b;
    dx9 = dot_c;
    
    dx10 = 1/(model.Ja) * model.l * (F1 - F2); 
    dx11 = 1/(model.Jb) * model.l * (F3 - F4);
    dx12 = 1/(model.Jc) * model.K_psi * (F1 + F2 - F3 - F4);
    
    ref(3) = 101;
    delta_z = ref(3) - state(3);
    
    % theta_ref = atan2(delta_x, 1);
    % phi_ref = atan2(-delta_y, 1);
    % psi_ref = 0; 
    
    theta_ref = pi/6*sin(t);
    phi_ref = pi/6*cos(t);
    psi_ref = t;
    
    dx13 = theta_ref - state(7);
    dx14 = phi_ref   - state(8);
    dx15 = psi_ref   - state(9);
    % dx16 = norm(ref-state(1:3));
    dx16 = delta_z;
    
    f = [dx1; dx2; dx3; dx4; dx5; dx6; dx7; dx8; dx9; dx10; dx11; dx12; dx13; dx14; dx15; dx16];

end






