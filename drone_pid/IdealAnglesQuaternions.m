time = 0:1e-3:9*pi;
quat = normalize(quaternion([0 pi/6 0], 'euler', 'ZYX', 'frame'));
theta_ref = zeros(length(time) , 1);
phi_ref = zeros(length(time) , 1);
psi_ref = zeros(length(time) , 1);
for i = 1: length(time)
    quat_delta = normalize(quaternion([0 0 time(i)], 'euler', 'ZYX', 'frame'));
    new_quat = quat_delta*quat;
    eul = quat2eul(new_quat);

    theta_ref(i) = eul(1);
    phi_ref(i) = eul(2);
    psi_ref(i) = eul(3);
end

figure(3)

subplot(3, 1, 1);
plot(time,theta_ref);
subplot(3, 1, 2);
plot(time,phi_ref);
subplot(3, 1, 3);
plot(time,psi_ref);