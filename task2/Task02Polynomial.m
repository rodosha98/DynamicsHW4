% Task 01
clear all; clc;
syms a0 a1 a2 a3 a4 a5 a6  t q
%Initial time
t1 = 0; t2 = 2; t3 = 4;
%Desired position and velocities
q1 = 1; dq1 = 0; ddq1 =0;
q2 = 2; 
q3 = 0; dq3 = 0; ddq3 = 0;
q_des = [q1 q2 q3 dq1 dq3 ddq1 ddq3]';
%%  Polynomial

% Polinomial expressions
q = a0 + a1*t + a2*t^2 + a3*t^3 + a4*t^4 + a5*t^5 + a6*t^6;
dq = diff(q, t);
%Polynomial matrices
P = [1, t, t^2, t^3, t^4, t^5, t^6];
dP = diff(P);
ddP = diff(dP);

%Polynomial matrices with initial conditions
p1 = subs(P, {t}, {t1});
p2 = subs(P, {t}, {t2});
p3 = subs(P, {t}, {t3});

dp1 = subs(dP, {t}, {t1});
dp2 = subs(dP, {t}, {t3});

ddp1 = subs(ddP, {t}, {t1});
ddp2 = subs(ddP, {t}, {t3});


%x = [a0, a1, a2, a3 ]';

A = [p1; p2; p3; dp1; dp2; ddp1; ddp2];

X = A\q_des;
qr = fliplr(double(X)');
dqr = polyder(qr);
ddqr = polyder(dqr);


st=0.01;
time = t1 : st : t3;
%Plot
figure(1);
plot(polyval(qr, time),'b',  'LineWidth', 2)
xlabel('Time of moving')
ylabel('Trajectory')
legend('q_{pol}', 'Fontsize', 10)
title('Trajectory polynomial')
grid on

figure(2);
plot(polyval(dqr, time), 'g', 'LineWidth', 2)
xlabel('Time of moving')
ylabel('Velocity ')
legend('dq_{pol}', 'Fontsize', 10)
title('Velocity polynomial')
grid on

figure(3);
plot(polyval(ddqr, time), 'r', 'LineWidth', 2)
xlabel('Time of moving')
ylabel('Acceleration ')
legend('ddq_{pol}', 'Fontsize', 10)
title('Acceleration polynomial')
grid on





