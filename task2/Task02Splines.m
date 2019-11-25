% Task 01
clear all; clc;
syms a01 a11 a21 a31 a41 a02 a12 a22 a32 a42 t 
%Initial time
ti = 0; tm = 2;  tf = 4;
%Desired position and velocities
q1 = 1; dq1 = 0; ddq1 = 0;
q2 = 2; 
q3 = 0; dq3 = 0; ddq3 = 0;

q_des1 = [q1 q2 dq1 ddq1]';

%%  Polynomial1s

% Polinomial expressions
q11 = a01 + a11*t + a21*t^2 + a31*t^3 + a41*t^4;
dq11 = diff(q11, t);
%Polynomial matrices
P1 = [1, t, t^2, t^3];
dP1 = diff(P1);
ddP1 = diff(dP1);


%Polynomial 1
p11 = subs(P1, {t}, {ti});
p21 = subs(P1, {t}, {tm});
dp11 = subs(dP1, {t}, {ti});
ddp11 = subs(ddP1, {t}, {ti});


A1 = [p11; p21; dp11; ddp11; ];
X1 = A1\q_des1;
qr1 = fliplr(double(X1)');
dqr1 = polyder(qr1);
ddqr1 = polyder(dqr1);

%%
dq2 = polyval(dqr1, tm); ddq2 = polyval(ddqr1, tm);
q_des2 = [q2 q3 dq2 dq3 ddq2 ddq3]';

%Polynomial 2 
P2 = [1, t, t^2, t^3, t^4, t^5];
dP2 = diff(P2);
ddP2 = diff(dP2);

%Polynomial 2 
p12 = subs(P2, {t}, {tm});
p22 = subs(P2, {t}, {tf});
dp12 = subs(dP2, {t}, {tm});
dp22 = subs(dP2, {t}, {tf});
ddp12 = subs(ddP2, {t}, {tm});
ddp22 = subs(ddP2, {t}, {tf});

%x = [a0, a1, a2, a3 ]';
A2 = [p12; p22; dp12; dp22;  ddp12; ddp22];
X2 = A2\q_des2;
qr2 = fliplr(double(X2)');
dqr2 = polyder(qr2);
ddqr2 = polyder(dqr2);


st=0.001;
time1 = ti : st : tm;
time2 = tm+st : st : tf;
time = [time1, time2];

qr = [polyval(qr1, time1), polyval(qr2, time2)];
dqr = [polyval(dqr1, time1), polyval(dqr2, time2)];
ddqr = [polyval(ddqr1, time1), polyval(ddqr2, time2)];

%Plot
figure(1);
plot(qr,'b',  'LineWidth', 2)
xlabel('Time of moving')
ylabel('Trajectory')
legend('q_{pol}', 'Fontsize', 10)
title('Trajectory splines')
grid on


figure(2);
plot(dqr,'g', 'LineWidth', 2)
xlabel('Time of moving')
ylabel('Velocity ')
legend('dq_{pol}', 'Fontsize', 10)
title('Velocity splines')
grid on

figure(3);
plot(ddqr, 'r', 'LineWidth', 2)
xlabel('Time of moving')
ylabel('Acceleration ')
legend('ddq_{pol}', 'Fontsize', 10)
title('Acceleration splins')
grid on
