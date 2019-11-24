% Task 01
clear all; clc;
syms a0 a1 a2 a3 t q
%Initial time
ti = 0; tf = 2;
%Desired position and velocities
qi = 1; dqi = 0;
qf = 4; dqf = 0;
q_des = [qi qf dqi dqf]';
%%  Polynomial

% Polinomial expressions
q = a0 + a1*t + a2*t^2 + a3*t^3;
dq = diff(q, t);
%Polynomial matrices
P = [1, t, t^2, t^3];
dP = diff(P);

%Polynomial matrices with initial conditions
p1 = subs(P, {t}, {ti});
p2 = subs(P, {t}, {tf});
dp1 = subs(dP, {t}, {ti});
dp2 = subs(dP, {t}, {tf});
%x = [a0, a1, a2, a3 ]';

A = [p1; p2; dp1; dp2];

X = A\q_des;
qr = fliplr(double(X)');
dqr = polyder(qr);
ddqr = polyder(dqr);

st=0.01;
time = ti : st : tf;
%Plot
figure(1);
plot(polyval(qr, time),'b',  'LineWidth', 2)
xlabel('Time of moving')
ylabel('Trajectory')
legend('q_{pol}', 'Fontsize', 10)
title('Trajectory polynomial')
grid on

figure(2);
plot(polyval(dqr, time),'g', 'LineWidth', 2)
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

%% Trapezoidal
%Preparation
% flag = 1 to use velocity method, flag = 0 for acceleration 
flag = 1;

if flag ==1
    %assign desired velocity
    dqc = 1.6;
    %constraint for desired acceleration
    v_cons = abs(qf - qi)/ (tf-ti);
    
    if (dqc < v_cons) || (dqc > 2*v_cons)
        disp('Error. Reassign dqc. It should be in limits'); 
        disp(v_cons);
        disp(2*v_cons);
    end
    % time c
    tc = (qi -qf +dqc*tf)/dqc;
    ddqc = dqc^2/(qi -qf +dqc*tf);
else
    %Constrain for accelearation in point c
    ddqc_star = 4*(qf - qi)/tf^2;
    %Assign desired acceleration
    ddqc = 4;
    if ddqc < ddqc_star
        disp('Error. Reassign ddqc. It should be greater than'); 
        disp(ddqc_star);
    end
    % time tc
    tc = tf/2 - 1/2*sqrt((tf^2*ddqc-4*(qf-qi))/ddqc);

end
tc = round(tc,3);

%polinomials for trapezoid
syms t2
tr1 = qi + 1/2*ddqc*t2^2;
tr2 = qi + ddqc*tc*(t2-tc/2);
tr3 = qf - 1/2*ddqc*(tf-t2)^2;

dtr1 = diff(tr1, t2);
dtr2 = diff(tr2, t2);
dtr3 = diff(tr3, t2);

ddtr1 = diff(dtr1, t2);
ddtr2 = diff(dtr2, t2);
ddtr3 = diff(dtr3, t2);

st=0.001;
time = ti : st : tf;
time1 = 0 : st : tc;
time2 = (tc+st) : st : (tf - tc);
time3 = (tf-tc+st) : st : (tf);

tr1_t = double(subs(tr1, {t2}, {time1}));
tr2_t = double(subs(tr2, {t2}, {time2}));
tr3_t = double(subs(tr3, {t2}, {time3}));

dtr1_t = double(subs(dtr1, {t2}, {time1}));
dtr2_t = double(subs(dtr2, {t2}, {time2}));
dtr3_t = double(subs(dtr3, {t2}, {time3}));

ddtr1_t = double(subs(ddtr1, {t2}, {time1}));
ddtr2_t = double(subs(ddtr2, {t2}, {time2}));
ddtr3_t = double(subs(ddtr3, {t2}, {time3}));

tr = [tr1_t, tr2_t, tr3_t];
dtr = [dtr1_t, dtr2_t, dtr3_t];
ddtr = [ddtr1_t, ddtr2_t, ddtr3_t];



figure(4)
plot(time, tr,'b',  'LineWidth', 2)
xlabel('Time of moving')
ylabel('Trajectory')
legend('q_{tr}', 'Fontsize', 10)
title('Trajectory trapezoidal')
grid on

figure(5)
plot(time, dtr, 'g',  'LineWidth', 2)
xlabel('Time of moving')
ylabel('Velocity')
legend('dq_{tr}', 'Fontsize', 10)
title('Velocity trapezoidal')
grid on

figure(6)
plot(time, ddtr, 'r',  'LineWidth', 2)
xlabel('Time of moving')
ylabel('Acceleration')
legend('ddq_{tr}', 'Fontsize', 10)
title('Acceleration trapezoidal')
grid on

%%
syms t3 a b c 
% assign frequency
k = 1;

as = pi*(2*k+1);
bs = 0;
cs = pi*(2*k+1)/2;

q3 = a*b*t3 - a/c*(cos(c*t3)-1);
dq3 = a*(b+sin(c*t3));
ddq3 = c*cos(c*t3);

st=0.001;
time = ti : st : tf;


q3_t = double(subs(q3, {t3, a, b, c}, {time, as, bs, cs}));
dq3_t = double(subs(dq3, {t3, a, b, c}, {time, as, bs, cs}));
ddq3_t = double(subs(ddq3, {t3, a, b, c}, {time, as, bs, cs}));


figure(7)
plot(time, q3_t,'b',  'LineWidth', 2)
xlabel('Time of moving')
ylabel('Trajectory')
legend('q_{sin}', 'Fontsize', 10)
title('Trajectory sin')
grid on

figure(8)
plot(time, dq3_t, 'g',  'LineWidth', 2)
xlabel('Time of moving')
ylabel('Velocity')
legend('dq_{sin}', 'Fontsize', 10)
title('Velocity sin')
grid on

figure(9)
plot(time, ddq3_t, 'r',  'LineWidth', 2)
xlabel('Time of moving')
ylabel('Acceleration')
legend('ddq_{sin}', 'Fontsize', 10)
title('Acceleration sin')
grid on





