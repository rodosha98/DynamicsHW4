% Task 02 Trapezoidal
clear all; clc;
syms a0 a1 a2 a3 a4 a5 a6  t q
%Initial time
ti = 0; tm = 2; tf = 4;
%Desired position and velocities
q1 = 1; dq1 = 0; ddq1 =0;
q2 = 2; 
q3 = 0; dq3 = 0; ddq3 = 0;
q_des = [q1 q3 q2 dq1 dq3 ddq1 ddq3]';

%% Trapezoidal

% flag = 1 to use velocity method, flag = 0 for acceleration 
flag = 0;

for i = 1:2
    if i==1
        qi = q1;
        qf = q2;
    else
        qi = q2;
        qf = q3;
    end
    %Preparation
    
    
    if flag == 0
        %assign desired velocity
        dqc1 = 0.8;
        dqc2 = 1.5;
        if i==1
            dqc = dqc1;
        else
            dqc = dqc2;
        end
        dqc = dqc*sign(qf-qi);
        %constraint for desired acceleration
        v_cons = abs(qf - qi)/ (tm-ti);

        if (abs(dqc) < v_cons) || (abs(dqc) > 2*v_cons)
            disp('Error. Reassign dqc'); 
            disp(i);
            disp('It should be in limits');
            disp(v_cons);
            disp(2*v_cons);
        end
        % time c
        tc = (qi - qf +dqc*tm)/dqc;
        ddqc = dqc^2/(qi - qf +dqc*tm);
    else
        %Assign desired accelerations
        ddqc1 = 3;
        ddqc2 = 3;
        if i==1
            ddqc = ddqc1;
        else
            ddqc = ddqc2;
        end
        %Constrain for accelearation in point c
        ddqc_star = 4*abs(qf - qi)/tm^2;
        ddqc = ddqc*sign(qf-qi);
        if abs(ddqc) < ddqc_star
            disp('Error. Reassign ddqc'); 
            disp(i);
            disp('It should be greater than');
            disp(ddqc_star);
        end
        % time tc
        tc = tm/2 - 1/2*sqrt((tm^2*ddqc-4*(qf-qi))/ddqc);

    end
    tc = round(tc,3);

    %polinomials for trapezoid
    syms t2
    tr1 = qi + 1/2*ddqc*t2^2;
    tr2 = qi + ddqc*tc*(t2-tc/2);
    tr3 = qf - 1/2*ddqc*(tm-t2)^2;

    dtr1 = diff(tr1, t2);
    dtr2 = diff(tr2, t2);
    dtr3 = diff(tr3, t2);

    ddtr1 = diff(dtr1, t2);
    ddtr2 = diff(dtr2, t2);
    ddtr3 = diff(dtr3, t2);

    st=0.001;
    if i==1 
        time1 = ti : st : tc;
    else
        time1 = ti+st : st : tc;
    end
    time2 = (tc+st) : st : (tm - tc);
    time3 = (tm-tc+st) : st : (tm);

    tr1_t = double(subs(tr1, {t2}, {time1}));
    tr2_t = double(subs(tr2, {t2}, {time2}));
    tr3_t = double(subs(tr3, {t2}, {time3}));

    dtr1_t = double(subs(dtr1, {t2}, {time1}));
    dtr2_t = double(subs(dtr2, {t2}, {time2}));
    dtr3_t = double(subs(dtr3, {t2}, {time3}));

    ddtr1_t = double(subs(ddtr1, {t2}, {time1}));
    ddtr2_t = double(subs(ddtr2, {t2}, {time2}));
    ddtr3_t = double(subs(ddtr3, {t2}, {time3}));
    if i==1
        tr_1 = [tr1_t, tr2_t, tr3_t];
        dtr_1 = [dtr1_t, dtr2_t, dtr3_t];
        ddtr_1 = [ddtr1_t, ddtr2_t, ddtr3_t];
        time_tr1 = [time1, time2, time3];
    else
        tr_2 = [tr1_t, tr2_t, tr3_t];
        dtr_2 = [dtr1_t, dtr2_t, dtr3_t];
        ddtr_2 = [ddtr1_t, ddtr2_t, ddtr3_t];
        time_tr2 = [time1, time2, time3];
    end
end

time = ti : st : tf;
figure(4)
tr = [tr_1, tr_2];
plot(time, tr,'b',  'LineWidth', 2)
xlabel('Time of moving')
ylabel('Trajectory')
legend('q_{tr}', 'Fontsize', 10)
title('Trajectory trapezoidal')
grid on

figure(5)
dtr = [dtr_1, dtr_2];

plot(time, dtr, 'g',  'LineWidth', 2)
xlabel('Time of moving')
ylabel('Velocity')
legend('dq_{tr}', 'Fontsize', 10)
title('Velocity trapezoidal')
grid on

figure(6)
ddtr = [ddtr_1, ddtr_2];

plot(time, ddtr, 'r',  'LineWidth', 2)
xlabel('Time of moving')
ylabel('Acceleration')
legend('ddq_{tr}', 'Fontsize', 10)
title('Acceleration trapezoidal')
grid on


