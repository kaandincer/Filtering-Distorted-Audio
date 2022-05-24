%% Part 1: Model an RC Circuit

%Plot using own h
h = 6e-6;
k = 1000;
B = h/(10^3*1e-6);
A = 1-B;
V_in = 5*ones([k+1,1]);
system1 = ss(A, B, [], [], []);
[~,~,V_c] = lsim(system1, V_in, 0:k);

figure();
hold on;
plot((0:k)*h, V_in);
plot((0:k)*h, V_c);
xlabel('Time (s)');
ylabel('Voltage (V)');
legend('V_i_n', 'V_c', 'location', 'east');
title('Voltage of Capacitor over time with Sampling Interval h=6e-6s');

%Plot using h'
hNew = 8e-4;
k = 4;
B = hNew/(10^3*1e-6);
A = 1-B;
V_in = 5*ones([k+1,1]);
system1 = ss(A, B, [], [], []);
[~,~,V_c] = lsim(system1, V_in, 0:k);

figure();
hold on;
plot((0:k)*h, V_in);
plot((0:k)*h, V_c);
xlabel('Time (s)');
ylabel('Voltage (V)');
legend('V_i_n', 'V_c', 'location', 'east');
title('Voltage of Capacitor over time with Sampling Interval h=8e-4s');

%Plot of theoretical curve
R = 1000;
C = 1e-6;
k = 1000;
t = (0:k)*h;
theoretical = 5*(1-exp(-t/(R*C)));
V_in = 5*ones([k+1,1]);

figure();
hold on;
plot((0:k)*h, V_in);
plot((0:k)*h, theoretical);
xlabel('Time (s)');
ylabel('Voltage (V)');
legend('V_i_n','V_c', 'location', 'east');
title('Theoretical Voltage of Capacitor Over Time');

%% Part 2: Compare two RC circuits

% At 50Hz (lower frequency)
h = 6e-6;
k = 5000;
t = (0:k)*h;
f = 50;
B = h/(10^3*1e-6);
A = 1-B;
csystem1 = ss(A, B, [], [], []);
[~,~,V_c] = lsim(system1, V_in, 0:k);

V_r = V_in - V_c';

figure();
hold on;
plot((0:k)*h, V_in);
plot((0:k)*h, V_c);
plot((0:k)*h, V_r);
xlabel('Time (s)');
ylabel('Voltage (V)');
legend('V_i_n', 'V_c', 'V_R', 'location', 'northeast');
title('Voltage over time at f = 50 HZ');

%At 1000Hz (higher frequency)
h = 6e-6;
k = 500;
t = (0:k)*h;
f = 1000;
B = h/(10^3*1e-6);
A = 1-B;
V_in = 5*sin(2*pi*f*t);
system1 = ss(A, B, [], [], []);
[~,~,V_c] = lsim(system1, V_in, 0:k);

V_r = V_in - V_c';

figure();
hold on;
plot((0:k)*h, V_in);
plot((0:k)*h, V_c);
plot((0:k)*h, V_r);
xlabel('Time (s)');
ylabel('Voltage (V)');
legend('V_i_n', 'V_c', 'V_R', 'location', 'northeast');
title('Voltage over time at f = 1000 HZ');

%Plotting transfer function
transferC = zeros([10000-10,1]);
transferR = zeros([10000-10,1]);

h = 6e-6;
k = 500;
B = h/(10^3*1e-6);
A = 1-B;

for f = 10:10000
    
   V_in = 5*sin(2*pi*f*t);    
   
   system1 = ss(A, B, [], [], []);
   [~,~,V_c] = lsim(system1, V_in, 0:k);
   
   V_r = V_in - V_c';
   transferC(f-9) = max(V_c)/max(V_in);    
   transferR(f-9) = max(V_r)/max(V_in);
    
end

figure();
hold on;
plot(log(10:10000), log(transferC));
plot(log(10:10000), log(transferR));
xlabel('Log of Frequency (Hz)');
ylabel('Log of Transfer Function');
legend('H1', 'H2', 'location', 'east');
title('Transfer Functions over Inupt Frequencies');

%% Part 3: Compare two cascaded RC circuits

% For Circuit "C"
% create an A matrix
R_1 = 330;
R_2 = 330;
R_4 = 330;
C_1 = 0.68e-6;
C_2 = 0.68e-6;
C_3 = 0.68e-6;
A = [1, -1, -1, 0, 0, 0;
     0, 330, 0, 0, -1, 0;
     0, 0, 330, 0, 0 , -1;
     0, 0, 0, 1, -1, 0;
     0, 0, 0, 0, 1, -1;
     0, 0, 0, 1, 0, 0];
 
f1 = 440;
f2 = 3000;
h = 1e-6;
k = 10000;
t = (0:k-1)*h;
Vin = 5*sin(2*pi*f1*t)+sin(2*pi*f2*t);

Vc1 = zeros(1,k);
Vc3 = zeros(1,k);
V_out = zeros(1,k);

% simulate the system
for i = 1:k
    b = [0, 0, 0, Vc1(i), Vc3(i), Vin(i)]';
    x = linsolve(A, b);
    V_out(i) = x(6);
    
    %use some values from x in the equations below
    Vc1(i+1) = Vc1(i)+(h/C_1)*x(1);
    Vc3(i+1) = Vc3(i)+(h/C_3)*x(3);
    % update b with the new Vc1 and Vc2 values
    
end

% plot Vout
figure();
hold on;
plot(t, Vin);
plot(t, V_out);
xlabel('Time (s)');
ylabel('Voltage (V)');
legend('V_i_n', 'V_o_u_t', 'location', 'northeast');
title('V_i_n and V_o_u_t over time for Circuit "C"');
hold off;



% For Circuit "D"

newA = [1, -1, -1, 0, 0, 0;
     330, 0, 0, -1, 1, 0;
     0, 0, 330, 0, 0 , -1;
     0, 0, 0, 1, 0, 0;
     0, 0, 0, 0, 1, 0;
     0, 0, 0, 0, 1, -1];

f1 = 440;
f2 = 3000;
h = 1e-6;
k = 5000;
t = (0:k-1)*h;
Vin = 5*sin(2*pi*f1*t)+sin(2*pi*f2*t);

Vc2 = zeros(1,k);
Vc3 = zeros(1,k);
V_out = zeros(1,k);

for i = 1:k
    b = [0, 0, 0, Vin(i), Vc2(i), Vc3(i)]';
    x = linsolve(newA, b);
    Vc2(i+1) = Vc2(i)+(h/C_2)*x(2);
    Vc3(i+1) = Vc3(i)+(h/C_3)*x(3);
    V_out(i) = x(6);
end

figure();
hold on;
plot(t, Vin);
plot(t, V_out);
xlabel('Time (s)');
ylabel('Voltage (V)');
legend('V_i_n', 'V_o_u_t', 'location', 'northeast');
title('V_i_n and V_o_u_t over time for Circuit "D"');
hold off;


%plotting transfer functions
transferC = zeros([10000-10,1]);
transferD = zeros([10000-10,1]);

h = 1e-6;
k = 1000;
t = (0:k-1)*h;

V_outC = zeros(1,k);
V_outD = zeros(1,k);

Vc1_C = zeros(1,k);
Vc2_C = zeros(1,k);
Vc3_C = zeros(1,k);
Vc1_D = zeros(1,k);
Vc2_D = zeros(1,k);
Vc3_D = zeros(1,k);

for f = 10:10000
    
    Vin = 5*sin(2*pi*f*t)+sin(2*pi*f*t);
    
    for i = 1:k
        
        % Circuit C
        b_C = [0, 0, 0, Vc1_C(i), Vc3_C(i), Vin(i)]';
        x_C = linsolve(A, b_C);
        V_outC(i) = x_C(6);
    
        Vc1_C(i+1) = Vc1_C(i)+(h/C_1)*x_C(1);
        Vc3_C(i+1) = Vc3_C(i)+(h/C_3)*x_C(3);    
    
        % Circuit D
        b_D = [0, 0, 0, Vin(i), Vc2_D(i), Vc3_D(i)]';
        x_D = linsolve(newA, b_D);
        Vc2_D(i+1) = Vc2_D(i)+(h/C_2)*x_D(2);
        Vc3_D(i+1) = Vc3_D(i)+(h/C_3)*x_D(3);
        V_outD(i) = x_D(6);
    end

    transferC(f-9) = max(V_outC)/max(Vin);
    transferD(f-9) = max(V_outD)/max(Vin);
    
end

figure();
hold on;
plot(log(10:10000), log(transferC));
plot(log(10:10000), log(transferD));
xlabel('Log of Frequency (Hz)');
ylabel('Log of Transfer Function');
legend('H1', 'H2', 'location', 'east');
title('Transfer Functions over Inupt Frequencies');
hold off;
