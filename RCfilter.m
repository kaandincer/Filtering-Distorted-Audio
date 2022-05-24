% function cascadedRCfilter(Vin,h) receives a time-series voltage sequence
% sampled with interval h, and returns the output voltage sequence produced
% by a circuit
%
% inputs:
% Vin - time-series vector representing the voltage input to a circuit
% h - scalar representing the sampling interval of the time series in
% seconds
%
% outputs:
% Vout - time-series vector representing the output voltage of a circuit

function Vout = RCfilter(Vin,h)
R1 = 510;
R4 = 16;

A = [1, -1, -1, 0, 0, 0;
     R1, 0, 0, -1, 1, 0;
     0, 0, R4, 0, 0 , -1;
     0, 0, 0, 1, 0, 0;
     0, 0, 0, 0, 1, 0;
     0, 0, 0, 0, 1, -1];

k=length(Vin);
C_2=.47e-6;
C_3=.47e-2;
Vc2 = zeros(1,k);
Vc3 = zeros(1,k);
Vout = zeros(1,k);

for i = 1:k
    b = [0, 0, 0, Vin(i), Vc2(i), Vc3(i)]';
    x = linsolve(A, b);
    Vc2(i+1) = Vc2(i)+(h/C_2)*x(1);
    Vc3(i+1) = Vc3(i)+(h/C_3)*x(3);
    Vout(i) = x(6);
end

end