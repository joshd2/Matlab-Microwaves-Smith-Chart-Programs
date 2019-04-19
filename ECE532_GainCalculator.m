% Josh Davis
% Microwaves II Gain Equations
clear,clc,close all
% Initializing Input Variables
Gamma_S = [0,0];
Gamma_L = [0,0];
S11 = [0,0];
S12 = [0,0];
S21 = [0,0];
S22 = [0,0];

% Input Prompts
Gamma_S(1,1) = input('Magnitude of Gamma_S: ');
Gamma_S(1,2) = input('Phase (deg) of Gamma_S: ');

Gamma_L(1,1) = input('Magnitude of Gamma_L: ');
Gamma_L(1,2) = input('Phase (deg) of Gamma_L: ');

S11(1,1) = input('Magnitude of S11: ');
S11(1,2) = input('Phase (deg) of S11: ');

S12(1,1) = input('Magnitude of S12: ');
S12(1,2) = input('Phase (deg) of S12: ');

S21(1,1) = input('Magnitude of S21: ');
S21(1,2) = input('Phase (deg) of S21: ');

S22(1,1) = input('Magnitude of S22: ');
S22(1,2) = input('Phase (deg) of S22: ');

% Redefining of input variables as complex numbers:
Gamma_S = Gamma_S(1,1)*exp(1i*Gamma_S(1,2)*pi/180);
Gamma_L = Gamma_L(1,1)*exp(1i*Gamma_L(1,2)*pi/180);
S11 = S11(1,1)*exp(1i*S11(1,2)*pi/180);
S12 = S12(1,1)*exp(1i*S12(1,2)*pi/180);
S21 = S21(1,1)*exp(1i*S21(1,2)*pi/180);
S22 = S22(1,1)*exp(1i*S22(1,2)*pi/180);

% Calculating Results
Gamma_IN = S11 + S12*S21*Gamma_L/(1-S22*Gamma_L);
Gamma_OUT = S22 + S12*S21*Gamma_S/(1-S11*Gamma_S);

% Transducer Gain
GT = (1-abs(Gamma_S)^2)/abs(1-Gamma_IN*Gamma_S)^2*abs(S21)^2*(1-abs(Gamma_L)^2)/abs(1-S22*Gamma_L)^2;
% Operating Gain
GP = 1/abs(1-abs(Gamma_IN)^2)*abs(S21)^2*(1-abs(Gamma_L)^2)/abs(1-S22*Gamma_L)^2;
% Available Gain
GA = (1-abs(Gamma_S)^2)/abs(1-S11*Gamma_S)^2*abs(S21)^2*1/(1-abs(Gamma_OUT)^2);

% Source Mismatch Factor:
MS = (1-abs(Gamma_S)^2)*(1-abs(Gamma_IN)^2)/abs(1-Gamma_S*Gamma_IN)^2;
% Load Mismatch Factor:
ML = (1-abs(Gamma_L)^2)*(1-abs(Gamma_OUT)^2)/abs(1-Gamma_L*Gamma_OUT)^2;



