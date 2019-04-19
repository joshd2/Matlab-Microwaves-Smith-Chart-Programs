% Josh Davis
% Microwaves II Stability Calculator
clear,clc,close all
set(0,'DefaultLegendAutoUpdate','off')
% Initializing Input Variables
S11 = [0,0];
S12 = [0,0];
S21 = [0,0];
S22 = [0,0];

%{
% Input Prompts
S11(1,1) = input('Magnitude of S11: ');
S11(1,2) = input('Phase (deg) of S11: ');

S12(1,1) = input('Magnitude of S12: ');
S12(1,2) = input('Phase (deg) of S12: ');

S21(1,1) = input('Magnitude of S21: ');
S21(1,2) = input('Phase (deg) of S21: ');

S22(1,1) = input('Magnitude of S22: ');
S22(1,2) = input('Phase (deg) of S22: ');
%}
% Example input S parameters
%{
% Example Input Data in place of prompts:
S11(1,1) = ;
S11(1,2) = ;

S12(1,1) = ;
S12(1,2) = ;

S21(1,1) = ;
S21(1,2) = ;

S22(1,1) = ;
S22(1,2) = ;
%}
%{
% Example 7b Input Data in place of prompts:
S11(1,1) = 0.385;
S11(1,2) = -55;

S12(1,1) = .045;
S12(1,2) = 90;

S21(1,1) = 2.7;
S21(1,2) = 78;

S22(1,1) = 0.89;
S22(1,2) = -26.5;
%}
%{
% Example 7c Input Data in place of prompts:
S11(1,1) = 0.7;
S11(1,2) = -50;

S12(1,1) = .27;
S12(1,2) = 75;

S21(1,1) = 5;
S21(1,2) = 120;

S22(1,1) = 0.6;
S22(1,2) = 80;
%}

% Example 8 Input Data in place of prompts:
S11(1,1) = 0.69;
S11(1,2) = -78;

S12(1,1) = 0.033;
S12(1,2) = 41.4;

S21(1,1) = 5.67;
S21(1,2) = 123;

S22(1,1) = 0.84;
S22(1,2) = -25;

% Redefining of input variables as complex numbers:
S11 = S11(1,1)*exp(1i*S11(1,2)*pi/180);
S12 = S12(1,1)*exp(1i*S12(1,2)*pi/180);
S21 = S21(1,1)*exp(1i*S21(1,2)*pi/180);
S22 = S22(1,1)*exp(1i*S22(1,2)*pi/180);

% Determinant of S matrix:
Delta = S11*S22-S12*S21;
Delta_Mag = abs(Delta);

% Stability Circles:

% OUTPUT (Load) Stability:
% Center Load Vector for Output Stability:
CL = conj(S22-Delta*conj(S11))/(abs(S22)^2-abs(Delta)^2);
xL = real(CL);
yL = imag(CL);
% Radius of |Gamma_IN| = 1 circle for Output Stability: 
rL = abs(S12*S21/(abs(S22)^2-abs(Delta)^2));

% INPUT (Source) Stability
% Center Source Vector for Input Stability:
CS = conj(S11-Delta*conj(S22))/(abs(S11)^2-abs(Delta)^2);
xS = real(CS);
yS = imag(CS);
% Radius of |Gamma_OUT| = 1 circle for Input Stability: 
rS = abs(S12*S21/(abs(S11)^2-abs(Delta)^2));

% Unconditional Stability Tests:
K = (1-abs(S11)^2-abs(S22)^2+abs(Delta)^2)/abs(2*S12*S21);
if K > 1 && Delta_Mag < 1
    disp('Unconditionally Stable by K and |\Delta|')
else
    disp('Potentially Unstable by K and |\Delta|')
end

% Edwards-Sinsky Criterion:
mu1 = (1-abs(S11)^2)/(abs(S22-Delta*conj(S11))+abs(S12*S21));
mu2 = (1-abs(S22)^2)/(abs(S11-Delta*conj(S22))+abs(S12*S21));
if mu1 > 1 || mu2 > 1
    disp('Unconditionally Stable by \mu')
else
    disp('Potentially Unstable by \mu')
end

% Displaying Information
fprintf('\nInfo:\n');
fprintf('C_L Mag = %4.4d\n',abs(CL));
fprintf('C_L Phase = %4.4d degrees\n',angle(CL));
fprintf('r_L = %4.4d\n',rL);
fprintf('C_S Mag = %4.4d\n',abs(CS));
fprintf('C_S Phase = %4.4d degrees\n',angle(CS));
fprintf('r_S = %4.4d\n',rS);
fprintf('K = %4.4d\n',K);
fprintf('Delta = %4.4d\n',Delta);
fprintf('mu_1 = %4.4d\n',mu1);
fprintf('mu_2 = %4.4d\n',mu2);



% Circle Plotting
% Input Stability
figure(1)
hold on
% Plot the unit circle (aka smith chart)
th = 0:pi/5000:2*pi;
plot(cos(th),sin(th),'k')
axis equal
axis([-1.1 2 -1.1 1.1])
title('Output Stability Circles')
% Plot figure range:
range = -1.5:0.001:1.5;
[x1, y1] = meshgrid(range);
condL1 = x1.^2 + y1.^2 < 1; % Condition for within unit circle
if abs(S22) < 1 && rL < abs(CL)
    % Condition for stability circle. S22 < 1 means center smith 
    % chart is stable circle, so whether that's contained within
    % stability circle or not depends on if the stability circle's
    % radius is less than the C center vector (implying outside 
    % stability circle) or not.
    condL2 = (x1 - xL).^2 + (y1 - yL).^2 > rL^2;
else
    condL2 = (x1 - xL).^2 + (y1 - yL).^2 < rL^2;
end
condL1 = double(condL1);  % convert to double for plotting
condL2 = double(condL2);
condL1(condL1 == 0) = NaN;  % set the 0s to NaN so they are not plotted
condL2(condL2 == 0) = NaN;
condL = condL1.*condL2;  % multiply the two condaces to keep only the common points
s1 = surf(x1,y1,condL,'FaceAlpha',0.25);
s1.EdgeColor = 'none';
s1.FaceColor = [0.3010 0.7450 0.9330];
view(0,90)
% Plot Load Output Limit
xunitL = rL * cos(th) + xL;
yunitL = rL * sin(th) + yL;
plot(xunitL, yunitL,'--r');
legend('Smith Chart', 'Stable Region','|\Gamma_{IN} = 1| circle')
% Plot horizontal line and smith chart curves
line([-1,1],[0,0],'Color','black') % line([xmin:-1,xmax:1],[ymin:0,ymax:0])
% R circles: 
% center: (Gamma_r,Gamma_i) = (R/(1+R),0)
% radius: 1/(1+R)
for R = [.25 .5 .75 1 2 3 4 5]
    xfunction = 1/(1+R)*cos(th) + R/(1+R);
    yfunction = 1/(1+R)*sin(th) + 0;
    plot(xfunction,yfunction,'r')
end
% X circles:
% center: (Gamma_r,Gamma_i) = (1,1/X)
% radius: 1/abs(X)
for X = [-5 -4 -3 -2 -1 -.75 -.5 -.25 0 .25 .5 .75 1 2 3 4 5]
    xfunction = 1/abs(X)*cos(th) + 1;
    yfunction = 1/abs(X)*sin(th) + 1/X;
    % Eliminate points that are outside the unit circle
    xfunction(xfunction.^2+yfunction.^2 > 1) = NaN;
    yfunction(xfunction.^2+yfunction.^2 > 1) = NaN;
    plot(xfunction,yfunction,'r')
end
% G circles: 
% center: (Gamma_r,Gamma_i) = (-G/(1+G),0)
% radius: 1/(1+G)
for G = [.25 .5 .75 1 2 3 4 5]
    xfunction = 1/(1+G)*cos(th) - G/(1+G);
    yfunction = 1/(1+G)*sin(th) + 0;
    plot(xfunction,yfunction,'g')
end
% B circles:
% center: (Gamma_r,Gamma_i) = (-1,1/B)
% radius: 1/abs(B)
for B = [-5 -4 -3 -2 -1 -.75 -.5 -.25 0 .25 .5 .75 1 2 3 4 5]
    xfunction = 1/abs(B)*cos(th) - 1;
    yfunction = 1/abs(B)*sin(th) + 1/B;
    % Eliminate points that are outside the unit circle
    xfunction(xfunction.^2+yfunction.^2 > 1) = NaN;
    yfunction(xfunction.^2+yfunction.^2 > 1) = NaN;
    plot(xfunction,yfunction,'g')
end
hold off



figure(2)
hold on
% Plot the unit circle (aka smith chart)
th = 0:pi/5000:2*pi;
plot(cos(th),sin(th),'k')
axis equal
axis([-1.1 2 -1.1 1.1])
title('Input Stability Circles')
% Plot figure range:
range = -1.5:0.001:1.5;
[x2, y2] = meshgrid(range);
condS1 = x2.^2 + y2.^2 < 1; % Condition for within unit circle
if abs(S11) < 1 && rS < abs(CS)
    % Condition for stability circle. S11 < 1 means center smith 
    % chart is stable circle, so whether that's contained within
    % stability circle or not depends on if the stability circle's
    % radius is less than the C center vector (implying outside 
    % stability circle) or not.
    condS2 = (x2 - xS).^2 + (y2 - yS).^2 > rS^2;
else
    condS2 = (x2 - xS).^2 + (y2 - yS).^2 < rS^2;
end
condS1 = double(condS1);  % convert to double for plotting
condS2 = double(condS2);
condS1(condS1 == 0) = NaN;  % set the 0s to NaN so they are not plotted
condS2(condS2 == 0) = NaN;
condS = condS1.*condS2;  % multiply the two condaces to keep only the common points
s2 = surf(x2,y2,condS,'FaceAlpha',0.25);
s2.EdgeColor = 'none';
s2.FaceColor = [0.3010 0.7450 0.9330];
view(0,90)
% Plot Source Input Limit
xunitS = rS * cos(th) + xS;
yunitS = rS * sin(th) + yS;
plot(xunitS, yunitS,'--r');
legend('Smith Chart', 'Stable Region','|\Gamma_{OUT} = 1| circle')
% Plot horizontal line and smith chart curves
line([-1,1],[0,0],'Color','black') % line([xmin:-1,xmax:1],[ymin:0,ymax:0])
% R circles: 
% center: (Gamma_r,Gamma_i) = (R/(1+R),0)
% radius: 1/(1+R)
for R = [.25 .5 .75 1 2 3 4 5]
    xfunction = 1/(1+R)*cos(th) + R/(1+R);
    yfunction = 1/(1+R)*sin(th) + 0;
    plot(xfunction,yfunction,'r')
end
% X circles:
% center: (Gamma_r,Gamma_i) = (1,1/X)
% radius: 1/abs(X)
for X = [-5 -4 -3 -2 -1 -.75 -.5 -.25 0 .25 .5 .75 1 2 3 4 5]
    xfunction = 1/abs(X)*cos(th) + 1;
    yfunction = 1/abs(X)*sin(th) + 1/X;
    % Eliminate points that are outside the unit circle
    xfunction(xfunction.^2+yfunction.^2 > 1) = NaN;
    yfunction(xfunction.^2+yfunction.^2 > 1) = NaN;
    plot(xfunction,yfunction,'r')
end
% G circles: 
% center: (Gamma_r,Gamma_i) = (-G/(1+G),0)
% radius: 1/(1+G)
for G = [.25 .5 .75 1 2 3 4 5]
    xfunction = 1/(1+G)*cos(th) - G/(1+G);
    yfunction = 1/(1+G)*sin(th) + 0;
    plot(xfunction,yfunction,'g')
end
% B circles:
% center: (Gamma_r,Gamma_i) = (-1,1/B)
% radius: 1/abs(B)
for B = [-5 -4 -3 -2 -1 -.75 -.5 -.25 0 .25 .5 .75 1 2 3 4 5]
    xfunction = 1/abs(B)*cos(th) - 1;
    yfunction = 1/abs(B)*sin(th) + 1/B;
    % Eliminate points that are outside the unit circle
    xfunction(xfunction.^2+yfunction.^2 > 1) = NaN;
    yfunction(xfunction.^2+yfunction.^2 > 1) = NaN;
    plot(xfunction,yfunction,'g')
end
hold off

% Old Scripts:
%{
% Circle Plotting
figure(1)
hold on
% Plot the unit circle (aka smith chart)
th = 0:pi/50:2*pi;
plot(cos(th),sin(th),'k')
axis square
% Plot Load Output Limit
xL = real(CL);
yL = imag(CL);
xunitL = rL * cos(th) + xL;
yunitL = rL * sin(th) + yL;
plot(xunitL, yunitL,'--r');
title('Output Stability Circles')
if abs(S22) < 1
    legend('Smith Chart', 'Unstable within the circle')
else
    legend('Smith Chart', 'Unstable outside of circle')
end
hold off


figure(2)
hold on
% Plot the unit circle (aka smith chart)
plot(cos(th),sin(th),'k')
axis square
% Plot Source Input Limit
xS = real(CS);
yS = imag(CS);
th = 0:pi/50:2*pi;
xunitS = rS * cos(th) + xS;
yunitS = rS * sin(th) + yS;
plot(xunitS, yunitS,'--r');
title('Input Stability Circles')
if abs(S11) < 1
    legend('Smith Chart', 'Unstable within the circle')
else
    legend('Smith Chart', 'Unstable outside of circle')
end
hold off
%}
