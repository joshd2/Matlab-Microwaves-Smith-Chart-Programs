% Josh Davis
% Microwaves II Plot Smith Chart
clear,clc,close all
% Circle Plotting
% Input Stability
figure(1)
hold on
% Plot the unit circle (aka smith chart)
th = 0:pi/5000:2*pi;
plot(cos(th),sin(th),'k')
axis square
line([-1,1],[0,0],'Color','black') % line([xmin:-1,xmax:1],[ymin:0,ymax:0])
% R circles: 
% center: (Gamma_r,Gamma_i) = (R/(1+R),0)
% radius: 1/(1+R)
for R = [.25 .5 1 2 5]
    xfunction = 1/(1+R)*cos(th) + R/(1+R);
    yfunction = 1/(1+R)*sin(th) + 0;
    plot(xfunction,yfunction,'r')
end
% X circles:
% center: (Gamma_r,Gamma_i) = (1,1/X)
% radius: 1/abs(X)
for X = [-5 -2 -1 -.5 -.25 0 .25 .5 1 2 5]
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
for G = [.25 .5 1 2 5]
    xfunction = 1/(1+G)*cos(th) - G/(1+G);
    yfunction = 1/(1+G)*sin(th) + 0;
    plot(xfunction,yfunction,'g')
end
% B circles:
% center: (Gamma_r,Gamma_i) = (-1,1/B)
% radius: 1/abs(B)
for B = [-5 -2 -1 -.5 -.25 0 .25 .5 1 2 5]
    xfunction = 1/abs(B)*cos(th) - 1;
    yfunction = 1/abs(B)*sin(th) + 1/B;
    % Eliminate points that are outside the unit circle
    xfunction(xfunction.^2+yfunction.^2 > 1) = NaN;
    yfunction(xfunction.^2+yfunction.^2 > 1) = NaN;
    plot(xfunction,yfunction,'g')
end
axis([-1.1 1.1 -1.1 1.1])
hold off