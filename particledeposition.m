% Particle dry deposition on surface
% 03/06/2020
% Jialei Shen

clc
clear

%% particle size distribution and concentration
% particle diameter range
dp = [0.5, 0.7, 1.0, 2.0, 5.0]; %um
% particle concentration at different diameters
Cn = [0.06 0.02 0.03 0.07 0.05; ...
     0.13 0.02 0.03 0.07 0.08; ...
     0.20 0.07 0.05 0.12 0.40; ...
     0.14 0.03 0.03 0.07 0.03; ...
     0.08 0.02 0.01 0.02 0.03; ...
     0.33 0.12 0.07 0.13 0.30; ...
     0.08 0.04 0.03 0.05 0.11
    ]; %ug/m3
C = mean(Cn);
%dpi = [0.5:0.1:5];
%Ci = interp1q(dp',C',dpi');
%dp = dpi;
%C = Ci';

%% specimen 
% inclination angle
theta = 0;
% specimen size
L = 0.1; %m
W = 0.01; %m
A = L.*W; %m2

%% surrounding condition
% friction velocity
% Source: Liu, D.L., 2010. Particle Deposition onto Enclosure Surfaces, in: Developments in Surface Contamination and Cleaning: Particle Deposition, Control and Removal. Elsevier Inc., pp. 1–56. 
% https://doi.org/10.1016/B978-1-4377-7830-4.10001-5
U = 0.5; %free flow velocity, m/s
rhoa = 1.204; %air density at 20degC, kg/m3
nu = 1.516e-5; %kinematic viscosity of air at 20degC, m2/s
dUdy = (0.074./(rhoa.*nu)).*(rhoa.*U.^2./2).*(U.*L./nu).^(-1./5);
u0 = sqrt(nu.*dUdy); %m/s
u0 = 100.*u0 %cm/s
%u0 = 1; %cm/s

%% calculate deposition velocity (vd)
% Empirical function: vd = vd_empirical(dp, theta, u0);
vd = vd_empirical(dp, theta, u0);

%% calculate cumulative deposited particle mass
% cumulative period [s]
h = 24*30;
t = 3600*h;
% calculate cumulative particle mass
md = vd.*A.*C.*t;

%% result presentation
figure;
stem(dp, md)
grid on
set(gca, 'Xscale', 'log')
xlim([0.01,10]);
title("Cumulative particle deposition mass after "+num2str(h)+" hours")
xlabel('Particle diameter dp [\mum]')
ylabel('Cumulative particle deposition mass [\mug]')
legend("Total cumulative particle deposition mass = "+num2str(sum(md))+"\mug");

fprintf('Total cumulative particle deposition mass after %.0f hours is %.2f ug\n', h, sum(md));

figure;
stem(dp,C)
grid on
set(gca, 'Xscale', 'log')
xlim([0.01,10]);
title("Particle concentration in data center")
xlabel('Particle diameter dp [\mum]')
ylabel('Particle concentration [\mug/m^3]')

dp0 = [0.01:0.01:10];
vd0 = vd_empirical(dp0, theta, u0);
figure;
loglog(dp0, vd0)
hold on
stem(dp, vd, 'm')
grid on
title('Particle deposition velocity (vd) for different particle diameters')
xlabel('Particle diameter [\mum]')
ylabel('Particle deposition velocity (vd) [m/s]')
legend("U="+num2str(U)+"m/s, u*="+num2str(u0)+"cm/s, theta="+num2str(theta),'Studied particle size')
