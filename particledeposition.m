% Particle dry deposition on surface
% 03/06/2020
% Jialei Shen

clc
clear

%% Important parameters
% exposure days
day = 30;
% mass concentration of PM10, ug/m3
Cpm10 = 16;

%% particle size distribution and concentration
% mass concentration for pm10 of sea salt
%Cpm10 = 150; % mass concentration for pm10 of sea salt [ug/m3]

% volume fraction of size distribution of sea salt particle
dx = 0.1;
x = 0.1:dx:10;

for i = 1:length(x)
    if x(i)>=0.1 & x(i)<0.4
        y(i) = -0.515./3.*x(i)+0.9450;
    elseif x(i)>=0.4 & x(i)<1.5
        y(i) = -8.247./11.*x(i)+1.1761;
    else
        y(i) = -0.515./85.*x(i)+5.15./85;
    end
end
Avolume = trapz(x, y);
%{
figure;
semilogx(x,y)
title('Probability density of volume fraction in different particle diameters')
ylabel('Probability density [\mum^{-1}]')
xlabel('Particle diameter [\mum]')
grid on
%}

% mass fraction of sea salt particle
for i = 1:length(x)
    z(i) = (x(i)./2).^3.*y(i);
end

z = y;

Amass = trapz(x, z);
z = z./Amass;
Amass = trapz(x, z);
figure;
semilogx(x,z)
title('Probability density of mass fraction in different particle diameters')
ylabel('Probability density [\mum^{-1}]')
xlabel('Particle diameter [\mum]')
grid on

cm = z.*Cpm10; % mass concentration of sea salt particles in differnt size
figure;
semilogx(x,cm)
title('Probability density of mass concentration in different particle diameters')
ylabel('Probability density of mass concentration [\mug/(m^3\mum)]')
xlabel('Particle diameter [\mum]')
Acm = trapz(x, cm);
grid on

dp = x;
C = cm;
CC = cm.*dx;

figure;
semilogx(x,CC)
%set(gca, 'Xscale', 'log')
title('Mass concentration in different particle diameters')
ylabel('Mass concentration [\mug/m^3]')
xlabel('Particle diameter [\mum]')
grid on
CCtotal = sum(CC);

%% 
%{
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
%}
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
%day = 30;
h = 24*day;
t = 3600*h;
% probability density of cumulative particle mass by particle diameter
md = vd.*A.*C.*t;

%% cumulative particle mass with time
for dd = 1:day
    tt = dd.*24.*3600;
    mdd = vd.*A.*C.*tt; % probability density of deposition mass
    Amdd(dd) = trapz(dp, mdd);
end

figure;
plot(Amdd)
title('Cumulative particle deposition mass on surface with time')
xlabel('Day [day]')
ylabel('Cumulative particle deposition mass [\mug]')
grid on


%% Particle cover area on surface 
ttt = 1; % time step: 1s
mddd = vd.*A.*C.*ttt; % probability density of deposition mass, ug/um
rhoNaCl = 2.165; % NaCl density, g/cm3
rhoNaCl = rhoNaCl .* 1000; % unit conversion, kg/m3
for j = 1:length(mddd)
    vparticle = 4./3.*pi.*((dp(j)./2).*1e-6).^3; % volume of a NaCl particle (ball), m3
    nparticle(j) = mddd(j).*1e-9./(vparticle.*rhoNaCl); % probability density of number of particles on surface, #/um
    aparticle(j) = nparticle(j).*pi.*((dp(j)./2).*1e-6).^2; % probability of particle coverage area, m2/um
end
%figure;
%semilogx(dp,nparticle)

apt = zeros(1,t);
for k = 1:t-1
    if apt(k)<A
        apt(k+1) = apt(k) + trapz(dp,aparticle); % calculate total particle coverage area with time, m2
    else
        apt(k+1) = A;
    end
end

tttt = A./trapz(dp,aparticle); % count how much time it need to fully cover the specimen, s
tttt = tttt./(24.*3600);
fprintf('It will take %.2f days for particles to fully cover the specimen.\n', tttt)

kk = 1:t;
kk = kk./(3600.*24);
figure;
plot(kk,apt.*10000);
title("Total particle coverage area change with time")
xlabel('Time [day]')
ylabel('Total particle coverage area on specimen [cm^2]')
legend("Specimen area is "+num2str(A*10000)+"cm^2; "+num2str(tttt)+" days to fully cover.");
grid on

figure;
plot(kk,100.*apt./A);
title("Total particle coverage percent on specimen")
xlabel('Time [day]')
ylabel('Total particle coverage percent [%]')
legend("Specimen area is "+num2str(A*10000)+"cm^2; "+num2str(tttt)+" days to fully cover.");
grid on


figure;
semilogx(dp,(t.*nparticle.*dx))
title("Cumulative particle numbers in different diameters after "+num2str(h)+" hours")
xlabel('Particle diameter dp [\mum]')
ylabel('Cumulative particle numbers [#]')
grid on

figure;
semilogx(dp,(t.*aparticle.*10000.*dx))
title("Cumulative particle coverage area in different diameters after "+num2str(h)+" hours")
xlabel('Particle diameter dp [\mum]')
ylabel('Cumulative particle coverage area [cm^2]')
grid on

%trapz(dp,(t.*aparticle.*10000))

%{
for j = 1:t
    for pp = 1:length(cm)
        cm(pp)
        x(pp)
    end
end
%}

%% result presentation
figure;
plot(dp, md.*dx)
grid on
set(gca, 'Xscale', 'log')
xlim([0.1,10]);
Amd = trapz(dp, md);
title("Cumulative particle deposition mass after "+num2str(h)+" hours")
xlabel('Particle diameter dp [\mum]')
ylabel('Cumulative particle deposition mass [\mug]')
legend("Total cumulative particle deposition mass = "+num2str(Amd)+"\mug");


fprintf('Total cumulative particle deposition mass after %.0f hours is %.2f ug, and covers %.2f%% area of specimen.\n', h, Amd, (100.*t.*trapz(dp,aparticle)./A));


vd0 = vd_empirical(dp, theta, u0);
figure;
loglog(dp, vd0)
grid on
title('Particle deposition velocity (vd) for different particle diameters')
xlabel('Particle diameter [\mum]')
ylabel('Particle deposition velocity (vd) [m/s]')
legend("U="+num2str(U)+"m/s, u*="+num2str(u0)+"cm/s, theta="+num2str(theta))
