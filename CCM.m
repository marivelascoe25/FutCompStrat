%% Chua corsage memristor (CCM)
% Voltage-controlled ideal generic memristor
% It can be derived from a flux-controlled memristor
% State-dependent Ohm's law state equation
% DAE set
% dx/dt = g(x,vM)*vM; state equation
% i = G(x)*v; ohm's law; G(x) = G0.x^2.vM
% G defines the memductance function: G(x) = dq(flux)/dflux
% g defines the morphing function g(x) = dx(flux)/dflux

clc,close all,clear all

fs = 100; %Hz
v0 = -10;
vf = 0;
vs = v0:1/fs:vf;
T = 90;
t0 = 0;
dt = (T-t0)/(length(vs)-1);
t = t0:dt:T;

v_M = -25;%[-25, -20, -10, 0, 5, 10]

%flux0 = -26.6060;
flux = cumtrapz(t,vs);%+flux0;
x_M = x(flux,v_M);

x0=v_M+5;
% simulation time
%tspan=[0 90];
% numerical integration tolerances
options = odeset('RelTol',1e-6,'AbsTol',1e-9);
%options = odeset('RelTol',1e-1,'AbsTol',1e-4);
[tt,x_M2] = ode15s(@gg,t,x0,options,v_M);
x_M2 = x_M2';

G_M = G(x_M);
G_M2 = G(x_M2);
g_M = g(x_M,v_M);
g_M2 = g(x_M2,v_M);
V_M = -g(x_M,v_M);
I = G_M.*V_M;
I2 = G_M2.*vs;
%i_M = G_M.*v_M;


subplot(3,2,1)
plot(x_M,g_M)
xlabel('x_M/Vs')
ylabel('g_M/V')
%yyaxis left
%plot(t,vs)
%ylabel('v/V')
%yyaxis right
%plot(t,I)
%ylabel('i/A')
%xlabel('t/s')

subplot(3,2,2)
plot(x_M2,g_M2)
xlabel('x_M2/Vs')
ylabel('g_M2/V')
%plot(vs,I)
%xlabel('v/V')
%ylabel('i/A')

subplot(3,2,3)
plot(vs,I)
%xlim([-10,0])
xlabel('v/V')
ylabel('i/A')
%box on;
%plot(vs,I)
%xlim([-10,0])
%grid on;

subplot(3,2,4)
plot(vs,I2)
xlim([-10,0])
xlabel('v/V')
ylabel('i2/A')

subplot(3,2,5)
plot(t,x_M)
xlabel('t/s')
ylabel('x_M1/Vs')

subplot(3,2,6)
plot(tt,x_M2)
xlabel('t/s')
ylabel('x_M2/Vs')

% One-to-one function
function out=x(flux,v_M)
    
    %out = 17/8*flux - 15/8*abs(flux);
    out = 30*flux - flux.^2/2 + sign(flux-20).*(flux-20).^2/2 - sign(flux-40).*(flux-40).^2/2 + v_M*flux;
    %out = (-4000 + 180*flux - 3*flux.^2 + (5600 - 260*flux + 3*flux.^2)*sign(40 - flux) - (-20 + flux)*sign(20 - flux)*(-40 + 3*flux + (-40 + flux)*sign(40 - flux)))/8;
end

function out=g(x,v_M) % morphing function

out = 30 - x + abs(x-20) - abs(x-40) + v_M;


end

function out=gg(t,x,v_M) % morphing function

out = 30 - x + abs(x-20) - abs(x-40) + v_M;

end

% function the memductance of memristor
% derivative of constitutive relation
function out=G(x) % memductance

G0 = 1;
out = G0*x.^2;

end


% function of the constitutive relation
function out=q(flux) % charge on memristor

G0 = 1;
out = G0*flux.^3/3;

end