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

fs = 1000; %Hz
v0 = -25;
vf = 10;
vs = v0:1/fs:vf;
T = 10;
t0 = 0;
dt = (T-t0)/(length(vs)-1);
t = t0:dt:T;

v_M = -25;
N = length(v_M);
%x_M = [];
%flux0 = -26.6060;
flux = cumtrapz(t,vs);%+flux0;

%for i = 1:N
x_M = x(flux,v_M);
g_M = g(x_M,v_M);
V_M = -g(x_M,0);
i_M = G(x_M).*V_M;

subplot(3,3,1)
hold on
plot(x_M,g_M)
xlim([-20 70])
ylabel('$g(x,v_{M})$/V','Interpreter','latex')
xlabel('x/Vs')
grid on

subplot(3,3,2)
hold on
plot(x_M,V_M)
xlim([-20 70])
ylabel('V/V')
xlabel('x/Vs')
grid on

subplot(3,3,4)
hold on
plot(x_M,i_M)
xlim([-20 70])
ylabel('I/A')
xlabel('x/Vs')
grid on

subplot(3,3,5)
hold on
plot(x_M,i_M)
xlim([-2 5])
ylabel('I/A')
xlabel('x/Vs')
grid on

subplot(3,3,6)
hold on
plot(x_M,i_M)
xlim([34 38])
ylim([-24700 -24400])
ylabel('I/A')
xlabel('x/Vs')
grid on

subplot(3,3,7)
hold on
plot(V_M,i_M)
xlim([-30 25])
ylim([-25000 23000])
ylabel('I/A')
xlabel('V/V')
grid on

subplot(3,3,8)
hold on
plot(V_M,i_M)
xlim([-6 -1])
ylim([-20 0])
ylabel('I/A')
xlabel('V/V')
grid on

subplot(3,3,9)
hold on
plot(V_M,i_M)
xlim([-20 -18])
ylim([-24650 -24500])
ylabel('I/A')
xlabel('V/V')
grid on

% One-to-one function
function out=x(flux,v_M)
    
    out = 30*flux - flux.^2/2 + sign(flux-10).*(flux-10).^2/2 - sign(flux-35).*(flux-35).^2/2 + v_M*flux;
    %out = (-4000 + 180*flux - 3*flux.^2 + (5600 - 260*flux + 3*flux.^2)*sign(40 - flux) - (-20 + flux)*sign(20 - flux)*(-40 + 3*flux + (-40 + flux)*sign(40 - flux)))/8;
end

function out=g(x,v_M) % morphing function

    out = 30 - x + abs(x-10) - abs(x-35) + v_M;%*ones(length(x));

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