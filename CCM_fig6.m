%% Edge of Chaos of the Chua corsage memristor (CCM)
% Voltage-controlled ideal generic memristor
% It can be derived from a flux-controlled memristor
% State-dependent Ohm's law state equation
% DAE set
% dx/dt = g(x,vM)*vM; state equation
% i = G(x)*v; ohm's law; G(x) = G0.x^2.vM
% G defines the memductance function: G(x) = dq(flux)/dflux
% g defines the morphing function g(x) = dx(flux)/dflux
% Small-signal Analysis and Edge of Caos

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

%% Loci of V vs X
subplot(2,2,1)
plot(x_M,V_M)
xlim([-20 70])
ylabel('V/V')
xlabel('X/Vs')
grid on

% Zero-pole diagrams of Y of the CCM over the voltage range
% -10 < V < -1 where X < 10
LAD (10,"lower",x_M,g_M,V_M,i_M);
% -20 < V < -10 where X > 35
LAD (35,"higher",x_M,g_M,V_M,i_M);


function LAD (xQ,caseQ,x_M,g_M,V_M,i_M)

N = length(x_M);
[val,idx]=min(abs(x_M-xQ));
if caseQ == "lower"
    x_M = x_M(idx+1:N);
    V_M = V_M(idx+1:N);
    g_M = g_M(idx+1:N);
    i_M = i_M(idx+1:N);
else
    x_M = x_M(1:idx+1);
    V_M = V_M(1:idx+1);
    g_M = g_M(1:idx+1);
    i_M = i_M(1:idx+1);
end

a11 = 2*x_M.*V_M;
a12 = x_M.^2;

for i=1:length(x_M)
    if or(x_M(i) < 10, x_M(i) > 35)  
        b11(i) = -1;
    else
        b11(i) = 1;
    end
end

b12 = 1;


Lx = 1./(b12.*a11);
Rx = -b11./(b12.*a11);
Ry = 1./a12;

p = -Rx./Lx;
z = -(Rx+Ry)./Lx;

% LAD1 = -5 < V < -1.6667 --> Edge of Chaos domain 1
% LAD2 = -20 < V < -18.3333 --> Edge of Chaos domain 2

%% Plot of zero vs V
subplot(2,2,2)
hold on
plot(V_M,z, "x")
xlim([-20 -2])
ylim([-30 30])
ylabel('Zero')
xlabel('V/V')
grid on

% Edge of caos domain 2
subplot(2,2,3)
hold on
plot(V_M,z)
xlim([-20 -16])
ylim([-0.3 0.3])
ylabel('Zero')
xlabel('V/V')
grid on

%% Plot of pole vs V
subplot(2,2,4)
hold on
plot(V_M,p)
xlim([-20 -2])
ylim([-2 1])
ylabel('Pole')
xlabel('V/V')
grid on

end

% One-to-one function
function out=x(flux,v_M)%,param)
    
    out = 30*flux - flux.^2/2 + sign(flux-10).*(flux-10).^2/2 - sign(flux-35).*(flux-35).^2/2 + v_M*flux;
    
end

function out=g(x,v_M)%,param) % morphing function

    out = 30 - x + abs(x-10) - abs(x-35) + v_M;
    
end

% function the memductance of memristor
% derivative of constitutive relation
function out=G(x)%,param) % memductance

    G0 = 1;
    out = G0*x.^2;

end


% function of the constitutive relation
function out=q(flux) % charge on memristor

    G0 = 1;
    out = G0*flux.^3/3;

end