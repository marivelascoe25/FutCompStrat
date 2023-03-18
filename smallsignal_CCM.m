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

fs = 500; %Hz
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

%param = [a0, a1, b2, c21, c22, c23, c24, c25, d0, d1, d2, d3, d4];
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

%LAD1 = -5 < V < -1.6667
%LAD2 = -20 < V < -18.3333

Lx = 1./(b12.*a11);
Rx = -b11./(b12.*a11);
Ry = 1./a12;

p = -Rx./Lx;
z = -(Rx+Ry)./Lx;

w0 = -10;
wf = 10;
dw = (wf-w0)/(length(Lx)-1);
w = -10:dw:10;
Y_Re = Re(w,Rx,Ry,Lx);
Y_Im = Im(w,Rx,Ry,Lx);

subplot(2,3,1)
plot(x_M,g_M)
xlim([-20 70])
ylabel('$g(x,v_{M})$/V','Interpreter','latex')
xlabel('x/Vs')
grid on

subplot(2,3,2)
plot(V_M,z)
xlim([-20 -2])
ylim([-30 30])
ylabel('Zero')
xlabel('V/V')
grid on

subplot(2,3,3)
plot(V_M,p)
xlim([-20 -2])
ylim([-2 1])
ylabel('Pole')
xlabel('V/V')
grid on

subplot(2,3,4)
yyaxis left
plot(w,Y_Im)
ylabel('ImY/S')
yyaxis right
plot(w,Y_Re)
ylabel('ReY/S')
xlabel('w/(rad/s)')
grid on

% subplot(2,3,5)
% plot(V_M,p)
% xlim([-20 -2])
% ylim([-2 1])
% ylabel('Pole')
% xlabel('V/V')
% grid on

%Local admittance real part
function out=Re(w,Rx,Ry,Lx)

    num = Rx.*Ry.*(Rx+Ry) + w.^2.*Lx.^2.*Ry;
    den = Rx.^2.*Ry.^2 + w.^2.*Lx.^2.*Ry.^2;
    out = num./den;

end

%Local admittance imaginry part
function out=Im(w,Rx,Ry,Lx)

    out = -w.*Lx./(Rx.^2 + w.^2.*Lx.^2);

end

% One-to-one function
function out=x(flux,v_M)%,param)
    
    out = 30*flux - flux.^2/2 + sign(flux-10).*(flux-10).^2/2 - sign(flux-35).*(flux-35).^2/2 + v_M*flux;
    %out = param(1)*flux + param(2)*flux.^2/2 + param(3)*v_M.^2 + param(4)*v_M.^2*flux.^2/2 + param(4)*v_M.^2*flux.^3/3 + param(5)*v_M.^2*flux.^4/4 + param(6)*v_M.^2*flux.^5/5 + param(7)*v_M.^2*flux.^6/6;

end

function out=g(x,v_M)%,param) % morphing function

    out = 30 - x + abs(x-10) - abs(x-35) + v_M;
    %out = param(1) + param(2)*x + param(3)*v_M.^2 + param(4)*v_M.^2*x + param(5)*v_M.^2*x.^2 + param(6)*v_M.^2*x.^3 + param(7)*v_M.^2*x.^4 + param(8)*v_M.^2*x.^5;

end

% function the memductance of memristor
% derivative of constitutive relation
function out=G(x,param) % memductance

    G0 = 1;
    out = G0*x.^2;
    %out = param(9) + param(10)*x + param(11)*x.^2 + param(12)*x.^3 + param(13)*x.^4;

end


% function of the constitutive relation
function out=q(flux) % charge on memristor

    G0 = 1;
    out = G0*flux.^3/3;

end