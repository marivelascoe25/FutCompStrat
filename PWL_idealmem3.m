
clc,close all,clear all

t0 = 1;
dt = 1.6;
tf = 4.5;
E = -1.25; % [V] Amplitud of pulse 
fs = 10e3; %ensure precision, if higher, slope is visible: triangle
t = 0:1/fs:tf;
v_M = E*rectpuls(t-t0-dt/2,dt);

flux0 = 2.5;
flux = cumtrapz(t,v_M)+flux0;
G_M = G(flux);
i_M = G_M.*v_M;
q_M = q(flux);

subplot(3,2,1)
plot(flux,q_M)
%xlim([-3 3])
xlabel('$\varphi$/Vs','Interpreter','latex')
ylabel('q/C')

subplot(3,2,2)
plot(flux,G_M)
%xlim([-3 3])
xlabel('$\varphi$/Vs','Interpreter','latex')
ylabel('$G(\varphi)$/S','Interpreter','latex')

subplot(3,2,3)
yyaxis left
plot(t,v_M)
ylabel('v/V')
yyaxis right
plot(t,flux)
ylabel('$\varphi$/Vs','Interpreter','latex')
xlabel('t/s')

subplot(3,2,4)
yyaxis left
plot(t,G_M)
ylabel('$G(\varphi)$/S','Interpreter','latex')
yyaxis right
plot(t,flux)
ylabel('$\varphi$/Vs','Interpreter','latex')
xlabel('t/s')

subplot(3,2,5)
plot(v_M,i_M)
xlabel('v/V')
ylabel('i/A')

subplot(3,2,6)
yyaxis left
plot(t,G_M)
ylabel('$G(\varphi)$/S','Interpreter','latex')
yyaxis right
plot(t,i_M)
ylabel('i/A')
xlabel('t/s')

% constitutive relation
function out=q(flux) % charge on memristor

out = 1 + 1.5*flux - 0.5*abs(flux-2);

end

% function the memductance of memristor
% derivative of constitutive relation
function out=G(flux) % memductance

out = 1.5 - 0.5*sign(flux-2);

end