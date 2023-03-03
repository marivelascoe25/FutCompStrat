
clc,close all,clear all

T = 4*pi;
w = 1;
vs = 1.2;
t = 0:0.001:T;
v_M = vs*sin(w*t);

flux0 = -0.4;
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
plot(t,i_M)
xlabel('t/s')
ylabel('i/A')
set(gca,'XTick',0:pi/2:4*pi)
set(gca,'XTickLabel',{'0','\pi/2','\pi','3\pi/2','2\pi','5\pi/2','3\pi','7\pi/2','4\pi'})

subplot(3,2,5)
plot(v_M,i_M)
xlabel('v/V')
ylabel('i/A')

subplot(3,2,4)
yyaxis left
plot(t,flux)
ylabel('$\varphi$/Vs','Interpreter','latex')
yyaxis right
plot(t,q_M)
ylabel('q/C')
xlabel('t/s','Interpreter','latex')

subplot(3,2,6)
yyaxis left
plot(t,G_M)
ylabel('$G(\varphi)$/S','Interpreter','latex')
yyaxis right
plot(t,i_M)
ylabel('i/A')
xlabel('t/s','Interpreter','latex')
set(gca,'XTick',0:pi/2:4*pi)
set(gca,'XTickLabel',{'0','\pi/2','\pi','3\pi/2','2\pi','5\pi/2','3\pi','7\pi/2','4\pi'})

% constitutive relation
function out=q(flux) % charge on memristor

out = -0.00375 + 0.025*flux + 0.04*abs(flux+0.25) - 0.025*abs(flux-0.25);

end

% function the memductance of memristor
% derivative of constitutive relation
function out=G(flux) % memductance

out = 0.025 + 0.04*sign(flux+0.25) - 0.025*sign(flux-0.25);

end