% Flux-controlled ideal memristor
% State-dependent Ohm's law state equation
% i = dq/dt = dq(flux)/d(flux) * d(flux)/dt = G(flux)*v
% d(flux)/dt = v
% G defines the memory-conductance or memductance at Q
% i and v are in memristor
% vs is supplied voltage = v
% q = ^q(flux): constitutive relation (CR)


clc,close all,clear all

n = input('Source 1) Sine; 2) Pulse: ');

switch n
    case 1
        T = 30;
        w = 1;
        vs = 1;
        t = 0:0.1:T;
        v_M = vs*sin(w*t);
    case 2
        t0 = 1;
        dt = 1;
        tf = 4;
        E = 1; % [V] Amplitud of pulse 
        fs = 10e3; %ensure precision, if higher, slope is visible: triangle
        t = 0:1/fs:tf;
        v_M = rectpuls(t-t0-dt/2,dt);
    otherwise
        disp('No valid entry')
end

flux0 = 0;
flux = cumtrapz(t,v_M)+flux0;
G_M = G(flux);
i_M = G_M.*v_M;
q_M = q(flux);

subplot(2,2,1)
yyaxis left
plot(t,v_M)
ylabel('v/V')
yyaxis right
plot(t,i_M)
ylabel('i/A')
xlabel('t/s')

subplot(2,2,2)
plot(v_M,i_M)
xlabel('v/V')
ylabel('i/A')

subplot(2,2,3)
yyaxis left
plot(t,G_M)
ylabel('$G(\varphi)$/S','Interpreter','latex')
yyaxis right
plot(t,flux)
hold off
ylabel('$\varphi$/Vs','Interpreter','latex')
xlabel('t/s')

subplot(2,2,4)
plot(flux,q_M)
%xlim([-3 3])
xlabel('$\varphi$/Vs','Interpreter','latex')
ylabel('q/C')

% function of the constitutive relation
function out=q(flux) % charge of memristor

    out = flux + 1/3*flux.^3;

end

% function the memductance of memristor
% derivative of constitutive relation
function out=G(flux) % memductance: dq(flux)/dflux

    out = 1 + flux.^2;

end