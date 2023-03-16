%% Chua corsage memristor (CCM)
%% Voltage-controlled ideal generic memristor
% It can be derived from a flux-controlled memristor
% State-dependent Ohm's law state equation
% DAE set
% dx/dt = g(x)*v; state equation
% i = G(x)*v; ohm's law
% G defines the memductance function: G(x) = dq(flux)/dflux
% g defines the morphing function g(x) = dx(flux)/dflux

clc,close all,clear all

n = input('Source 1) Sine; 2) Pulse: ');

switch n
    case 1
        T = 4*pi;
        w = 1;
        vs = 1.2;
        t = 0:0.01:T;
        v_M = vs*sin(w*t);

    case 2
        t0 = 1;
        dt = 1;
        tf = 4;
        E = 1; % [V] Amplitud of pulse 
        fs = 10e3; %ensure precision, if higher, slope is visible: triangle
        t = 0:1/fs:tf;
        v_M = rectpuls(t-t0-dt/2,dt);

    case 3        
        fs = 100; %%Hz
        v0 = -10;
        vf = -3.334;
        v_M = v0:1/fs:vf;

    otherwise
        disp('No valid entry')
end
syms xx
%G_M = diff(q(xx),xx);
%g_M = diff(x(xx),xx);

flux0 = -0.4;
flux = cumtrapz(t,v_M)+flux0;
x_M = x(flux);

G_M = G(x_M);
g_M = g(x_M);
i_M = G_M.*v_M;

%d = diff(f);
%val = subs(d,x,pi/2);
%G_M = diff(q(flux),flux);
%i_M = G_M(x_M).*v_M;
%g_M = diff(x_M,flux);
%q_M = q(flux);

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
%yyaxis left
plot(t,G_M)
ylabel('$G(x)$/S','Interpreter','latex')
%yyaxis right
%plot(t,flux)
%hold off
%ylabel('$\varphi$/Vs','Interpreter','latex')
xlabel('t/s')

subplot(2,2,4)
plot(t,x_M)
%xlim([-3 3])
xlabel('$\varphi$/Vs','Interpreter','latex')
ylabel('q/C')

% One-to-one function
function out=x(flux)
    
    out = 17/8*flux - 15/8*abs(flux);

end

function out=g(x) % morphing function

out = 30 - x + abs(x-20) - abs(x-40);


end

% function the memductance of memristor
% derivative of constitutive relation
function out=G(x) % memductance

G0 = 1;
out = G0*x.^2;

end


% function of the constitutive relation
function out=q(flux) % charge on memristor

out = -0.00375 + 0.025*flux + 0.04*abs(flux+0.25) - 0.025*abs(flux-0.25);

end

% function the memductance of memristor
% derivative of constitutive relation
function out=Gf(flux) % memductance

out = 0.025 + 0.04*sign(flux+0.25) - 0.025*sign(flux-0.25);

end