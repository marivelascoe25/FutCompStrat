%% Current-controlled ideal generic memristor
% It can be derived from a flux-controlled memristor
% State-dependent Ohm's law state equation
% DAE set
% dx/dt = g(x)*v; state equation
% i = G(x)*v; ohm's law
% G defines the memductance function: G(x) = dq(flux)/dflux
% g defines the morphing function g(x) = dx(flux)/dflux

clc,close all,clear all

fs = 100; %Hz
T = 2*pi;
w = 1;
is = 1;
t = 0:1/fs:T;
i_M = is*sin(w*t);

q0 = -1;
q = cumtrapz(t,i_M)+q0;
x_M = x(q,i_M);

R_M = R(x_M,i_M);
f_M = f(x_M,i_M);
v_M = R_M.*i_M;

subplot(2,2,1)
yyaxis left
plot(t,v_M)
ylabel('v/V')
yyaxis right
plot(t,i_M)
ylabel('i/A')
xlabel('t/s')

subplot(2,2,2)
plot(i_M,v_M)
xlabel('i/A')
ylabel('v/V')

subplot(2,2,3)
yyaxis left
plot(t,R_M)
ylabel('$R(x)$/Ohm','Interpreter','latex')
yyaxis right
plot(t,f_M)
%hold off
ylabel('$f(x)$/A','Interpreter','latex')
xlabel('t/s')

subplot(2,2,4)
plot(q,v_M)
%xlim([-3 3])
xlabel('q/C')
ylabel('V/V')

% One-to-one function
function out=x(q,i)
    
    out =i.*q;

end

% constitutive relation
function out=f(x,i) % morphing function

%out = 17/8 - 15/8*sign(17/8*x+15/8*abs(x));
out = i;

end

% function the memductance of memristor
% derivative of constitutive relation
function out=R(x,i) % memductance

%out = 0.025 + 0.04*sign(17/8*x+15/8*abs(x)+0.25) - 0.025*sign(17/8*x+15/8*abs(x)-0.25);
out = x./i;

end