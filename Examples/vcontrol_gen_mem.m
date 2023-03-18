%% Voltage-controlled ideal generic memristor
% It can be derived from a flux-controlled memristor
% State-dependent Ohm's law state equation
% DAE set
% dx/dt = g(x)*v; state equation
% i = G(x)*v; ohm's law
% G defines the memductance function: G(x) = dq(flux)/dflux
% g defines the morphing function g(x) = dx(flux)/dflux

clc,close all,clear all

T = 4*pi;
w = 1;
vs = 1.2;
t = 0:0.1:T;
v_M = vs*sin(w*t);
x0 = -1.6;

% simulation time
tspan=[0 T]; 


% numerical integration tolerances
options = odeset('RelTol',1e-6,'AbsTol',1e-9);

% call the ODE solver 
[t,x] = ode15s(@vc_gm,tspan,x0,options,v_M); 
% ode15s, stiff solver (able to take much larger steps)
% options to use odeset function AbsTol and RelTol (tolerances)


%x = cumtrapz(t,g(x)*v_M)+x0
%xaux = x0
%for jj=1:length(v_M)
    %xi = cumtrapz(t,g(xaux)*v_M)+xaux;
    %x(jj) = xi
    %xaux = xi
%end

G_M = G(x);
i_M = G_M.*v_M;


subplot(3,2,1)
plot(x,g(x))
%xlim([-3 3])
xlabel('$\varphi$/Vs','Interpreter','latex')
ylabel('q/C')

subplot(3,2,2)
plot(x,G_M)
%xlim([-3 3])
xlabel('$\varphi$/Vs','Interpreter','latex')
ylabel('$G(\varphi)$/S','Interpreter','latex')

subplot(3,2,3)
yyaxis left
plot(t,v_M)
ylabel('v/V')
yyaxis right
plot(t,x)
ylabel('$\varphi$/Vs','Interpreter','latex')
xlabel('t/s')

subplot(3,2,4)
yyaxis left
plot(t,G_M)
ylabel('$G(\varphi)$/S','Interpreter','latex')
yyaxis right
plot(t,x)
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



% function specifying the ODE under numerical solution
% state equations
function out=vc_gm(t,x,v_M)

out = g(x)*v_M;

%[1/L*v_C;     1/C*(-i_L-g(v_C,param))];  % [dI_L/dt ; dV_c/dt]

end

% constitutive relation
function out=g(x) % morphing function

out = 17/8 - 15/8*sign(17/8*x+15/8*abs(x));


end

% function the memductance of memristor
% derivative of constitutive relation
function out=G(x) % memductance

out = 0.025 + 0.04*sign(17/8*x+15/8*abs(x)+0.25) - 0.025*sign(17/8*x+15/8*abs(x)-0.25);

end