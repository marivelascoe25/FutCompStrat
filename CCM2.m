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
v0 = -20;
vf = 15;
vs = v0:1/fs:vf;
T = 10;
t0 = 0;
dt = (T-t0)/(length(vs)-1);
t = t0:dt:T;

v_M = [-25, -20, -10, 0, 5, 10];

for i = 1:length(v_M)
%flux0 = -26.6060;
flux = cumtrapz(t,vs);%+flux0;
x_M = x(flux,v_M(i));
g_M = g(x_M,v_M(i));

hold on
plot(x_M,g_M,'DisplayName',int2str(v_M(i)))
xlim([-20 70])
ylabel('$g(x,v_{M})$/V','Interpreter','latex')
xlabel('x/Vs')
grid on
legend

end

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