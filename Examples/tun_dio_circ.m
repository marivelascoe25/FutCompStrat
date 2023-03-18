% Tunnel diode circuit has solution if there is a parasitic capacitor 
% parallel to diode and inductor
% I_L describes the current through inductor, I_c through capacitor and I
% through diode
% State equations (2nd order system):
% dI_L/dt = +V_c(t)/L
% dV_c/dt = 1/C_p [-I_L - g(V_c)]
% I = g(V) : constitutive relation (CR)

clc,close all,clear all

% fixed parameters
alpha=51; % [A/V] 
beta=24; % [A]
gamma=11; % [V^(-2)] autovalues
L=.1e-3; %inductor
C=10e-3; %capacitor

param=[alpha,beta,gamma,L,C];

%fig, ax = plt.subplots(nrows=1, ncols=3, sharex=True, sharey=True, figsize=(6, 6))
%fig.text(0.5, 0.04, '$v_C$/V', ha='center')
%fig.text(0.04, 0.5, 'common Y', va='center', rotation='vertical')

% initial conditions

 init_v_C=[-0.15 -.175 .2 .15 .425 .225 .575 -.1 -.1 .5];
 init_i_L=[-2 2.5 3 -2.75 -2.75 0 -2.875 0 0.5 -3.0343];
            


% run one simulation per initial condition
for jj=1:length(init_v_C)

% assign one of the aforespecified initial conditions
v_0=init_v_C(jj);
i_0=init_i_L(jj);
y_init=[i_0,v_0];

% simulation time
tspan=[0 50e-3]; 


% numerical integration tolerances
options = odeset('RelTol',1e-6,'AbsTol',1e-9);
%options = odeset('RelTol',1e-1,'AbsTol',1e-4);


% call the ODE solver 
[t,y] = ode15s(@tun_dio_circ_eqs,tspan,y_init,options,param); 
% ode15s, stiff solver (able to take much larger steps)
% options to use odeset function AbsTol and RelTol (tolerances)


% post-processing the numerical solutions
i_L=y(:,1);
v_C=y(:,2);


subplot(1,4,2)
plot(v_C,i_L)
hold on % superposing all 10 plots
xlabel('$v_C$/V','Interpreter','latex')
ylabel('$i_L$/A','Interpreter','latex')
grid on

vv = v_C(1:2:end); %% remove even elements
vvv = vv(1:2:end); 
ii = i_L(1:2:end);
iii = ii(1:2:end); 
subplot(1,4,1)
q = quiver(vvv,iii); %% otherwise too overwelmed
hold on
xlabel('$v_C$/V','Interpreter','latex')
ylabel('$i_L$/A','Interpreter','latex')
grid on


end


subplot(1,4,3)
plot(v_C,g(v_C,param))
xlabel('V(V)','Interpreter','latex')
ylabel('I(A)','Interpreter','latex')
grid on

subplot(1,4,4)
yyaxis left
plot(t,i_L)
ylabel('I(A)')
yyaxis right
plot(t,v_C)
ylabel('V(V)')
xlabel('t(s)')
grid on


% function specifying the ODE under numerical solution
% state equations
function out=tun_dio_circ_eqs(t,y,param)
alpha=param(1);
beta=param(2);
gamma=param(3);
L=param(4);
C=param(5);

i_L=y(1);
v_C=y(2);

out=[1/L*v_C;
     1/C*(-i_L-g(v_C,param))];  % [dI_L/dt ; dV_c/dt]

end

% function specifying the tunnel diode characteristic
% constitutive relation
function out=g(v,param) %current on diode

alpha=param(1);
beta=param(2);
gamma=param(3);
out=alpha*v-beta*(1-exp(-gamma*v.^2));

end