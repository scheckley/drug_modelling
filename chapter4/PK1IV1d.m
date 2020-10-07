function PK1IV1d
% calculation of drug concentrations using Mass option. One-compartment PK
% model with IV input and 1st order elimination. The procedure is based on
% an implcit ODE with the Mass option.

clear all;
close all;
clc;

p.V = 5.0; % L
p.CL10 = 5.0; % L/h
dose = 100; % mg
timeEnd = 10; % h

options = odeset;
options.Mass = p.V;

[t, y] = ode15s(@derivalgeb, [0 timeEnd], dose/p.V, options, p);

c = y;
a = c*p.V;
plot(t, a, '-r', t, c, 'b-', 'lineWidth', 2)
legend('Drug amounts [mg]', 'Drug Concentrations [mg/L]')
xlabel('Time [h]')
ylabel('Value')
end

function dydt = derivalgeb(~, y, p)
%compute the ode's

c1 = y(1);
dc1 = -p.CL10 * c1;
dydt = dc1;
end