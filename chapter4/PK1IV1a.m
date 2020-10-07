function PK1IV1a

% sequential calculation of drug amounts and concentrations using
% one-compartment PK model with IV input and first-order elimination.

clear all;
clc;
close all;

p.V = 5.0; % L
p.k10 = 1.0; % mg/h
dose = 100; % mg
timeEnd = 10; % h

[t,y] = ode45(@derivatives, [0 timeEnd], dose, [], p);

a = y;
c = a/p.V;
plot(t, a, '-r', t, c, 'b-', 'lineWidth', 2)
legend('Drug amounts [mg]', 'Drug concentrations [mg/l]')
xlabel('Time [h]')
ylabel('Value')
%print('-dtiff', '-r900', 'PK1IV1a.tif')

end

function dydt = derivatives(~, y, p)
% compute the right hand side of the ode

a1 = y(1);
da1 = -p.k10 * a1;
dydt = da1;

end